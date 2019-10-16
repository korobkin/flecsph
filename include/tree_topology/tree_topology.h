/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  //
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_topology_tree_topology_h
#define flecsi_topology_tree_topology_h

/*!
  \file tree_topology.h
  \authors nickm@lanl.gov
  \date Initial file creation: Apr 5, 2016
 */

/*
  Tree topology is a statically configured N-dimensional hashed tree for
  representing localized entities, e.g. particles. It stores entities in a
  configurable branch type. Inserting entities into a branch can cause that
  branch to be refined or coarsened correspondingly. A client of tree topology
  defines a policy which defines its branch and entity types and other
  compile-time parameters. Specializations can define a policy and default
  branch types which can then be specialized in a simpler fashion
  (see the basic_tree specialization).
*/

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <float.h>
#include <functional>
#include <iostream>
#include <map>
#include <math.h>
#include <mpi.h>
#include <mutex>
#include <omp.h>
#include <set>
#include <stack>
#include <thread>
#include <unordered_map>
#include <vector>

#include "flecsi/data/data_client.h"
#include "flecsi/geometry/point.h"

#include "tree_geometry.h"
#include "tree_types.h"
#include "hashtable.h"


namespace flecsi {
namespace topology {

/*!
  The tree topology is parameterized on a policy P which defines its branch and
  entity types.
 */
template <class P> class tree_topology : public P, public data::data_client_t {

public:
  using Policy = P; // Tree policy defined by the user

  static const size_t dimension = Policy::dimension; // Current dimension: 1,2,3
  using element_t = typename Policy::element_t; // Type of element either F or D
  using point_t = point_u<element_t, dimension>;
  using range_t = std::array<point_t, 2>;
  using key_t = typename Policy::key_t;
  using branch_id_t = key_t;
  //using branch_t = typename Policy::branch_t;
  using entity_t = typename Policy::entity_t;
  using tree_entity_t = tree_entity<dimension,element_t,key_t,entity_t>; 
  using geometry_t = tree_geometry<element_t, dimension>;
  
  using cofm_t = cofm_u<dimension,element_t,key_t>; 
  using hcell_t = hcell<dimension,key_t,cofm_t,entity_t>;  

  // Hasher for the branch id used in the unordered_map data structure
  template <class KEY> struct branch_id_hasher__ {
    size_t operator()(const KEY &k) const noexcept {
      return k.value() & ((1 << 22) - 1);
    }
  };

  /*!
    Constuct a tree topology with unit coordinates, i.e. each coordinate
    dimension is in range [0, 1].
   */
  tree_topology() {
    //max_depth_ = 0;
    //ghosts_entities_.resize(max_traversal);
    //current_ghosts = 0;
  }

  /**
   * @brief Destroy the tree: empty the hash-table and destroy the entities
   * lists
   */
  ~tree_topology() {
  }

  /**
   * Clean the tree topology but not the local bodies
   */
  void clean() {
    cofm_.clear(); 
    htable_.clear();
    // Reset the root in the table 
    //htable_.emplace(key_t::root(),key_t::root());
    //root_ = htable_.find(key_t::root()); 
  }



  /**
   * Reset the ghosts local information for the next tree traversal
   */
  void reset_ghosts(bool do_share_edge = true) {}

  void share_edge() {}

  /**
   * \brief Change the range of the tree topology
   */
  void set_range(const range_t &range) { range_ = range; }

  /**
   * @brief Get the range
   */
  const std::array<point_t, 2> &range() { return range_; }

  hcell_t* child(hcell_t* h, size_t ci){
    // Use the hash table
    key_t bid = h->key();         // Branch id
    bid.push(ci);                       // Add child number
    auto child = htable_.find(bid); // Search for the child
    // If it does not exists, return nullptr
    if (child == htable_.end())
      return nullptr;
    return &htable_->second;
  }

  /**
   * @ brief Return a reference to the vector of the entities
   */
  std::vector<entity_t> &entities() { return entities_; }

  /**
   * @brief Return an entity by its id
   */
  template <typename E> entity_t &entity(E e) {
    return entities_[static_cast<int>(e)];
  }

  /**
   * @brief Find all the center of mass of the tree up to the
   * maximum sub particles criterion.
   *
   * @param b The starting branch for the search, usually root
   * @param criterion The maximum of subparticles for the COMs
   * @param search_list The extracted COMs
   */
  void find_sub_cells(hcell_t *b, const uint64_t &criterion,
                      std::vector<hcell_t *> &search_list) {}

  /**
` * @brief Apply a function ef to the sub_cells using asynchronous comms.
  * @param [in] b The starting branch for the traversal
  * @param [in] ncritical The number of subparticles for the branches
  * @param [in] b The starting branch for the traversal
  * @param [in] b The starting branch for the traversal
  * @return <return_description>
  * @details <details>
  */
  template <typename EF, typename... ARGS>
  void traversal_sph(EF &&ef, ARGS &&... args) {
    // Perform a tree traversal applying the specified function 
    // on the neighbors of the entities 
    entities_w_ = entities_; 
    
    // Loop for all the entities \TODO change to group them
    #pragma omp parallel for   
    for(int i = 0 ; i < entities_w_.size(); ++i){
      entity_t& ent = entities_w_[i]; 
      point_t center = ent.coordinates(); 
      element_t radius = ent.radius(); 
      std::vector<entity_t*> neighbors;
      neighbors.reserve(60); 
      key_t nkey;  

      std::stack<hcell_t*> stk; 
      stk.push(root()); 
      hcell_t* daughters[nchildren_]; 
            
      while(!stk.empty()){
        hcell_t* cur = stk.top(); 
        stk.pop();
        int children = 0;
        nkey = cur->key(); 
        for(int j = 0 ; j < nchildren_; ++j){
          if(cur->get_child(j)){
            key_t ckey = nkey; ckey.push(j); 
            auto it = htable_.find(ckey); 
            daughters[children++] = &(htable_.find(ckey)->second); 
          } // if 
        } // for 
        // Loop on the children and remove the non-used ones 
        for(int j = 0 ; j < children; ++j){
        //for(int j = children-1; j >= 0; --j){
          if(daughters[j]->is_node()){
            element_t dist2 = 0.; 
            point_t d = cofm_[daughters[j]->node_idx()].coordinates(); 
            dist2 = (d[0]-center[0])*(d[0]-center[0]);
            if constexpr (dimension == 2){
              dist2 += (d[1]-center[1])*(d[1]-center[1]);
            }else if constexpr (dimension == 3){
              dist2 += (d[1]-center[1])*(d[1]-center[1]);
              dist2 += (d[2]-center[2])*(d[2]-center[2]);
            }
            element_t extent = std::max(radius,
              cofm_[daughters[j]->node_idx()].lap())+
              cofm_[daughters[j]->node_idx()].radius(); 
            if(dist2 <= extent*extent){
              stk.push(daughters[j]); 
            } // if
          }else{
            element_t dist2 = 0.; 
            point_t d = entities_[daughters[j]->entity_idx()].coordinates(); 
            dist2 = (d[0]-center[0])*(d[0]-center[0]);
            if constexpr (dimension == 2){
              dist2 += (d[1]-center[1])*(d[1]-center[1]);
            }else if constexpr (dimension == 3){
              dist2 += (d[1]-center[1])*(d[1]-center[1]);
              dist2 += (d[2]-center[2])*(d[2]-center[2]);
            }
            element_t extent = std::max(radius, 
              entities_[daughters[j]->entity_idx()].radius()); 
            if(dist2 <= extent*extent){
              neighbors.push_back(&entities_[daughters[j]->entity_idx()]); 
            } // if
          } 
        } // for 
      } // while
      ef(ent,neighbors,std::forward<ARGS>(args)...); 
    } // for  
    entities_ = entities_w_; 
  } // traversal_sph

  /**
   * @brief Perform a tree traversal in parallel using omp threads for
   * each working branch. If a branch is not local, make a request to the
   * thread handler.
   * @param [in] do_square True is varible smoothing length false otherwize
   * @param [in] work_branch The set of branch on which to compute the inter.
   * @param [in] non_local_branches The branches identified as non local
   * during this traversal, they will be computed after this traversal when the
   * distant particles will be gathered
   * @param [in] assert_local In the case of local run assert that no distant
   * particles are reached.
   * @param [in] ef The function to apply to the particles in the work_branch
   * @param [in] args Arguments of the function ef
   * @return void
   * @details
   */
  template <typename EF, typename... ARGS>
  void traverse_sph(std::vector<hcell_t *> &working_branches,
                    std::vector<hcell_t *> &non_local_branches,
                    const bool assert_local, EF &&ef, ARGS &&... args) {}   // traverse_sph

  /**
   * @brief Compute the interaction list for all the sub-particles in b
   * @param [in] b Branch on which to propagate the force
   * @param [in] inter_list The interaction list of the sub-particles of b
   * @param [in] ef The function to apply on the particles
   * @param [in] args The arguments of the function to apply
   * @return void
   * @details This function apply a tree traversal because the branch b can
   * be different than a leaf.
   */
  bool interactions_branches(hcell_t *work_branch,
                             std::vector<hcell_t *> &inter_list,
                             std::vector<hcell_t *> &non_local) {
  }

  void interactions_particles(hcell_t *working_branch,
                              const std::vector<hcell_t *> &inter_list,
                              std::vector<std::vector<entity_t *>> &neighbors) {
  }

#if 0 
  template <typename FC, typename DFCDR, typename DFCDRDR, typename C2P>
  void traversal_fmm(hcell_t *b, double maxmasscell, const double MAC,
                     FC &&f_fc, DFCDR &&f_dfcdr, DFCDRDR &&f_dfcdrdr,
                     C2P &&f_c2p) {

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    clog_one(trace) << "FMM : maxmasscell: " << maxmasscell << " MAC: " << MAC
                    << std::endl;

    entities_w_ = entities_;

    std::vector<branch_t *> work_branch;
    std::vector<branch_t *> remaining_branches;
    find_sub_cells(b, 32, work_branch);

    if (size != 1) {
      std::thread handler(&tree_topology::handle_requests, this);
      // Start tree traversal
      traverse_fmm(work_branch, remaining_branches, MAC, f_fc, f_dfcdr,
                    f_dfcdrdr, f_c2p);
      MPI_Send(NULL, 0, MPI_INT, rank, MPI_DONE, MPI_COMM_WORLD);
      // Wait for communication thread
      handler.join();
    } else {
      traverse_fmm(work_branch, remaining_branches, MAC, f_fc, f_dfcdr,
                    f_dfcdrdr, f_c2p);
    }

    // Check if no message remainig
#ifdef DEBUG
    int flag = 0;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
               &status);
    if(flag != 0){
      std::cerr<<rank<<" TAG:"<<status.MPI_TAG<<" SOURCE: "<<
        status.MPI_SOURCE<<std::endl;
    }
    assert(flag == 0);
#endif
    for (size_t i = 0; i < ghosts_entities_[current_ghosts].size(); ++i) {
      entity_t &g = ghosts_entities_[current_ghosts][i];
      auto id = make_entity(g.key(), g.coordinates(), nullptr, g.owner(),
                            g.mass(), g.id(), g.radius());
      insert(id);
      auto nbi = get(id);
      nbi->set_entity_ptr(&g);
      // Set the parent to local for the search
      find_parent(g.key()).set_ghosts_local(true);
    }
    if(ghosts_entities_[current_ghosts].size() > 0){
      ++current_ghosts;
      assert(current_ghosts < max_traversal);
   
      // compute the particle to particle interaction on each sub-branches
      if (size != 1) {
        std::vector<branch_t *> ignore;
        traverse_fmm(remaining_branches, ignore, MAC, f_fc, f_dfcdr, f_dfcdrdr,
                      f_c2p);
      } else {
        assert(remaining_branches.size() == 0);
      }
    }
    entities_ = entities_w_;


    // Synchronize the threads to be sure they dont start to share edges
    // Critical
    MPI_Barrier(MPI_COMM_WORLD);

  } // apply_sub_cells

  /**
   * @brief Similar to traverse_sph but using the MAC criterion and applying
   * the C2C and C2P computations
   */
  template <typename FC, typename DFCDR, typename DFCDRDR, typename C2P>
  void traverse_fmm(std::vector<branch_t *> &work_branch,
                     std::vector<branch_t *> &remaining_branches,
                     const double MAC, FC &&f_fc, DFCDR &&f_dfcdr,
                     DFCDRDR &&f_dfcdrdr, C2P &&f_c2p) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t total_done = 0;

    size_t nelem = work_branch.size();

#pragma omp parallel for
    for (size_t i = 0; i < nelem; ++i) {
      std::vector<branch_t *> inter_list;
      std::vector<branch_t *> requests_branches;

      total_done += work_branch[i]->sub_entities();
      // Sub traversal and compute the MAC

      if (traversal_c2c_c2p_p2p(work_branch[i], requests_branches, MAC, f_fc,
                                f_dfcdr, f_dfcdrdr, f_c2p)) {
        if(size == 1 ) assert(false);
        std::vector<key_t> send;
#pragma omp critical
        {
          remaining_branches.push_back(work_branch[i]);
          // Send branch key to request handler
          for (auto b : requests_branches) {
            assert(b->owner() < size && b->owner() >= 0);
            assert(b->owner() != rank);
            if (!b->requested()) {
              send.push_back(b->key());
              b->set_requested(true);
            }
          }
        }
        MPI_Request request;
        MPI_Send(&(send[0]), send.size() * sizeof(key_t), MPI_BYTE, rank,
                 LOCAL_REQUEST, MPI_COMM_WORLD);
      }
    }
  }

  /**
   * @brief Center of Mass (COM) to COM interaction and the COM to the particles
   * in this COM.
   */
  template <typename FC, typename DFCDR, typename DFCDRDR, typename F_C2P>
  bool traversal_c2c_c2p_p2p(branch_t *b, std::vector<branch_t *> &non_local,
                             const double MAC, FC &&f_fc, DFCDR &&f_dfcdr,
                             DFCDRDR &&f_dfcdrdr, F_C2P &&f_c2p) {
    point_t fc = 0.;
    double dfcdr[9] = {0.};
    double dfcdrdr[27] = {0.};
    point_t coordinates = b->coordinates();

    std::vector<branch_t *> queue;
    std::vector<branch_t *> new_queue;
    std::vector<point_t> c2c_coordinates;
    std::vector<double> c2c_mass;

    // -------- 1. C2C interactions, keep track of the leaves non interacted
    std::vector<branch_t *> interactions_leaves;
    queue.push_back(root());
    while (!queue.empty()) {
      new_queue.clear();
      const int queue_size = queue.size();
      for (int i = 0; i < queue_size; ++i) {
        branch_t *q = queue[i];
        for (int d = 0; d < (1 << dimension); ++d) {
          if (!q->as_child(d))
            continue;
          branch_t *c = child(q, d);
          if (geometry_t::box_MAC(coordinates, c->coordinates(), c->bmin(),
                                  c->bmax(), MAC)) {
            c2c_coordinates.push_back(c->coordinates());
            c2c_mass.push_back(c->mass());
          } else if (!c->is_leaf()) {
            new_queue.push_back(c);
          } else if (!c->is_local() && !c->ghosts_local()) {
            non_local.push_back(c);
          } else {
            interactions_leaves.push_back(c);
          }
        }
      }
      queue.clear();
      queue = new_queue;
    }

    // Compute the C2C
    for (int i = 0; i < c2c_coordinates.size(); ++i) {
      // Compute the matrices
      f_fc(fc, coordinates, c2c_coordinates[i], c2c_mass[i]);
      f_dfcdr(dfcdr, coordinates, c2c_coordinates[i], c2c_mass[i]);
      f_dfcdrdr(dfcdrdr, coordinates, c2c_coordinates[i], c2c_mass[i]);
    }

    // If all the sub particles are present
    if (non_local.size() == 0) {
      // Propagate this information to the sub-particles for C2P
      for (int i = b->begin_tree_entities(); i <= b->end_tree_entities(); ++i) {
        if (tree_entities_[i].is_local()) {
          f_c2p(fc, dfcdr, dfcdrdr, b->coordinates(), &(entities_w_[i]));
        }
      }
      // Apply the P2P for the local particles
      for (int l = 0; l < interactions_leaves.size(); ++l) {
        for (auto k : *(interactions_leaves[l])) {
          for (int i = b->begin_tree_entities(); i <= b->end_tree_entities();
               ++i) {
            if (tree_entities_[k].id() != tree_entities_[i].id()) {
              // N square computation
              point_t fc;
              entities_w_[i].setAcceleration(
                  entities_w_[i].getAcceleration() +
                  f_fc(fc, entities_w_[i].coordinates(),
                       tree_entities_[k].coordinates(),
                       tree_entities_[k].mass()));
            }
          }
        }
      }
    }
    return non_local.size() > 0;
  }
#endif 
  void find_level(hcell_t *start, const int &level,
                  std::vector<hcell_t *> &find) {}

  /*!
    Return an index space containing all entities within the specified
    spheroid.
   */
  template <typename EF>
  std::vector<entity_t*> find_in_radius(const point_t &center, element_t radius,
                                    EF &&ef) {}

  /*!
      Return an index space containing all entities within the specified
      Box
     */
  template <typename EF>
  std::vector<entity_t*> find_in_box(const point_t &min, const point_t &max,
                                 EF &&ef) {}

  void get_leaves(std::vector<hcell_t *> &leaves) {}

  void remove_non_local() {}

  /**
   * Insert directly a branch (certainly remote) in the tree
   */
  void insert_branch(const point_t &coordinates, const element_t &mass,
                     const point_t &bmin, const point_t &bmax, const key_t &key,
                     const int &owner, const size_t &sub_entities) {
#if 0 
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(rank != owner);
    // Check if this key already exists
    auto itr = branch_map_.find(key);
    // Case 1, branch does not exists localy
    if (itr == branch_map_.end()) {
      // Add the missing parents
      key_t pk = key;
      int last_bit = pk.last_value();
      pk.pop();
      while (branch_map_.find(pk) == branch_map_.end()) {
        branch_map_.emplace(pk, pk);
        itr = branch_map_.find(pk);
        itr->second.set_ghosts_local(false);
        itr->second.set_coordinates(coordinates);
        itr->second.set_mass(mass);
        itr->second.set_bmin(bmin);
        itr->second.set_bmax(bmax);
        itr->second.set_owner(owner);
        itr->second.set_sub_entities(sub_entities);
        itr->second.set_locality(branch_t::NONLOCAL);
        itr->second.set_leaf(false);
        itr->second.add_bit_child(last_bit);
        last_bit = pk.last_value();
        pk.pop();
      }
      // Set upper level not to leave
      branch_map_.find(pk)->second.set_leaf(false);
      branch_map_.find(pk)->second.add_bit_child(last_bit);

      branch_map_.emplace(key, key);
      itr = branch_map_.find(key);
      itr->second.set_ghosts_local(false);
      itr->second.set_coordinates(coordinates);
      itr->second.set_mass(mass);
      itr->second.set_bmin(bmin);
      itr->second.set_bmax(bmax);
      itr->second.set_owner(owner);
      itr->second.set_sub_entities(sub_entities);
      itr->second.set_locality(branch_t::NONLOCAL);
      itr->second.set_leaf(true);

      size_t depth = key.depth();
      // Set the new depth of the tree
      max_depth_ = std::max(max_depth_,depth);
    } else {
      if (itr->second.owner() == rank)
        assert(itr->second.is_shared());
      else {
        // DO NOTHING, this branch have already been updated
        assert(!itr->second.is_local());
      }
    }
    // Add this branch if does not exists
#endif 
  }


  /**
   * @brief Compute the keys of all the entities present in the structure
   */
  void compute_keys() {
    //key_t::set_range(range_);
#pragma omp parallel for
    for (size_t i = 0; i < entities_.size(); ++i) {
      entities_[i].set_key(key_t(range_,entities_[i].coordinates()));
    }
  }

  /*!
    Return the tree's current max depth.
   */
  size_t max_depth() const { return max_depth_; }

  /*!
    Get the root branch (depth 0).
   */
  hcell_t *root() { return &root_->second; }

  /**
   * @brief Generic information for the tree topology
   */
  friend std::ostream &operator<<(std::ostream &os, tree_topology &t) {
    auto r = t.htable_.find(key_t::root());
    cofm_t* root_ptr = &(t.cofm_[r->second.node_idx()]); 
    os << "Tree: "
       << "#node: " << t.htable_.size()-t.entities_.size();
    os << " depth: " << t.max_depth_;
    os << " #root_subents: " << root_ptr->sub_entities();
    os << " center: "<<root_ptr->coordinates(); 
    os << " mass: "<< root_ptr->mass(); 
    os << " radius: "<< root_ptr->radius(); 
    return os;
  }

  size_t pos = 0; 

  /**
   * Loop over the bodies to insert them in the tree and construct the 
   * branches 
   **/
  void build_tree(){


    int size, rank; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size); 

    /* Exchange high and low bound */
    key_t hikey = entities_[entities_.size()-1].key(); 
    key_t lokey = entities_[0].key(); 

    key_t hibound, lobound; 
    exchange_boundaries_(hikey,lokey,hibound,lobound); 

    max_depth_ = 0; 
    // Add the root 
    htable_.emplace(key_t::root(),key_t::root());
    root_ = htable_.find(key_t::root()); 

    size_t nnodes = 0; 
    size_t current_depth = key_t::max_depth(); 
    // Entity keys, last and current 
    key_t lastekey = key_t(0); 
    if(rank != 0) lastekey = lobound; 
    key_t ekey; 
    // Node keys, last and Current 
    key_t lastnkey = key_t::root(); 
    key_t nkey, loboundnode, hiboundnode; 
    // Current parent and value 
    hcell_t * parent = nullptr;
    int oldidx = -1; 

    bool iam0 = rank==0; 
    bool iamlast = rank==size-1; 

    assert(lobound <= lokey);
    assert(hibound >= hikey); 

    for(int i = 0; i <= entities_.size(); ++i){
      if(i < entities_.size()){
        ekey = entities_[i].key(); 
        // Compute the current node key 
      }else{
        ekey = hibound; 
      }
      nkey = ekey; nkey.pop(current_depth);
      bool loopagain = false;
      // Loop while there is a difference in the current keys 
      while(nkey != lastnkey || i == entities_.size()){
        loboundnode = lobound; loboundnode.pop(current_depth); 
        hiboundnode = hibound; hiboundnode.pop(current_depth); 
        if(loopagain && 
        (iam0 || lastnkey > loboundnode) && 
        (iamlast || lastnkey < hiboundnode)){
          // This node is done, we can compute CoFM
          finish_(lastnkey); 
        }
        if(lastnkey == key_t::root()) break;
        loopagain = true; 
        current_depth++;
        nkey = ekey; nkey.pop(current_depth);  
        lastnkey = lastekey; lastnkey.pop(current_depth); 
        
      }

      if(i == entities_.size()) break;

      parent = &(htable_.find(lastnkey)->second); 
      oldidx = parent->entity_idx(); 
      // Insert the eventual missing parents in the tree 
      // Find the current parent of the two entities 
      while(1){
        current_depth--; 
        lastnkey = lastekey; lastnkey.pop(current_depth);
        nkey = ekey; nkey.pop(current_depth);
        if(nkey != lastnkey) break; 
        // Add a children 
        int bit = nkey.last_value(); 
        parent->add_child(bit);
        parent->set_entity_idx(-1); 
        htable_.emplace(nkey,nkey);
        ++nnodes; 
        parent = &(htable_.find(nkey)->second); 
      }

      // Recover deleted entity 
      if(oldidx != -1){
        int bit = lastnkey.last_value(); 
        parent->add_child(bit); 
        parent->set_entity_idx(-1);
        ++nnodes; 
        htable_.emplace(lastnkey,hcell_t(lastnkey,i-1));  
      }

      // Insert new entity
      int bit = nkey.last_value(); 
      parent->add_child(bit); 
      htable_.emplace(nkey,hcell_t(nkey,i)); 
      
      // Prepare next loop  
      lastekey = ekey; 
      lastnkey = nkey; 
      max_depth_ = std::max(max_depth_,current_depth); 
    }

    share_nodes_(); 
  }

  /**
  * Complete the CoFM in the tree with new entities and branches
  */
  void cofm(hcell_t *current) {
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    key_t nkey = current->key(); 
    // Just do something for the nodes, if ptr not set yet
    if(current->is_unset()){
      std::vector<hcell_t*> daughters; 
      daughters.reserve(nchildren_);
      for(int i = 0 ; i < nchildren_; ++i){
        if(current->get_child(i)){
          key_t ckey = nkey; ckey.push(i); 
          auto it = htable_.find(ckey); 
          assert(it != htable_.end()); 
          daughters.push_back(&(htable_.find(ckey)->second));
        } // if
      } // for
      // Loop over daughters first 
      for(int i = 0 ; i < daughters.size(); ++i){
        cofm(daughters[i]); 
      }
      // Add new cofm 
      current->set_shared(); 
      current->set_node_idx(cofm_.size());  
      cofm_.push_back(nkey);
      cofm_children_(&cofm_[current->node_idx()],daughters);
    } // if
  }

  /**
   * Traverse the tree and draw by levels in tikz  
   **/
  void tikz_draw(const char* prefix, int level){
    // Create the file 
    char filename[64];
    sprintf(filename,"%s_%05d.tex",prefix,level);
    std::ofstream output; 
    output.open(filename); 

    // Create the file header 
    output<<"\\documentclass{standalone}"<<std::endl;
    output<<"\\usepackage{tikz}"<<std::endl;
    output<<"\\begin{document}"<<std::endl; 
    output<<"\\begin{tikzpicture}"<<std::endl;

    // Output the tree 
    std::stack<hcell_t*> stk;
    stk.push(root());

    std::vector<hcell_t*> queue; 
    std::vector<hcell_t*> nqueue; 
    queue.push_back(root()); 
    int clevel = -1; 
    while(!queue.empty()){
      for(hcell_t* e: queue){
        hcell_t * cur = e;
        key_t nkey = cur->key();  
        if(cur->is_node()){
          point_t c = cofm_[cur->node_idx()].coordinates(); 
          element_t r = cofm_[cur->node_idx()].radius();
          element_t l = cofm_[cur->node_idx()].lap(); 
          if(clevel == level-1){ 
            output<<"\\draw[blue] ("<<c[0]<<","<<c[1]<<") circle (0.005cm);"<<std::endl;
            //output<<"\\draw[blue] ("<<c[0]<<","<<c[1]<<") circle ("<<r<<"cm);"<<std::endl;
            output<<"\\draw[blue!60!white] ("<<c[0]<<","<<c[1]<<") circle ("<<l+r<<"cm);"<<std::endl;
          }else{
            output<<"\\draw[blue,opacity=0.2] ("<<c[0]<<","<<c[1]<<") circle (0.005cm);"<<std::endl;
            //output<<"\\draw[blue,opacity=0.2] ("<<c[0]<<","<<c[1]<<") circle ("<<r<<"cm);"<<std::endl;
            output<<"\\draw[blue!60!white,opacity=0.2] ("<<c[0]<<","<<c[1]<<") circle ("<<l+r<<"cm);"<<std::endl;
          }
          for(int i = 0 ; i < nchildren_; ++i){
            if(cur->get_child(i)){
              key_t ckey = nkey; ckey.push(i); 
              auto it = htable_.find(ckey); 
              nqueue.push_back(&(htable_.find(ckey)->second));
            } // if
          } 
        }else{
            point_t c = entities_[cur->entity_idx()].coordinates(); 
            element_t r = entities_[cur->entity_idx()].radius(); 
          if(clevel == level -1){
            output<<"\\draw[red] ("<<c[0]<<","<<c[1]<<") circle (0.005cm);"<<std::endl;
            output<<"\\draw[red] ("<<c[0]<<","<<c[1]<<") circle ("<<r<<"cm);"<<std::endl;
          }else{
            output<<"\\draw[red,opacity=0.2] ("<<c[0]<<","<<c[1]<<") circle (0.005cm);"<<std::endl;
            output<<"\\draw[red,opacity=0.2] ("<<c[0]<<","<<c[1]<<") circle ("<<r<<"cm);"<<std::endl;
          }
        }
      }
      queue = nqueue; 
      nqueue.clear(); 
      clevel++; 
    } // while 

    // Finish the file 
    output<<"\\end{tikzpicture}"<<std::endl;
    output<<"\\end{document}"<<std::endl;
    output.close(); 

  }

private:

  /**
   * Share the branches with neighbors 
   * Find the branches that are not allocated yet 
   * They are on the limit of my domain 
   */
  void share_nodes_(){

    MPI_Barrier(MPI_COMM_WORLD); 
    clog_one(trace)<<"Share nodes/entities "<<std::endl; 
    MPI_Barrier(MPI_COMM_WORLD); 

    int size, rank; 
    MPI_Comm_size(MPI_COMM_WORLD,&size); 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    MPI_Status status;

    // Do the hypercube communciation to share the branches
    // Add them in the tree in the same time
    int dim = log2(size);
    bool non_power_2 = false;
    // In case of non power two, consider dim + 1
    // In this case all rank also take size + 1 rank
    int ghosts_rank = -1;
    if(1<<dim < size){
      non_power_2 = true;
      dim++;
      ghosts_rank = rank + (1<<dim-1);
      if(ghosts_rank > (1<<dim)-1)
        ghosts_rank = rank;
    }


    std::vector<std::pair<key_t,entity_t>> ghosts_entities, r_ghosts_entities;
    std::vector<std::pair<key_t,cofm_t>> ghosts_nodes, r_ghosts_nodes; 
    int sz_pair_entities = sizeof(std::pair<key_t,entity_t>); 
    int sz_pair_nodes = sizeof(std::pair<key_t,cofm_t>); 

    int s_ge_size, s_gn_size;
    int r_ge_size, r_gn_size;  
    
    // Communication for all channels 
    for(int i = 0; i < dim; ++i){
      ghosts_entities.clear(); 
      ghosts_nodes.clear(); 
      // Find branches or entities not marked yet  
      find_nodes_(ghosts_nodes,ghosts_entities); 
      int partner = rank ^ (1<<i);
      assert(partner != rank);
      // Send entities
      s_ge_size = ghosts_entities.size()*sz_pair_entities;  
      MPI_Sendrecv(&s_ge_size,1,MPI_INT,partner,0,
                   &r_ge_size,1,MPI_INT,partner,0,MPI_COMM_WORLD,&status); 
      r_ghosts_entities.resize(r_ge_size/sz_pair_entities); 
      MPI_Sendrecv(&ghosts_entities[0],s_ge_size,MPI_BYTE,partner,0,
                    &r_ghosts_entities[0],r_ge_size,MPI_BYTE,partner,0,
                    MPI_COMM_WORLD,&status); 
      // Send nodes 
      s_gn_size = ghosts_nodes.size()*sz_pair_nodes;  
      MPI_Sendrecv(&s_gn_size,1,MPI_INT,partner,0,
                   &r_gn_size,1,MPI_INT,partner,0,MPI_COMM_WORLD,&status); 
      r_ghosts_nodes.resize(r_gn_size/sz_pair_nodes); 
      MPI_Sendrecv(&ghosts_nodes[0],s_gn_size,MPI_BYTE,partner,0,
                    &r_ghosts_nodes[0],r_gn_size,MPI_BYTE,partner,0,
                    MPI_COMM_WORLD,&status); 
      // Handle the non power two cases
      if(non_power_2)
      {
        assert(false); 
        // If this ghosts_rank exists for a real rank don't use it
        if(rank != ghosts_rank && ghosts_rank <= size-1)
          continue;
        int ghosts_partner = ghosts_rank ^ (1<<i);
        partner = ghosts_partner;
        // Case already handled before
        if(ghosts_rank == rank && partner < size)
          continue;
        if(partner >= size){
          partner -= (1<<dim-1);
        }
        if(partner == rank){
          // Add into the buffer, no communication needed
        }else{
          assert(partner != rank);
          if(ghosts_rank == rank){
            // MPI_Sendrecv(); 
            //mpi_one_to_one(rank,partner,branches,nsend,last);
          }else{
            // MPI_Sendrecv(); 
            //mpi_one_to_one(rank,partner,ghosts_branches,ghosts_nsend,
            //  ghosts_last);
          } // if 
        } // if
      } // if
      // Insert the nodes/entities in the tree 
      if(r_ghosts_entities.size() > 0){
        for(int i = 0 ; i < shared_entities_.size(); ++i){
          shared_entities_.push_back(r_ghosts_entities[i].second); 
          load_shared_entity_(i,r_ghosts_entities[i].first); 
        }
      } // if 
      if(r_ghosts_nodes.size() > 0){
        for(int i = 0 ; i < shared_nodes_.size(); ++i){
          shared_nodes_.push_back(r_ghosts_nodes[i].second); 
          load_shared_node_(i,r_ghosts_nodes[i].first); 
        }
      } // if 
      cofm(root()); 
    } // for 
    MPI_Barrier(MPI_COMM_WORLD);
    clog_one(trace)<<".done"<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
  }

  void find_nodes_(
    std::vector<std::pair<key_t,cofm_t>>& nodes, 
    std::vector<std::pair<key_t,entity_t>>& entities)
  {
    std::vector<hcell_t*> queue; 
    std::vector<hcell_t*> nqueue; 
    queue.push_back(root()); 
    while(!queue.empty()){
      for(hcell_t* e: queue){
        hcell_t * cur = e;
        key_t nkey = cur->key();
        // If this node does not have a "pointer", ship it 
        if(cur->is_unset()){
          // Add children to queue
          assert(cur->type() != 0);
          for(int j = 0 ; j < nchildren_; ++j){
            if(cur->get_child(j)){
              key_t ckey = nkey; ckey.push(j); 
              auto it = htable_.find(ckey); 
              nqueue.push_back(&(it->second)); 
            } // if 
          } // for 
        }else{
          if(cur->is_node()){
            nodes.push_back(std::make_pair(cur->key(),cofm_[cur->node_idx()])); 
          }else{
            // Add to the output 
            entities.emplace_back(cur->key(),entities_[cur->entity_idx()]); 
          } // if
        } // else  
      } // for
      queue = nqueue; 
      nqueue.clear(); 
    } // while 
  }

  void load_shared_entity_(const int& entity_idx, key_t key){
    assert(htable_.find(key) == htable_.end());
    // inserting entity 
    htable_.emplace(key,hcell_t(key,entity_idx));
    htable_.find(key)->second.set_shared();  
    // Shift key 
    int lastbit = key.pop_value(); 
    add_parent_(key,lastbit); 
  }

  void load_shared_node_(const int& node_idx, key_t key){
    assert(htable_.find(key) == htable_.end());
    // inserting node 
    htable_.emplace(key,key);
    hcell_t* cur = &(htable_.find(key)->second); 
    cur->set_shared(); 
    cur->set_node_idx(node_idx);  
    // Shift key 
    int lastbit = key.pop_value(); 
    add_parent_(key,lastbit); 
  }

  void add_parent_(key_t key, int child){
    auto parent = htable_.end(); 
    while((parent = htable_.find(key)) == htable_.end()){
      // Add parent 
      //cofm_.emplace_back(key); 
      htable_.emplace(key,key);
      parent = htable_.find(key); 
      parent->second.set_shared(); 
      parent->second.add_child(child); 
      //parent->set_node_idx(cofm_.size()-1); 
      child = key.pop_value(); 
    }// while
    // Set child in the parent found 
    //parent = &(htable_.find(key)->second); 
    assert(parent->second.node_idx() == -1); 
    parent->second.set_shared(); 
    parent->second.add_child(child); 
  }

  /**
   * Finish a branch during the creation of the tree 
   * Associate the CofM index and compute CofM.  
   */
  void finish_(const key_t& key){
    hcell_t* n = &(htable_.find(key)->second); 
    if(n->entity_idx() != -1)
      return;
    assert(n->node_idx() == -1 && n->entity_idx() == -1); 
    n->set_node_idx(cofm_.size()); 
    cofm_.emplace_back(key); 
    // Compute CoFM data 
    std::vector<hcell_t*> daughters;
    daughters.reserve(nchildren_); 
    for(int j = 0 ; j < nchildren_; ++j){
      if(n->get_child(j)){
        key_t ckey = key; ckey.push(j); 
        auto it = htable_.find(ckey); 
        daughters.push_back(&(htable_.find(ckey)->second)); 
      } // if 
    } // for 
    cofm_children_(&cofm_[n->node_idx()],daughters);
  }

  void cofm_children_(
    cofm_t* cofm, 
    std::vector<hcell_t*> daughters)
  {
    // Then compute the CoFM
    point_t coordinates = point_t{}; 
    element_t radius = 0; // bmax
    element_t mass = 0; 
    size_t sub_entities = 0; 
    element_t lap = 0; 
    // Compute the center of mass and mass  
    for(int i = 0 ; i < daughters.size(); ++i){
      if(daughters[i]->type() == 0){
        // This correspond to a body 
        int idx = daughters[i]->entity_idx(); 
        assert(idx != -1); 
        coordinates += entities_[idx].mass() * entities_[idx].coordinates(); 
        mass += entities_[idx].mass();
        ++sub_entities;       
      }else{
        // This correspond to another node 
        int idx = daughters[i]->node_idx(); 
        assert(idx != -1); 
        coordinates += cofm_[idx].mass() * cofm_[idx].coordinates(); 
        mass += cofm_[idx].mass(); 
        sub_entities += cofm_[idx].sub_entities();
      } // if 
    } // for 
    assert(mass != 0.);
    // Compute the radius 
    coordinates /= mass; 
    for(int i = 0 ; i < daughters.size(); ++i){
      if(daughters[i]->type() == 0){
        int idx = daughters[i]->entity_idx(); 
        element_t dist = distance(coordinates,entities_[idx].coordinates()); 
        radius = std::max(radius,dist);
        lap = std::max(lap,dist+entities_[idx].radius()); 
      }else{
        int idx = daughters[i]->node_idx();
        element_t dist = distance(coordinates,cofm_[idx].coordinates());
        radius = std::max(radius,dist+cofm_[idx].radius());  
        lap = std::max(lap,dist+cofm_[idx].radius()+cofm_[idx].lap());  
      }
    }// for
    // Register and quit this node 
    cofm->set_coordinates(coordinates); 
    cofm->set_radius(radius); 
    cofm->set_mass(mass); 
    cofm->set_sub_entities(sub_entities); 
    cofm->set_lap(lap); 
  }

  void exchange_boundaries_(
    const key_t& hikey, const key_t& lokey,
    key_t& hibound,key_t& lobound)
  {
    int size, rank; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Status status; 

    int up = rank==size-1?-1:rank+1; 
    int down = rank==0?-1:rank-1;

    if(rank & 1){
      if(up >= 0){
        MPI_Sendrecv(&hikey,sizeof(key_t),MPI_BYTE,up,0,
                    &hibound,sizeof(key_t),MPI_BYTE,
                    up,0,MPI_COMM_WORLD,&status); 
      }else{
        hibound = key_t::max(); 
      }
      if(down >= 0){
        MPI_Sendrecv(&lokey,sizeof(key_t),MPI_BYTE,down,0,
                    &lobound,sizeof(key_t),MPI_BYTE,
                    down,0,MPI_COMM_WORLD,&status); 
      }else{
        lobound = key_t::min(); 
      }
    }else{
      if(up >= 0){
        MPI_Sendrecv(&hikey,sizeof(key_t),MPI_BYTE,up,0,
                    &hibound,sizeof(key_t),MPI_BYTE,
                    up,0,MPI_COMM_WORLD,&status); 
      }else{
        hibound = key_t::max(); 
      }
      if(down >= 0){
        MPI_Sendrecv(&lokey,sizeof(key_t),MPI_BYTE,down,0,
                    &lobound,sizeof(key_t),MPI_BYTE,
                    down,0,MPI_COMM_WORLD,&status); 
      }else{
        lobound = key_t::min(); 
      }
    }
  } 

  size_t max_depth_;
  using umap_t = hashtable<key_t,hcell_t>; 

  typename umap_t::iterator root_; 
  umap_t htable_;

  range_t range_;

  std::vector<cofm_t> cofm_; 

  std::vector<entity_t> entities_;
  std::vector<entity_t> entities_w_;

  //const size_t max_traversal = 5;
  //std::vector<std::vector<entity_t>> ghosts_entities_;
  //size_t current_ghosts = 0;

  // Entities shared during load shared entities
  std::vector<entity_t> shared_entities_;
  std::vector<cofm_t> shared_nodes_; 

  //const int ncritical = 32;

  static constexpr int nchildren_ = (1<<dimension); 

};

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_topology_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
