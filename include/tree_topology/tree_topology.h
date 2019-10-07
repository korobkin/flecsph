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
  
  using cofm_t = cofm_u<dimension,element_t>; 
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
    // Add the root in the htable 
    htable_.emplace(key_t::root(),key_t::root());
    root_ = htable_.find(key_t::root()); 

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
    htable_.emplace(key_t::root(),key_t::root());
    root_ = htable_.find(key_t::root()); 
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
        //std::cout<<"Exploring: "<<nkey<<std::endl; 
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
            point_t d = daughters[j]->node_ptr()->coordinates(); 
            d[0] -= center[0]; 
            dist2 = d[0]*d[0]; 
            if constexpr (dimension == 2){
              d[1] -= center[1]; 
              dist2 += d[1]*d[1]; 
            }else if constexpr (dimension == 3){
              d[1] -= center[1]; 
              d[2] -= center[2]; 
              dist2 += d[1]*d[1];  
              dist2 += d[2]*d[2];
            }
            element_t extent = std::max(radius,
              daughters[j]->node_ptr()->lap())+daughters[j]->node_ptr()->radius(); 
            if(dist2 <= extent*extent){
              stk.push(daughters[j]); 
            } // if
          }else{
            element_t dist2 = 0.; 
            point_t d = daughters[j]->entity_ptr()->coordinates(); 
            d[0] -= center[0]; 
            dist2 = d[0]*d[0]; 
            if constexpr (dimension == 2){
              d[1] -= center[1]; 
              dist2 += d[1]*d[1]; 
            }else if constexpr (dimension == 3){
              d[1] -= center[1]; 
              d[2] -= center[2]; 
              dist2 += d[1]*d[1];  
              dist2 += d[2]*d[2];
            }
            element_t extent = std::max(radius, 
              daughters[j]->entity_ptr()->radius()); 
            if(dist2 <= extent*extent){
              neighbors.push_back(daughters[j]->entity_ptr()); 
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
    cofm_t* root_ptr = r->second.node_ptr(); 
    os << "Tree: "
       << "#node: " << t.htable_.size()-t.entities_.size();
    os << " depth: " << t.max_depth_;
    os << " #root_subents: " << root_ptr->sub_entities();
    //os << " center: "<<root_ptr->coordinates(); 
    //os << " mass: "<< root_ptr->mass(); 
    //os << " radius: "<< root_ptr->radius(); 
    return os;
  }

  size_t pos = 0; 


  /**
   * Loop over the bodies to insert them in the tree and construct the 
   * branches 
   **/
  void build_tree(){
    size_t nnodes = 0; 
    size_t current_depth = key_t::max_depth(); 
    // Entity keys, last and current 
    key_t lastekey = key_t(0); 
    key_t ekey; 
    // Node keys, last and Current 
    key_t lastnkey = key_t::root(); 
    key_t nkey; 
    // Current parent and value 
    hcell_t * parent = nullptr; 
    entity_t * oldptr = nullptr; 
    for(int i = 0; i < entities_.size(); ++i){
      ekey = entities_[i].key(); 
      // Compute the current node key 
      nkey = ekey; nkey.pop(current_depth);
      // While there is a difference in the current keys 
      while(nkey != lastnkey){  
        current_depth++;
        nkey = ekey; nkey.pop(current_depth);  
        lastnkey = lastekey; lastnkey.pop(current_depth); 
      } 
      parent = &(htable_.find(lastnkey)->second); 
      oldptr = parent->entity_ptr(); 
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
        parent->set_entity_ptr(nullptr); 
        htable_.emplace(nkey,nkey);
        ++nnodes; 
        parent = &(htable_.find(nkey)->second); 
      }

      // Recover deleted entity 
      if(oldptr){
        int bit = lastnkey.last_value(); 
        parent->add_child(bit); 
        parent->set_entity_ptr(nullptr);
        ++nnodes; 
        htable_.emplace(lastnkey,hcell_t(lastnkey,&(entities_[i-1])));  
      }
      // Insert new entity
      int bit = nkey.last_value(); 
      parent->add_child(bit); 
      htable_.emplace(nkey,hcell_t(nkey,&(entities_[i]))); 

      // Prepare next loop  
      lastekey = ekey; 
      lastnkey = nkey; 

    }
    cofm_.resize(nnodes+1); 
    // Call the cofm 
    pos = 0;
    cofm(root()); 
    auto r = htable_.find(key_t::root()); 

    if(entities_.size() <= 200 && dimension == 2){
      for(int i = 0 ; i < 6; ++i){
        tikz_draw("latex", i);
      }
    }

  }


  void cofm(hcell_t *current) {
    uint children = 0; 
    key_t nkey = current->key(); 
    // Just do something for the nodes, if ptr not set yet
    if(current->node_ptr() == nullptr && current->entity_ptr() == nullptr){
      // Create the cofm data 
      current->set_node_ptr(&(cofm_[pos++]));  
      hcell_t* daughters[nchildren_];
      for(int i = 0 ; i < nchildren_; ++i){
        if(current->get_child(i)){
          key_t ckey = nkey; ckey.push(i); 
          auto it = htable_.find(ckey); 
          assert(it != htable_.end()); 
          daughters[children++] = &(htable_.find(ckey)->second); 
        }
      }
      // Loop over daughters first 
      for(int i = 0 ; i < children; ++i){
        cofm(daughters[i]); 
      }
      // Then compute the CoFM
      point_t coordinates = point_t{}; 
      element_t radius = 0; // bmax
      element_t mass = 0; 
      size_t sub_entities = 0; 
      element_t lap = 0; 

      // Compute the center of mass and mass  
      for(int i = 0 ; i < children; ++i){
        if(daughters[i]->type() == 0){
          // This correspond to a body 
          entity_t* d = daughters[i]->entity_ptr(); 
          assert(d != nullptr); 
          coordinates += d->mass() * d->coordinates(); 
          mass += d->mass();
          ++sub_entities;       
        }else{
          // This correspond to another node 
          cofm_t* d = daughters[i]->node_ptr(); 
          assert(d != nullptr); 
          coordinates += d->mass() * d->coordinates(); 
          mass += d->mass(); 
          sub_entities += d->sub_entities();
        }
      } // for 
      assert(mass != 0.);
      // Compute the radius 
      coordinates /= mass; 
      for(int i = 0 ; i < children; ++i){
        if(daughters[i]->type() == 0){
          entity_t* d = daughters[i]->entity_ptr(); 
          element_t dist = distance(coordinates,d->coordinates()); 
          radius = std::max(radius,dist);
          lap = std::max(lap,dist+d->radius()); 
        }else{
          cofm_t* d = daughters[i]->node_ptr();
          element_t dist = distance(coordinates,d->coordinates());
          radius = std::max(radius,dist+d->radius());  
          lap = std::max(lap,dist+d->radius()+d->lap());  
        }
      }// for
      lap -= radius; 
      assert(lap >= 0); 
      // Register and quit this node 
      current->node_ptr()->set_coordinates(coordinates); 
      current->node_ptr()->set_radius(radius); 
      current->node_ptr()->set_mass(mass); 
      current->node_ptr()->set_sub_entities(sub_entities); 
      current->node_ptr()->set_lap(lap); 
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
          
            point_t c = cur->node_ptr()->coordinates(); 
            element_t r = cur->node_ptr()->radius();
            element_t l = cur->node_ptr()->lap(); 
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
            point_t c = cur->entity_ptr()->coordinates(); 
            element_t r = cur->entity_ptr()->radius(); 
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

  //using branch_map_t = hashtable<key_int_t,branch_t>;
  //using branch_map_t =
  //    std::unordered_map<branch_id_t, branch_t, branch_id_hasher__<key_t>>;
  //branch_map_t branch_map_;
  size_t max_depth_;
  //typename std::unordered_map<branch_id_t, branch_t,
  //                            branch_id_hasher__<key_t>>::iterator root_;

  //using umap_t =
  //    std::unordered_map<key_t, hcell_t, branch_id_hasher__<key_t>>;
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

  //std::vector<entity_t> shared_entities_;

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
