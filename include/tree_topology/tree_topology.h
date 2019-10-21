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

#define _DEBUG_TREE_
#ifdef _DEBUG_TREE_ 
#warning "Tree in debug mode with assert"
#endif 


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

private: 
  struct share_entity_t{
    share_entity_t(){}
    share_entity_t(const int& o,const key_t& k,const entity_t& e): 
      owner(o), key(k), entity(e) {};
    int owner;
    key_t key; 
    entity_t entity;
  };
  struct share_node_t{
    share_node_t(){}
    share_node_t(const int& o,const key_t& k,const cofm_t& n): 
      owner(o), key(k), node(n) {};
    int owner;
    key_t key; 
    cofm_t node; 
  };
  enum COMMS: int { REQUEST = 10, REPLY_NODE = 11, REPLY_ENTITY = 12, DONE_COMMS = 13 }; 


public: 

  /*!
    Constuct a tree topology with unit coordinates, i.e. each coordinate
    dimension is in range [0, 1].
   */
  tree_topology() {}

  /**
   * @brief Destroy the tree: empty the hash-table and destroy the entities
   * lists
   */
  ~tree_topology() {}

  /**
   * Clean the tree topology but not the local bodies
   */
  void clean() {
    cofm_.clear(); 
    htable_.clear();
    shared_entities_.clear(); 
    shared_nodes_.clear(); 
    requests_keys_.clear(); 
    mpi_requests_.clear(); 
    mpi_replies_.clear(); 
    nodes_replies_.clear(); 
    entities_replies_.clear(); 
    comms_done_.clear(); 
  }



  /**
   * @brief Reset the ghosts, clean the tree and reconstruct it. 
   * Do not share the particles again, use the current version of the keys 
   */
  void reset_ghosts(bool do_share_edge = true) {
    clean(); 
    build_tree(); 
  }

  /**
   * \brief Change the range of the tree topology
   */
  void set_range(const range_t &range) { range_ = range; }

  /**
   * @brief Get the range
   */
  const std::array<point_t, 2> &range() { return range_; }

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
` * @brief Apply a function ef to the sub_cells using asynchronous comms.
  */
  template <typename EF, typename... ARGS>
  void traversal_sph(EF &&ef, ARGS &&... args) {
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    MPI_Comm_size(MPI_COMM_WORLD,&size); 

    // Perform a tree traversal applying the specified function 
    // on the neighbors of the entities 
    entities_w_ = entities_; 

    // Create a traversal queue 
    std::stack<int> stk_nonlocal; 

    comms_done_.resize(size); 
    std::fill(comms_done_.begin(),comms_done_.end(),false); 
    requests_keys_.resize(1);
    requests_keys_[0].reserve(requests_keys_max_); 
    current_ = 0; 

    int i = 0;
    bool alternate = true;  
    while(i < entities_w_.size() || !stk_nonlocal.empty()){
      //clog_one(trace)<<"STK: "<<stk_nonlocal.size()<<std::endl;
      int curid = -1; 
      if(i >= entities_w_.size()){
        alternate = false; 
      }
      if(alternate){
        curid = i++; 
        alternate = false;
      }else{
        if(!stk_nonlocal.empty()){
          curid = stk_nonlocal.top(); stk_nonlocal.pop(); 
        }else{
          if( i < entities_w_.size())
            curid = i++;
          else
            break;
        }  
        alternate = true; 
      }
      assert(curid != -1); 
      bool non_local = false; 

      // Check for requests or replies 
      check_comms_(); 

      entity_t& ent = entities_w_[curid];
      point_t center = ent.coordinates();
      element_t radius = ent.radius();
      std::vector<entity_t*> neighbors;
      neighbors.reserve(100);
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
            assert(it != htable_.end()); 
            if(it->second.is_empty_node())
            {
              non_local = true; 
              if(!it->second.requested()){
                assert(it->second.owner() != rank);
                it->second.set_requested();
                request_(it->first,it->second.owner()); 
              }
            }else{
              // Check that all children are local 
              for(int k = 0; k < nchildren_; ++k){
                if(it->second.get_child(k)){
                  key_t cckey = ckey; cckey.push(k); 
                   auto it1 = htable_.find(cckey); 
                  assert(it1 != htable_.end()); 
                  assert(!it1->second.is_unset()); 
                }
              } // for
            } // if
            daughters[children++] = &(htable_.find(ckey)->second); 
          } // if
        } // for
        if(non_local){
          stk_nonlocal.push(curid); 
          break; 
        } 
        for(int j = children-1; j >= 0; --j){
          if(daughters[j]->is_node()){
            auto c = get_node(daughters[j]); 
            element_t dist2 = 0.; 
            point_t d = c->coordinates(); 
            dist2 = (d[0]-center[0])*(d[0]-center[0]);
            if constexpr (dimension == 2){
              dist2 += (d[1]-center[1])*(d[1]-center[1]);
            }else if constexpr (dimension == 3){
              dist2 += (d[1]-center[1])*(d[1]-center[1]);
              dist2 += (d[2]-center[2])*(d[2]-center[2]);
            }
            element_t extent = std::max(radius,
              c->lap())+
              c->radius(); 
            if(dist2 <= extent*extent){
              stk.push(daughters[j]); 
            } // if
          }else{
            auto e = get_entity(daughters[j]); 
            element_t dist2 = 0.; 
            point_t d = e->coordinates(); 
            dist2 = (d[0]-center[0])*(d[0]-center[0]);
            if constexpr (dimension == 2){
              dist2 += (d[1]-center[1])*(d[1]-center[1]);
            }else if constexpr (dimension == 3){
              dist2 += (d[1]-center[1])*(d[1]-center[1]);
              dist2 += (d[2]-center[2])*(d[2]-center[2]);
            }
            element_t extent = std::max(radius, 
              e->radius()); 
            if(dist2 <= extent*extent){
              neighbors.push_back(e); 
            } // if
          } 
        } // for 
      } // while
      if(!non_local)
        ef(ent,neighbors,std::forward<ARGS>(args)...); 
    } // while
    comms_all_done_ = false;  
    MPI_Request request;
    for(int i = 0 ; i < size; ++i){
      MPI_Isend(nullptr,0,MPI_INT,i,DONE_COMMS,MPI_COMM_WORLD,&request); 
    }
    // Handle communications 
    while(!comms_all_done_){
      check_comms_(); 
    }
    entities_ = entities_w_; 
    MPI_Barrier(MPI_COMM_WORLD); 
  } // traversal_sph

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

  /**
   * @brief Compute the keys of all the entities present in the structure
   */
  void compute_keys() {
    for (size_t i = 0; i < entities_.size(); ++i) {
      entities_[i].set_key(key_t(range_,entities_[i].coordinates()));
    } // for 
  }

  /*!
    @brief eturn the tree's current max depth.
   */
  size_t max_depth() const { return max_depth_; }

  /*!
    @brief Get the root branch (depth 0).
   */
  hcell_t *root() { return &root_->second; }

  /**
   * @brief Generic information for the tree topology
   */
  friend std::ostream &operator<<(std::ostream &os, tree_topology &t) {
    auto r = t.htable_.find(key_t::root());
    cofm_t* root_ptr = r->second.is_shared()?
      &t.shared_nodes_[r->second.node_idx()]:
      &(t.cofm_[r->second.node_idx()]); 
    os << "Tree: "
       << "#node: " << t.htable_.size()-t.entities_.size();
    os << " depth: " << t.max_depth_;
    os << " #root_subents: " << root_ptr->sub_entities();
    os << " center: "<<root_ptr->coordinates(); 
    os << " mass: "<< root_ptr->mass(); 
    os << " radius: "<< root_ptr->radius(); 
    return os;
  }

  /**
   * @brief Loop over the bodies to insert them in the tree and construct the
   * nodes. 
   * 1. Exchange the boundaries with the neighbors 
   * 2. insert the entities and create the branches 
   * 2.a. If a branch is between lo-hi key, the cofm can be computed 
   * 3. The tree is ready to share entities/nodes with the neighbors 
   **/
  void build_tree(){

    int size, rank; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size); 

    /* Exchange high and low bound */
    key_t lokey = entities_[0].key(); 
    key_t hikey = entities_[entities_.size()-1].key(); 
    exchange_boundaries_(hikey,lokey,hibound_,lobound_);

    max_depth_ = 0; 
    // Add the root 
    htable_.emplace(key_t::root(),key_t::root());
    root_ = htable_.find(key_t::root()); 

    size_t current_depth = key_t::max_depth(); 
    // Entity keys, last and current 
    key_t lastekey = key_t(0); 
    if(rank != 0) lastekey = lobound_; 
    key_t ekey; 
    // Node keys, last and Current 
    key_t lastnkey = key_t::root(); 
    key_t nkey, loboundnode, hiboundnode; 
    // Current parent and value 
    hcell_t * parent = nullptr;
    int oldidx = -1; 

    bool iam0 = rank==0; 
    bool iamlast = rank==size-1; 

    assert(lobound_ <= lokey);
    assert(hibound_ >= hikey); 

    // The extra turn in the loop is to finish the missing 
    // parent of the last entity 
    for(int i = 0; i <= entities_.size(); ++i){
      if(i < entities_.size()){
        ekey = entities_[i].key(); 
        // Compute the current node key 
      }else{
        ekey = hibound_; 
      }
      nkey = ekey; nkey.pop(current_depth);
      bool loopagain = false;
      // Loop while there is a difference in the current keys 
      while(nkey != lastnkey || (iamlast && i == entities_.size())){
        loboundnode = lobound_; loboundnode.pop(current_depth); 
        hiboundnode = hibound_; hiboundnode.pop(current_depth); 
        if(loopagain && (iam0 || lastnkey > loboundnode) && 
          (iamlast || lastnkey < hiboundnode))
        {
          // This node is done, we can compute CoFM
          finish_(lastnkey); 
        }
        if(iamlast && lastnkey == key_t::root()) break;
        loopagain = true; 
        current_depth++;
        nkey = ekey; nkey.pop(current_depth);  
        lastnkey = lastekey; lastnkey.pop(current_depth); 
      } // while

      if(iamlast && i == entities_.size()) break;

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
        parent = &(htable_.find(nkey)->second); 
      } // while

      // Recover deleted entity 
      if(oldidx != -1){
        int bit = lastnkey.last_value(); 
        parent->add_child(bit); 
        parent->set_entity_idx(-1);
        htable_.emplace(lastnkey,hcell_t(lastnkey,i-1));  
      } // if


      if(i < entities_.size()){
        // Insert the new entity
        int bit = nkey.last_value(); 
        parent->add_child(bit); 
        htable_.emplace(nkey,hcell_t(nkey,i)); 
      } // if
      
      // Prepare next loop  
      lastekey = ekey; 
      lastnkey = nkey; 
      max_depth_ = std::max(max_depth_,current_depth); 
    } // for 
    share_nodes_();
  }

  /**
   * @brief Traverse the tree and draw by levels in tikz
   **/
  void tikz_draw(const char* prefix, int level){
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);  
    // Create the file 
    char filename[64];
    sprintf(filename,"%s_%05d_%05d.tex",prefix,level,rank);
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
            output<<"\\draw[blue!60!white] ("<<c[0]<<","<<c[1]<<") circle ("<<l+r<<"cm);"<<std::endl;
          }else{
            output<<"\\draw[blue,opacity=0.2] ("<<c[0]<<","<<c[1]<<") circle (0.005cm);"<<std::endl;
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

  /**
   * @brief      Export to a file the current tree in memory
   * This is useful for small number of particles to see the tree 
   * representation
   */
  void graphviz_draw(int num) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    clog_one(trace) << rank << " outputing tree file #" << num << std::endl;

    char fname[64];
    sprintf(fname, "output_graphviz_%02d_%02d.gv", rank, num);
    std::ofstream output;
    output.open(fname);
    output << "digraph G {" << std::endl << "forcelabels=true;" << std::endl;

    // Add the legend
    //output << "branch [label=\"branch\" xlabel=\"sub_entities,owner\"]"
    //       << std::endl;


    std::stack<hcell_t*> stk;
    // Get root
    stk.push(root());

    while (!stk.empty()) {
      hcell_t *cur = stk.top();
      stk.pop();
      if(cur->is_unset()){
        output << std::oct << cur->key() << std::dec << " [label=\"" 
          << std::oct << cur->key() << std::dec 
          << "\", xlabel=\"\"];" << std::endl;
        output << std::oct << cur->key() << std::dec
          << " [shape=circle,color=black]" << std::endl;
        
        // Add the child to the stack and add for display
        for (size_t i = 0; i < nchildren_; ++i) {
          if(cur->get_child(i)){
            key_t ckey = cur->key(); ckey.push(i); 
            auto it = htable_.find(ckey); 
            if(it != htable_.end()){
              stk.push(&it->second);
            }else{
              continue; 
            }
            output << std::oct << cur->key() << "->" << it->second.key() << std::dec
                 << std::endl;
          } // if
        }
      } else if (cur->is_node()) {
        int sub_ent = 0; 
        int idx = cur->node_idx(); 
        cofm_t* c = cur->is_shared()?&shared_nodes_[idx]:&cofm_[idx];
        output << std::oct << cur->key() << std::dec << " [label=\"" 
          << std::oct << cur->key() << std::dec 
          << "\", xlabel=\""<<cur->nchildren()<<","<<c->sub_entities()<<","<<cur->owner() << "\"];" << std::endl;
        if(cur->is_shared()){
          output << std::oct << cur->key() << std::dec
            << " [shape=circle,color=green]" << std::endl;
        }else{
          output << std::oct << cur->key() << std::dec
            << " [shape=circle,color=blue]" << std::endl;
        }
        
        // Add the child to the stack and add for display
        for (size_t i = 0; i < nchildren_; ++i) {
          if(cur->get_child(i)){
            key_t ckey = cur->key(); ckey.push(i); 
            auto it = htable_.find(ckey); 
            if(it != htable_.end()){
              stk.push(&it->second);
            }else{
              continue; 
            }
            output << std::oct << cur->key() << "->" << it->second.key() << std::dec
                 << std::endl;
          } // if
        }
      } else {
        output << std::oct << cur->key() << std::dec << " [label=\"" 
          << std::oct << cur->key() << std::dec
           << "\", xlabel=\""<<cur->owner()<<"\"];" << std::endl;
        if(cur->is_shared()){
          output << std::oct << cur->key() << std::dec << 
           " [shape=circle,color=grey]" << std::endl;
        }else{
          output << std::oct << cur->key() << std::dec << 
           " [shape=circle,color=red]" << std::endl;
        }
      } // if 
    } // while
    output << "}" << std::endl;
    output.close();
  }

private:

  void check_comms_(){
    int flag = 1, size, rank; 
    MPI_Status status;
    static int tree_num = 1 ;
    MPI_Comm_size(MPI_COMM_WORLD,&size); 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    bool updated_tree = false;
    // Handle all current requests 
    while(flag == 1){
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status); 
      if(flag){
        int source = status.MPI_SOURCE; 
        int tag = status.MPI_TAG; 
        int nrecv = 0; 
        #ifdef _DEBUG_TREE_
        if(tag != DONE_COMMS)
          assert(source != rank); 
        #endif 
        MPI_Get_count(&status,MPI_BYTE,&nrecv); 
        switch(tag){
          case REQUEST:
            recv_requests_(source,nrecv); 
            break; 
          case REPLY_NODE: 
            updated_tree = true; 
            recv_node_replies_(source,nrecv); 
            break; 
          case REPLY_ENTITY:
            updated_tree = true; 
            recv_entity_replies_(source,nrecv); 
            break; 
          case DONE_COMMS: 
            MPI_Recv(nullptr,0,MPI_INT,source,DONE_COMMS,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
            comms_done_[source] = true; 
            comms_all_done_ = true;
            for(int i = 0 ; i < size; ++i){
              if(!comms_done_[i]){
                comms_all_done_ = false;
                break;
              }
            }
            break; 
          default: 
            std::cerr<<"Unknown message type: "<<tag<<" source: "<<source<<std::endl;
            MPI_Finalize(); 
            exit(1);   
        } // switch 
      } // if
    } // while 
    //if(updated_tree){
    //  graphviz_draw(tree_num++); 
    //}
  }

  /**
   * @brief Request a specific cell 
   */
  void request_(const key_t& key, const int& partner){
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    mpi_requests_.push_back(MPI_Request{}); 
    requests_keys_[current_].push_back(key);
    MPI_Isend(&requests_keys_[current_][requests_keys_[current_].size()-1]
      ,sizeof(key_t),MPI_BYTE,partner,REQUEST,
      MPI_COMM_WORLD,&mpi_requests_[mpi_requests_.size()]);
    if(requests_keys_[current_].size() >= requests_keys_max_-10){
      current_++; 
      requests_keys_.resize(requests_keys_.size()+1);
      requests_keys_[current_].reserve(requests_keys_max_); 
    } // if
  }

  /** 
   * @brief Check if another rank requested nodes. 
   **/
  void recv_requests_(const int& partner, const int& nrecv){
    bool found = false; 
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    key_t key; 
    MPI_Recv(&key,sizeof(key_t),MPI_BYTE,partner,REQUEST,
      MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
    hcell_t* cur = &(htable_.find(key)->second); 
    assert(cur->is_node());
    std::vector<share_node_t> tmp_nodes_replies; 
    std::vector<share_entity_t> tmp_entities_replies; 
    for(int i = 0 ; i < nchildren_; ++i){
      if(cur->get_child(i)){
        key_t ckey = cur->key(); ckey.push(i);
        auto child = htable_.find(ckey); 
        assert(child != htable_.end());
        if(child->second.is_node()){
          tmp_nodes_replies.emplace_back(
            child->second.owner(),child->second.key(),
            *get_node(&child->second)); 
        } else if(child->second.is_entity()){
          tmp_entities_replies.emplace_back(
            child->second.owner(),child->second.key(),
            *get_entity(&child->second)); 
        } else {
          assert(false); 
        } // if
      } // if
    } // for
    if(tmp_nodes_replies.size() != 0){
      mpi_replies_.push_back(MPI_Request{}); 
      nodes_replies_.push_back(tmp_nodes_replies); 
      MPI_Isend(&nodes_replies_[nodes_replies_.size()-1][0],
        sizeof(share_node_t)*tmp_nodes_replies.size(),
        MPI_BYTE,partner,REPLY_NODE,
        MPI_COMM_WORLD,&mpi_replies_[mpi_replies_.size()]);
      found = true; 
    } // if
    if(tmp_entities_replies.size() != 0 ){
      mpi_replies_.push_back(MPI_Request{}); 
      entities_replies_.push_back(tmp_entities_replies); 
      MPI_Isend(&entities_replies_[entities_replies_.size()-1][0],
        sizeof(share_entity_t)*tmp_entities_replies.size(),
        MPI_BYTE,partner,REPLY_ENTITY,
        MPI_COMM_WORLD,&mpi_replies_[mpi_replies_.size()]);
      found = true; 
    } // if
    assert(found); 
  }

  /**
   * @brief 
   */
  void recv_entity_replies_(const int& partner, const int& nrecv){
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    int nentities = nrecv/sizeof(share_entity_t);
    std::vector<share_entity_t> recv_entities(nentities);
    MPI_Recv(&recv_entities[0],sizeof(share_entity_t)*nentities,
      MPI_BYTE,partner,REPLY_ENTITY,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    key_t pkey = recv_entities[0].key; 
    pkey.pop(); 
    auto parent = htable_.find(pkey);
    assert(parent != htable_.end()); 
    for(int i = 0 ; i < nentities; ++i){
      shared_entities_.push_back(recv_entities[i].entity); 
      #ifdef _DEBUG_TREE_
      assert(htable_.find(recv_entities[i].key) == htable_.end());  
      #endif 
      htable_.emplace(recv_entities[i].key,
        hcell_t(recv_entities[i].key,shared_entities_.size()-1));
      auto it = htable_.find(recv_entities[i].key);
      it->second.set_shared();
      it->second.set_owner(recv_entities[i].owner);
      // Change parent
      int child = recv_entities[i].key.last_value(); 
      #ifdef _DEBUG_TREE_
      key_t ckey = recv_entities[i].key; ckey.pop(); 
      assert(parent->first == ckey); 
      #endif 
      parent->second.add_child(child);
    } // for
    if(parent->second.nchildren() == get_node(&parent->second)->sub_entities())
      parent->second.unset_requested(); 
  }

  void recv_node_replies_(const int& partner, const int& nrecv){
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    int nnodes = nrecv/sizeof(share_node_t);
    std::vector<share_node_t> recv_nodes(nnodes);
    MPI_Recv(&recv_nodes[0],sizeof(share_node_t)*nnodes,
      MPI_BYTE,partner,REPLY_NODE,MPI_COMM_WORLD,MPI_STATUS_IGNORE);  
    key_t pkey = recv_nodes[0].key; 
    pkey.pop(); 
    auto parent = htable_.find(pkey);
    for(int i = 0 ; i < nnodes; ++i){
      shared_nodes_.push_back(recv_nodes[i].node); 
      #ifdef _DEBUG_TREE_
      assert(htable_.find(recv_nodes[i].key) == htable_.end());  
      #endif 
      htable_.emplace(recv_nodes[i].key,recv_nodes[i].key);
      auto it = htable_.find(recv_nodes[i].key);
      it->second.set_shared();  
      it->second.set_node_idx(shared_nodes_.size()-1); 
      it->second.set_owner(recv_nodes[i].owner); 
      // Change parent 
      int child = recv_nodes[i].key.last_value(); 
      #ifdef _DEBUG_TREE_
      key_t ckey = recv_nodes[i].key; ckey.pop(); 
      assert(parent->first == ckey); 
      #endif
      parent->second.add_child(child);    
    } // for
    if(parent->second.nchildren() == get_node(&parent->second)->sub_entities())
      parent->second.unset_requested(); 
  }

  /**
   * @brief Share the entities/nodes with neighbors 
   * Find the branches that are not allocated yet, hey are on the limit 
   * of the domain. 
   */
  void share_nodes_(){
    clog_one(trace)<<"  Sharing nodes/entities "<<std::endl; 

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

    std::vector<share_entity_t> ghosts_entities, r_ghosts_entities;
    std::vector<share_node_t> ghosts_nodes, r_ghosts_nodes; 
    const int sz_entities = sizeof(share_entity_t); 
    const int sz_nodes = sizeof(share_node_t); 

    int s_ge_size, s_gn_size;
    int r_ge_size, r_gn_size;
    
    // Communication for all channels 
    for(int i = 0; i < dim; ++i){
      ghosts_entities.clear(); 
      ghosts_nodes.clear(); 
      // Find branches or entities not marked yet  
      find_nodes_(ghosts_nodes,ghosts_entities,rank); 
      int partner = rank ^ (1<<i);
      assert(partner >= 0 && partner != rank && partner < size);
      // Send lobound and hibound and bytes for nodes/entities
      s_ge_size = ghosts_entities.size()*sz_entities;  
      s_gn_size = ghosts_nodes.size()*sz_nodes;  
      std::pair<int[2],key_t[2]> s_keys; 
      s_keys.first[0] = s_ge_size;
      s_keys.first[1] = s_gn_size;
      s_keys.second[0] = lobound_; s_keys.second[1] = hibound_;
      std::pair<int[2],key_t[2]> s_rkeys; 
      MPI_Sendrecv(&s_keys,sizeof(std::pair<int,key_t[2]>),MPI_BYTE,partner,0,
                   &s_rkeys,sizeof(std::pair<int,key_t[2]>),MPI_BYTE,partner,0,
                   MPI_COMM_WORLD,&status); 
      lobound_ = std::min(s_rkeys.second[0],lobound_);  
      hibound_ = std::max(s_rkeys.second[1],hibound_); 
      // Send entities
      r_ghosts_entities.resize(s_rkeys.first[0]/sz_entities); 
      MPI_Sendrecv(&ghosts_entities[0],s_ge_size,MPI_BYTE,partner,0,
                    &r_ghosts_entities[0],s_rkeys.first[0],MPI_BYTE,partner,0,
                    MPI_COMM_WORLD,&status); 
      // Send nodes 
      r_ghosts_nodes.resize(s_rkeys.first[1]/sz_nodes); 
      MPI_Sendrecv(&ghosts_nodes[0],s_gn_size,MPI_BYTE,partner,0,
                    &r_ghosts_nodes[0],s_rkeys.first[1],MPI_BYTE,partner,0,
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
      for(int j = 0 ; j < r_ghosts_entities.size(); ++j){
        // \TODO add check for local particle back 
        shared_entities_.push_back(r_ghosts_entities[j].entity); 
        load_shared_entity_(shared_entities_.size()-1,
          r_ghosts_entities[j].key,
          r_ghosts_entities[j].owner); 
      }
      for(int j = 0 ; j < r_ghosts_nodes.size(); ++j){     
        // \TODO add check for local particle back 
        shared_nodes_.push_back(r_ghosts_nodes[j].node); 
        load_shared_node_(shared_nodes_.size()-1,
          r_ghosts_nodes[j].key,
          r_ghosts_nodes[j].owner); 
      }
      cofm_update_(root()); 
    } // for 

    #ifdef _DEBUG_TREE_ 
    assert(root()->is_node()); 
    graphviz_draw(0);
    #endif 
    clog_one(trace)<<"Sharing nodes/entities.done"<<std::endl;
  }


  /**
  * @brief Complete the CoFM in the tree with new entities and branches
  * This function is called during the sharing of entities/nodes. 
  * The inserted parents in the tree needs to be associated with a 
  * cofm and it needs to be computed. 
  */
  void cofm_update_(hcell_t *current) {
    key_t nkey = current->key(); 
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
      for(int i = 0 ; i < daughters.size(); ++i){
        cofm_update_(daughters[i]);  
      } // for
      key_t min_key, max_key; 
      key_boundary_(nkey,min_key,max_key); 
      if(min_key >= lobound_ && max_key <= hibound_){
        current->set_shared(); 
        current->set_node_idx(shared_nodes_.size());  
        shared_nodes_.push_back(nkey);
        cofm_children_(&shared_nodes_[current->node_idx()],daughters);
      } // if
    } // if
  }

  /** 
   * @brief Find the nodes or entities to be shared with other ranks. 
   * This searches for the first nodes/entities that have an index. 
   */
  void find_nodes_(
    std::vector<share_node_t>& nodes, 
    std::vector<share_entity_t>& entities,
    const int rank)
  {
    nodes.clear(); 
    entities.clear(); 
    std::vector<hcell_t*> queue; 
    std::vector<hcell_t*> nqueue; 
    queue.push_back(root()); 
    while(!queue.empty()){
      for(hcell_t* cur: queue){
        key_t nkey = cur->key();
        if(cur->is_unset()){
          #ifdef _DEBUG_TREE_
          assert(cur->type() != 0);
          #endif 
          for(int j = 0 ; j < nchildren_; ++j){
            if(cur->get_child(j)){
              key_t ckey = nkey; ckey.push(j); 
              auto it = htable_.find(ckey); 
              nqueue.push_back(&(it->second)); 
            } // if
          } // for
        }else{
          if(cur->is_node()){
            cofm_t* cofm = get_node(cur); 
            nodes.emplace_back(cur->owner(),cur->key(),*cofm);  
          }else{
            entity_t* ent = get_entity(cur); 
            entities.emplace_back(cur->owner(),cur->key(),*ent); 
          } // if
        } // else  
      } // for
      queue.clear(); 
      queue = nqueue; 
      nqueue.clear(); 
    } // while 
  }

  /**
   * @brief Load an entity in the tree from a distant process 
   * Call the add_parent_ function to link this entity to 
   * the tree.
   **/
  void load_shared_entity_(const int& entity_idx, key_t key, const int& owner){
    #ifdef _DEBUG_TREE_
    assert(htable_.find(key) == htable_.end());
    #endif 
    htable_.emplace(key,hcell_t(key,entity_idx));
    hcell_t* cur = &(htable_.find(key)->second); 
    cur->set_shared();  
    cur->set_owner(owner);
    int lastbit = key.pop_value(); 
    add_parent_(key,lastbit,owner); 
  }

  /** 
   * @brief Load a node in the tree from a distant process 
   * Call the add_parent_ function to link this entity to 
   * the tree
   **/
  void load_shared_node_(const int& node_idx, key_t key,const int& owner){
    #ifdef _DEBUG_TREE_
    assert(htable_.find(key) == htable_.end());
    #endif 
    htable_.emplace(key,key);
    hcell_t* cur = &(htable_.find(key)->second); 
    cur->set_shared(); 
    cur->set_node_idx(node_idx); 
    cur->set_owner(owner);  
    int lastbit = key.pop_value(); 
    add_parent_(key,lastbit,owner); 
  }

  /**
   * @brief Add missing parent in the tree from a distant 
   * node or entity insertion. 
   * The new parents are empty and we add the children in
   * the types for the tree traversal 
   */
  void add_parent_(key_t key, int child, const int& owner){
    auto parent = htable_.end(); 
    while((parent = htable_.find(key)) == htable_.end()){
      htable_.emplace(key,key);
      parent = htable_.find(key); 
      parent->second.set_shared(); 
      parent->second.add_child(child); 
      parent->second.set_owner(owner); 
      child = key.pop_value(); 
    } // while
    assert(parent->second.node_idx() == -1); 
    parent->second.add_child(child); 
  }

  /**
   * @brief Finish a branch during the creation of the tree 
   * This branch is done and this rank is the only one that 
   * have information for it (middle of its local tree).
   * Add node cofm + compute cofm data with cofm_children_
   */
  void finish_(const key_t& key){
    hcell_t* n = &(htable_.find(key)->second); 
    if(n->entity_idx() != -1)
      return;
    #ifdef _DEBUG_TREE_
    assert(n->node_idx() == -1 && n->entity_idx() == -1); 
    #endif 
    n->set_node_idx(cofm_.size()); 
    cofm_.emplace_back(key);  
    std::vector<hcell_t*> daughters;
    daughters.reserve(nchildren_); 
    for(int j = 0 ; j < nchildren_; ++j){
      if(n->get_child(j)){
        key_t ckey = key; ckey.push(j); 
        auto it = htable_.find(ckey); 
        #ifdef _DEBUG_TREE_
        assert(it != htable_.end());
        #endif  
        daughters.push_back(&(htable_.find(ckey)->second)); 
      } // if 
    } // for 
    cofm_children_(&cofm_[n->node_idx()],daughters);
  }

  /**
   * @brief Return an entity linked to a cell 
   * This takes care of the local/shared entity 
   */ 
  entity_t* get_entity(const hcell_t* hc){
    #ifdef _DEBUG_TREE_
    assert(hc->is_entity()); 
    #endif
    int idx = hc->entity_idx(); 
    #ifdef _DEBUG_TREE_ 
    assert(hc->is_shared()?idx<shared_entities_.size():
      idx<entities_.size() ); 
    #endif
    return hc->is_shared()?&shared_entities_[idx]:&entities_[idx];
  }

  /**
   * @brief Return a node linked to a cell 
   * This takes care of the local/shared node 
   */ 
  cofm_t* get_node(const hcell_t* hc){
    #ifdef _DEBUG_TREE_
    assert(hc->is_node());
    #endif  
    int idx = hc->node_idx(); 
    #ifdef _DEBUG_TREE_ 
    assert(hc->is_shared()?idx<shared_nodes_.size():
      idx<cofm_.size() ); 
    #endif
    return hc->is_shared()?&shared_nodes_[idx]:&cofm_[idx];
  }

  /**
   * @brief Compute the CofM data based on the daughters of the node.
   */
  void cofm_children_(
    cofm_t* cofm, 
    const std::vector<hcell_t*>& daughters)
  {
    // Then compute the CoFM
    point_t coordinates = point_t{}; 
    element_t radius = 0; // bmax
    element_t mass = 0; 
    size_t sub_entities = 0; 
    element_t lap = 0; 
    // Compute the center of mass and mass  
    for(int i = 0 ; i < daughters.size(); ++i){
      if(daughters[i]->is_entity()){
        entity_t* ent = get_entity(daughters[i]); 
        // This correspond to a body 
        coordinates += ent->mass() * ent->coordinates(); 
        mass += ent->mass();
        ++sub_entities;
      }else{
        // This correspond to another node 
        cofm_t* c = get_node(daughters[i]); 
        coordinates += c->mass() * c->coordinates(); 
        mass += c->mass(); 
        sub_entities += c->sub_entities();
      } // if 
    } // for 
    assert(mass != 0.);
    // Compute the radius 
    coordinates /= mass; 
    for(int i = 0 ; i < daughters.size(); ++i){
      if(daughters[i]->is_entity()){
        entity_t* ent = get_entity(daughters[i]); 
        element_t dist = distance(coordinates,ent->coordinates()); 
        radius = std::max(radius,dist);
        lap = std::max(lap,dist+ent->radius()); 
      }else{
        cofm_t* c = get_node(daughters[i]); 
        element_t dist = distance(coordinates,c->coordinates());
        radius = std::max(radius,dist+c->radius());  
        lap = std::max(lap,dist+c->radius()+c->lap());  
      }
    }// for
    // Register and quit this node 
    cofm->set_coordinates(coordinates); 
    cofm->set_radius(radius); 
    cofm->set_mass(mass); 
    cofm->set_sub_entities(sub_entities); 
    cofm->set_lap(lap); 
  }

  /**
   * @brief Exchange the keys boundries with direct neighbors 
   * This will be used during the sharing of nodes and entities.
   * The first and last process have: 100.. and 177.. for low or high. 
   **/
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
    } // if
  } // exchange_boundaries_

  /**
   * @brief Find the min (1XX00..) and max (1XX77..) keys around a key 
   * in the tree. This key might not be at the max depth.
   */ 
  void key_boundary_(const key_t& key, key_t& min_key, key_t& max_key){
    key_t stop; 
    stop = key_t::min(); 
    max_key = key; 
    min_key = key; 
    while(min_key < stop){
      max_key.push((1<<dimension)-1); 
      min_key.push(0); 
    } // while
  } // key_boundary

  // Tree topology
  size_t max_depth_;
  using umap_t = hashtable<key_t,hcell_t>; 
  typename umap_t::iterator root_; 
  umap_t htable_;

  range_t range_;
  std::vector<cofm_t> cofm_; 
  std::vector<entity_t> entities_;
  std::vector<entity_t> entities_w_;
  std::vector<entity_t> shared_entities_;
  std::vector<cofm_t> shared_nodes_; 
  static constexpr int nchildren_ = (1<<dimension); 
  key_t hibound_, lobound_; 
  // Communication 
  std::vector<std::vector<key_t>> requests_keys_;
  int current_;
  std::vector<MPI_Request> mpi_requests_; 
  std::vector<MPI_Request> mpi_replies_; 
  std::vector<std::vector<share_node_t>> nodes_replies_; 
  std::vector<std::vector<share_entity_t>> entities_replies_; 
  std::vector<bool> comms_done_; 
  bool comms_all_done_; 
  const int requests_keys_max_ = 100; 
};

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_topology_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
