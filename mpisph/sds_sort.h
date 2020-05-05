/*
 *
 * Generally, this file sort the dataset inside a HDF5 group
 *
 */

 #pragma once 

#include "hdf5.h"
#include "stdlib.h"
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "tree.h"

const key_type lowest_key = key_type::min();
const key_type higest_key = key_type::max();

MPI_Status Stat;

// This is the struct type for orginal data
MPI_Datatype OPIC_DATA_TYPE;
MPI_Datatype OPIC_KEY_TYPE;


#if 0 
int
pivots_replicated(
  std::vector<key_type>& pivots,
  int dest,
  int * dest_pivot_replicated_size,
  int * dest_pivot_replicated_rank,
  int mpi_size,
  int mpi_rank,
  key_type * p_value_head) {
  

  std::cout<<mpi_rank<<": pivots_replicated"<<std::endl;

  int replicated = 0;
  int replicated_size = 1; // At least have one
  key_type dest_pivot;
  if(dest != (mpi_size - 1)) {
    dest_pivot = pivots[dest];
  }
  else {
    dest_pivot = higest_key;
  }

  key_type next_pivot;
  int i;
  i = dest - 1;
  while(i >= 0) {
    if(i != (mpi_size - 1)) {
      next_pivot = pivots[i];
    }
    else {
      next_pivot = higest_key;
    }

    if(next_pivot == dest_pivot) {
      i--;
      replicated = 1;
      replicated_size++;
    }
    else {
      break;
    }
  }

  if(i >= 0) {
    *p_value_head = next_pivot;
  }
  else {
    *p_value_head = lowest_key;
  }

  // Rank starts from 0
  *dest_pivot_replicated_rank = replicated_size - 1;

  i = dest + 1;
  while(i < mpi_size) {
    if(i != (mpi_size - 1)) {
      next_pivot = pivots[i];
    }
    else {
      next_pivot = higest_key;
    }
    if(next_pivot == dest_pivot) {
      i++;
      replicated = 1;
      replicated_size++;
    }
    else {
      break;
    }
  }

  *dest_pivot_replicated_size = replicated_size;

  return replicated;
}
#endif 

// Main of the parallel sampling sort
template<typename TYPE, typename COMPARE> 
int
sds_sort(std::vector<TYPE> & data, COMPARE comp, int * dist_in, 
  const int64_t& nparticles, int64_t& localnparticles) {
  int size, rank;

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  MPI_Type_contiguous(sizeof(TYPE), MPI_BYTE, &OPIC_DATA_TYPE);
  MPI_Type_commit(&OPIC_DATA_TYPE);

  MPI_Type_contiguous(sizeof(key_type), MPI_BYTE, &OPIC_KEY_TYPE);
  MPI_Type_commit(&OPIC_KEY_TYPE);

  std::vector<key_type> pivots((size - 1));
  std::vector<key_type> all_samp; // Only used by master 

  // Phase 1 ------------------------------------------------------------------
  int pass, temp_index;

  MPI_Barrier(MPI_COMM_WORLD); 
  std::sort(data.begin(),data.end(), comp);
  
  // Choose sample
  std::vector<key_type> my_sample(size);
  pass = data.size() / size;
  for(int i = 0; i < size; i++) {
    my_sample[i] = data[pass*i].key();
  }

  if(rank != 0) {
    MPI_Send(&my_sample[0], size, OPIC_KEY_TYPE, 0, 2, MPI_COMM_WORLD);
  }else{
    all_samp.insert(all_samp.end(),my_sample.begin(),my_sample.end());
  }

  // Master part 
  if(rank == 0) {

    std::vector<key_type> temp_samp(size);
    for(int i = 1; i < size; i++) {
      MPI_Recv(&temp_samp[0], size, OPIC_KEY_TYPE, MPI_ANY_SOURCE, 2,
        MPI_COMM_WORLD, &Stat);
        all_samp.insert(all_samp.end(),temp_samp.begin(),temp_samp.end());
    }

    assert(all_samp.size() == size*size); 

    std::sort(all_samp.begin(),all_samp.end());

    int rou = size/2.;
//    std::cout<<"Pivots : ";
    for(int i = 0; i < (size - 1); i++) {
      pivots[i] = all_samp[(size * (i + 1) + rou - 1)];
//      std::cout<<pivots[i]<<"; ";
    }
//    std::cout<<std::endl;
  }

  MPI_Bcast(&pivots[0], (size - 1), OPIC_KEY_TYPE, 0, MPI_COMM_WORLD);

  // Phase 2 ------------------------------------------------------------------
  int dest = 0;
  key_type cur_value;

  std::vector<TYPE> final_buff; 

  std::vector<int> scount(size * sizeof(int),0), rcount(size * sizeof(int),0);
  std::vector<int> sdisp(size * sizeof(int),0), rdisp(size * sizeof(int),0);

  int previous_ii = 0;
  for(int k = 0; k < size; k++) {
    key_type dest_pivot = (dest != (size - 1))? pivots[dest]:higest_key;
    int jj = 0;

    for(int ii = previous_ii; ii < data.size(); ii++) {
      cur_value = data[ii].key();
      previous_ii = ii;

      if(cur_value <= dest_pivot) {
        jj++;
      }
      else {
        break;
      }
    } 
    scount[dest] = jj;
    dest = (dest + 1) % size;
    if(previous_ii + 1 >= data.size())
      break;
  } // End of for(k

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);
 
  sdisp[0] = 0;
  for(int i = 1; i <= size; i++) {
    sdisp[i] = sdisp[i - 1] + scount[i - 1];
//    std::cout<<size<<": scount("<<i-1<<"): "<<scount[i-1]<<" rcount("<<i-1<<"): "<<rcount[i - 1]<<std::endl;
  }

  rdisp[0] = 0;
  for(int i = 1; i <= size; i++)
    rdisp[i] = rdisp[i - 1] + rcount[i - 1];

  unsigned long long rsize = 0;
  unsigned long long ssize = 0;
  for(int i = 0; i < size; i++) {
    rsize += (unsigned long long)rcount[i];
    ssize += (unsigned long long)scount[i];
  }

  final_buff = std::vector<TYPE>(rsize);

  MPI_Alltoallv(&data[0], &scount[0], &sdisp[0], OPIC_DATA_TYPE, &final_buff[0], &rcount[0],
    &rdisp[0], OPIC_DATA_TYPE, MPI_COMM_WORLD);

  data = final_buff; 

  std::sort(data.begin(),data.end(), comp);

  MPI_Barrier(MPI_COMM_WORLD);

  //std::cout<<rank<<"SDS_SORT done: "<<data.size()<<": "<<data.begin()->key()
  //  <<"("<<(data.begin()+1)->key()<<") - ("<<(data.end()-2)->key()<<")"<<data.back().key()<<std::endl;
  

  // Balance ------------------------------------------------------------------
  MPI_Barrier(MPI_COMM_WORLD);

  // Gather what is on everyrank 
  std::vector<int> pcount(size); 
  pcount[rank] = data.size(); 
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, &pcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  //if(rank == 0){
  //  std::cout<<"Count: -"; 
  //  for(int i = 0 ; i < size; ++i){
  //    std::cout<<pcount[i]<<"-";
  //  }
  //  std::cout<<std::endl;
  //}

  MPI_Barrier(MPI_COMM_WORLD); 

  int64_t totalparticles = 0; 
  for(int i = 0 ; i < size; ++i){
    totalparticles += pcount[i]; 
  }
  // Count the desire size
  std::vector<int> wdist(size,0);
  int lparticles = totalparticles/size; 
  int mparticles = totalparticles%size; 
  for(int i = 0 ; i < size; ++i){
    wdist[i] = lparticles; 
    if(i < mparticles){
      ++wdist[i]; 
    }
  }

  // Compute matrix of balance 
  std::vector<std::vector<int>> msend(size); 
  std::vector<std::vector<int>> mrecv(size);  

  for(int i = 0 ; i < size ; ++i){
    msend[i].resize(size); 
    mrecv[i].resize(size); 
  }

  for(int i = 0 ; i < size-1 ; ++i){
    int need = pcount[i]-wdist[i];
    if(need < 0){
      need = -need; 
      // Need more
      assert(pcount[i+1] > need);
      mrecv[i][i+1] = need;
      msend[i+1][i] = need; 
      pcount[i+1] -= need;
      pcount[i] +=  need;
    }else if(need > 0){
      // Need less
      mrecv[i+1][i] = need;
      msend[i][i+1] = need;
      pcount[i] -= need;  
      pcount[i+1] += need; 
    }
  }

  //MPI_Barrier(MPI_COMM_WORLD);
  // Display matrix 
  //if(rank == 0){
  //  for(int i = 0 ; i < size; ++i){
  //    printf("%d(%05d) | ",i,wdist[i]); 
  //    for(int j = 0 ; j < size; ++j){
  //      printf("%05d(%05d)  ",mrecv[i][j],msend[i][j]); 
  //    }
  //    printf("\n"); 
  //  }
  //}
  //MPI_Barrier(MPI_COMM_WORLD);

  MPI_Request req; 

  // Do the send/receive to balance with Isend Irecv
  for(int i = 0 ; i < size; ++i){
    if(msend[rank][i] != 0){
      // Send lower or higher part of my particles 
      if(rank < i){
        MPI_Isend(
          &data[data.size()-msend[rank][i]],
          msend[rank][i],OPIC_DATA_TYPE,i,3,MPI_COMM_WORLD,
          &req); 
      }else{
        MPI_Isend(
          &data[0],
          msend[rank][i],OPIC_DATA_TYPE,i,3,MPI_COMM_WORLD,
          &req); 
      }
    }
  }

  std::vector<TYPE> temp; 
  // Receive in a special buffer 
  for(int i = 0 ; i < size ; ++i){
    if(mrecv[rank][i] != 0){
      temp.resize(temp.size()+mrecv[rank][i]); 
      MPI_Recv(&temp[temp.size()-mrecv[rank][i]],mrecv[rank][i],OPIC_DATA_TYPE,i,3,MPI_COMM_WORLD,&Stat);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Clean local data and add temp 
  for(int i = 0 ; i < size; ++i){
    if(msend[rank][i] != 0){
      // Send lower or higher part of my particles 
      if(rank < i){
        data.erase(data.end()-msend[rank][i],data.end());
      }else{
        data.erase(data.begin(),data.begin()+msend[rank][i]);
      }
    }
  }

  data.insert(data.end(),temp.begin(),temp.end());
  std::sort(data.begin(),data.end(),comp);

  localnparticles = data.size(); 
  //std::cout<<rank<<"SDS_SORT done: "<<data.size()<<": "<<data.begin()->key()
  //  <<"("<<(data.begin()+1)->key()<<") - ("<<(data.end()-2)->key()<<")"<<data.back().key()<<std::endl;
  
  MPI_Barrier(MPI_COMM_WORLD); 

  return 0; 
}


// Phase 2 of the parallel sampling sorting
// 1, Reseive the pilots from the master
// 2, Exchange the data
// 3, Sort the data again
template<typename TYPE, typename COMPARE> 
int
phase2(int mpi_rank,
  int mpi_size,
  std::vector<TYPE>& data,
  int64_t my_data_size,
  COMPARE comp, 
  std::vector<key_type>& pivots) {

  int dest = 0;
  key_type cur_value, t_value;

  std::vector<TYPE> send_buff; 
  std::vector<TYPE> recv_buff; 
  std::vector<TYPE> final_buff; 
  std::vector<TYPE> temp_final_buff; 

  std::vector<int> scount(mpi_size * sizeof(int),0);
  std::vector<int> rcount(mpi_size * sizeof(int),0);
  std::vector<int> sdisp(mpi_size * sizeof(int),0);
  std::vector<int> rdisp(mpi_size * sizeof(int),0);

  // Start exchanging the data
  MPI_Barrier(MPI_COMM_WORLD);

  //int dest_pivot_replicated = 0; //, dest_pivot_has_same_pivot_as_me = 0;
  //int dest_pivot_replicated_size;
  //int dest_pivot_replicated_rank;
  //key_type p_value_head;

  // my_pivot starts from zero
  key_type my_pivot, dest_pivot;
  if(mpi_rank != (mpi_size - 1)) {
    my_pivot = pivots[mpi_rank];
  }
  else {
    my_pivot = higest_key; // This is the maximum value
  }

  int rank_equal, rank_less, partition_size;

  int previous_ii = 0;
  for(int k = 0; k < mpi_size; k++) {
    if(dest != (mpi_size - 1)) {
      dest_pivot = pivots[dest];
    }
    else {
      dest_pivot = higest_key;
    }
    std::cout<<mpi_rank<<" work: "<<dest_pivot<<std::endl;
    int jj = 0;

#ifdef SKEW
    // Does the data is skewed
    //if(skew_data == 1) {
    // Send data (<< dest_pivot) to my self
    if(dest_pivot == my_pivot && mpi_rank == dest) {
      // Find the data starting from previous_ii;
      for(int ii = previous_ii; ii < my_data_size; ii++) {
        cur_value = data[ii].key();
        // The data is sorted and we could skip when it bigger values appear.
        previous_ii = ii;
        if(cur_value <= dest_pivot) {
          jj++;
        }
        else {
          break;
        }
      }
    }
    else if(dest_pivot == my_pivot && mpi_rank != dest) {
      // Send nothing to the one who has same pivot as me
      jj = 0;
    }
    else {
      //"dest_pivot != my_pivot";
      dest_pivot_replicated = pivots_replicated(pivots, dest,
        &dest_pivot_replicated_size, &dest_pivot_replicated_rank, mpi_size,
        mpi_rank, &p_value_head);
      if(dest_pivot_replicated == 0) {
        // Find the data starting from previous_ii;
        if(verbose && mpi_rank == 1) //|| mpi_rank  == (mpi_size -1))
          printf("I am here: %d , %f \n", previous_ii, dest_pivot);

        for(int ii = previous_ii; ii < my_data_size; ii++) {
          cur_value = data[ii].key();
          previous_ii = ii;
          // The data is sorted and we could skip when it bigger values
          // appear.
          if(cur_value <= dest_pivot) {
            // memcpy(send_buff+jj*row_size, data+ii*row_size, row_size);
            jj++;
          }
          else {
            break;
          }
        }
      }
      else {
        rank_pivot(pivots, data, my_data_size, dest, &rank_less, &rank_equal,
          mpi_rank, p_value_head);

        if(rank_less == 0) {
          // replciated pivots are the smallest values. do equally parition
          // the data
          if(rank_equal % dest_pivot_replicated_size == 0) {
            partition_size = rank_equal / dest_pivot_replicated_size;
          }
          else {
            if(dest_pivot_replicated_rank ==
                (dest_pivot_replicated_size - 1)) {
              partition_size = rank_equal / dest_pivot_replicated_size +
                                rank_equal % dest_pivot_replicated_size;
            }
            else {
              partition_size = rank_equal / dest_pivot_replicated_size;
            }
          }
        }
        else {
          if(dest_pivot_replicated_rank == 0) {
            if(rank_less < rank_equal / (dest_pivot_replicated_size - 1)) {
              partition_size =
                (rank_equal + rank_less) / (dest_pivot_replicated_size);
            }
            else {
              partition_size = rank_less;
            }
          }
          else {
            if(dest_pivot_replicated_rank !=
                (dest_pivot_replicated_size - 1)) {
              if(rank_less < rank_equal / (dest_pivot_replicated_size - 1)) {
                partition_size =
                  (rank_equal + rank_less) / (dest_pivot_replicated_size);
              }
              else {
                partition_size = rank_equal / dest_pivot_replicated_size;
              }
            }
            else {
              if(rank_less < rank_equal / (dest_pivot_replicated_size - 1)) {
                partition_size =
                  (rank_equal + rank_less) / (dest_pivot_replicated_size) +
                  (rank_equal + rank_less) % (dest_pivot_replicated_size);
              }
              else {
                partition_size = rank_equal / dest_pivot_replicated_size +
                                  rank_equal % dest_pivot_replicated_size;
              }
            }
          }
        }
        // if(mpi_rank == 0)
        //  printf("rank less %d , rank qual %d, partition_size %d, (my rank
        //  %d, dest = %d,  p_value_head %f, dest_pivot %f ) \n", rank_less,
        //  rank_equal, partition_size, mpi_rank, dest, p_value_head,
        //  dest_pivot);

        for(ii = previous_ii; ii < my_data_size; ii++) {
          cur_value = data[ii].key();
          previous_ii = ii;

          if(jj > partition_size)
            break;

          // When the dest pivot is one of replicated keys
          if(cur_value <= dest_pivot) {
            jj++;
          }
          else {
            break;
          }
        }
      } // end of replicated = 1
    }
#endif 

    for(int ii = previous_ii; ii < my_data_size; ii++) {
      cur_value = data[ii].key();
      previous_ii = ii;

      if(cur_value <= dest_pivot) {
        jj++;
      }
      else {
        break;
      }
    } 
    scount[dest] = jj;
    dest = (dest + 1) % mpi_size;
    if(previous_ii + 1 >= my_data_size)
      break;
  } // End of for(k

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);
 
  sdisp[0] = 0;
  for(int i = 1; i <= mpi_size; i++) {
    sdisp[i] = sdisp[i - 1] + scount[i - 1];
    std::cout<<mpi_rank<<": scount("<<i-1<<"): "<<scount[i-1]<<" rcount("<<i-1<<"): "<<rcount[i - 1]<<std::endl;
  }

  rdisp[0] = 0;
  for(int i = 1; i <= mpi_size; i++)
    rdisp[i] = rdisp[i - 1] + rcount[i - 1];

  unsigned long long rsize = 0;
  unsigned long long ssize = 0;
  for(int i = 0; i < mpi_size; i++) {
    rsize += (unsigned long long)rcount[i];
    ssize += (unsigned long long)scount[i];
  }

  if(ssize != my_data_size) {
    printf("At rank %d, Size mismatch after assigning on pivots: send size "
           "(%d), my_data_size (%llu) ! \n",
      mpi_rank, ssize, my_data_size);
    exit(-1);
  }

  final_buff = std::vector<TYPE>(rsize);

  MPI_Alltoallv(&data[0], &scount[0], &sdisp[0], OPIC_DATA_TYPE, &final_buff[0], &rcount[0],
    &rdisp[0], OPIC_DATA_TYPE, MPI_COMM_WORLD);

  data = final_buff; 

  std::sort(data.begin(),data.end(), comp);

  return 0;
}

#if 0 
template<typename TYPE> 
void
rank_pivot(
  std::vector<key_type>& pivots,
  std::vector<TYPE>& data,
  int64_t my_data_size,
  int dest_pivot,
  int * rank_less,
  int * rank_equal,
  int mpi_size,
  key_type p_value_head) {
  
  int rank = 0; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  std::cout<<rank<<": rank_pivot"<<std::endl;

  key_type dest_pivot_value;
  if(dest_pivot != (mpi_size - 1)) {
    dest_pivot_value =
      pivots[dest_pivot];
  }
  else {
    dest_pivot_value = higest_key;
  }
  int less_count = 0, equal_count = 0;
  key_type cur_value;

  int i;
  for(i = 0; i < my_data_size; i++) {
    cur_value = data[i].key();

    if(cur_value < dest_pivot_value && cur_value > p_value_head) {
      less_count++;
    }

    if(cur_value == dest_pivot_value) {
      equal_count++;
    }

    if(cur_value > dest_pivot_value) {
      break;
    }
  }

  *rank_less = less_count;
  *rank_equal = equal_count;
}
#endif 