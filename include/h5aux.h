/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

/*~--------------------------------------------------------------------------~*
 *
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file io.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Auxiliary functions for handling HDF5
 */

#pragma once

#include <cstdlib>
#include <dirent.h>
#include <iostream>
#include <libgen.h>
#include <vector>
#include <hdf5.h>


namespace h5aux {

hid_t IO_group_id; // Group id to keep track of the current step
// Data for hyperslab
hsize_t IO_offset;
hsize_t IO_count;
const int MAX_FNAME_LEN = 256;
// TODO: overload ostream instead, i.e.smth like, log_exit << "ERROR!"

MPI_Comm comm_ = MPI_COMM_WORLD;

template<typename T>
hid_t
H5P_getType(T * data) {
  hid_t type = H5T_NATIVE_INT;
  if(typeid(T) == typeid(int)) {
  }
  else if(typeid(T) == typeid(double)) {
    type = H5T_NATIVE_DOUBLE;
  }
  else if(typeid(T) == typeid(int64_t)) {
    type = H5T_NATIVE_LLONG;
  }
  else if(typeid(T) == typeid(uint64_t)) {
    type = H5T_NATIVE_ULLONG;
  }
  else {
    std::cout << "Unknown type: " << typeid(T).name() << std::endl;
    MPI_Barrier(comm_);
    MPI_Finalize();
  }
  return type;
}

hid_t
H5P_openFile(const char * filename, unsigned int flags) {
  MPI_Comm comm = comm_;
  MPI_Info info = MPI_INFO_NULL;
  /* Set up file access property list with parallel I/O access */
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);

  hid_t file_id = 0;
  if(access(filename, F_OK) != -1) {
    file_id = H5Fopen(filename, flags, plist_id);
  }
  else {
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  }

  H5Pclose(plist_id);
  return file_id;
}

void
H5P_closeFile(hid_t & file_id) {
  H5Gclose(IO_group_id);
  H5Fclose(file_id);
}

// Callback function for H5Literate
herr_t list_datasets(hid_t group_id, const char *name, const H5L_info_t *info, void *op_data) {
    H5O_info_t object_info;
    H5Oget_info_by_name(group_id, name, &object_info, H5O_INFO_BASIC, H5P_DEFAULT);

    if (object_info.type == H5O_TYPE_DATASET) {
        // Cast op_data to a vector of strings and add the dataset name
        std::vector<std::string>* datasets = static_cast<std::vector<std::string>*>(op_data);
        datasets->emplace_back(name);
    }
    return 0; // Continue iteration
}


//
// Get dimensions of a dataset in HDF5 data file
//
herr_t H5D_getDimensions(hid_t file_id, const char *dataset_name, int *ndims, hsize_t *dims) {

    // Assume that ndim and dimensions have been preallocated

    // Open the dataset
    hid_t dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
    if (dataset_id < 0) {
        std::cerr << "Error opening dataset: " << dataset_name << std::endl;
        H5Fclose(file_id);
        return -1;
    }

    // Get the dataspace of the dataset
    hid_t dataspace_id = H5Dget_space(dataset_id);
    if (dataspace_id < 0) {
        std::cerr << "Error getting dataspace for dataset: " << dataset_name << std::endl;
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        return -2;
    }

    // Get the number of dimensions in the dataspace
    *ndims = H5Sget_simple_extent_ndims(dataspace_id);
    if (*ndims < 0) {
        std::cerr << "Error getting number of dimensions for dataset: " << dataset_name << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        return -3;
    }

    // Get the size of each dimension
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

    // Close the dataspace and dataset
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    return 0;
}

//
// Read dataset; assume that *data has been allocated already
//
herr_t H5D_readDataset(const hid_t file_id, const char *dataset_name, double *data) {

    // Open the dataset
    hid_t dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
    if (dataset_id < 0) {
        std::cerr << "Error opening dataset: " << dataset_name << std::endl;
        H5Fclose(file_id);
        return -1;  // Indicate failure
    }

    // Get the dataspace of the dataset
    hid_t dataspace_id = H5Dget_space(dataset_id);
    if (dataspace_id < 0) {
        std::cerr << "Error getting dataspace for dataset: " << dataset_name << std::endl;
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        return -2;  // Indicate failure
    }

    // Get the number of dimensions and the size of each dimension
    int ndims = H5Sget_simple_extent_ndims(dataspace_id);
    if (ndims != 3) {  // Ensure the dataset is 3D
        std::cerr << "Dataset is not 3D." << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        return -3;  // Indicate failure
    }

    // Read the dataset into the preallocated array
    herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0) {
        std::cerr << "Error reading dataset: " << dataset_name << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        return -4;  // Indicate failure
    }

    // Close the dataspace and dataset
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    return 0;  // Indicate success
}




} // namespace h5aux

