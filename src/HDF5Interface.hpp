//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of the Motocyte package
//
//        Copyright (C) 2019 Satya N.V. Arjunan
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Motocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Motocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Motocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya N. V. Arjunan <satya.arjunan@gmail.com>
// based on Fenics source by Chris N. Richardson and Garth N. Wells
//


#ifndef __HDF5Interface_hpp
#define __HDF5Interface_hpp

#include <boost/filesystem.hpp>
#include <hdf5.h>
#include <Common.hpp>

#define HDF5_FAIL -1

class HDF5Interface {
public:
  static hid_t open_file(const std::string filename, const std::string mode);
  static void close_file(const hid_t hdf5_file_handle);
  static void flush_file(const hid_t hdf5_file_handle);
  static std::string get_filename(hid_t hdf5_file_handle);
  template <typename T>
  static void add_attribute(const hid_t hdf5_file_handle,
                            const std::string dataset_path,
                            const std::string attribute_name,
                            const T& attribute_value);
  template <typename T>
  static void write_dataset(const hid_t file_handle,
                            const std::string dataset_path,
                            const std::vector<T>& data,
                            const std::pair<std::int64_t, std::int64_t> range,
                            const std::vector<std::int64_t> global_size,
                            bool use_chunking);
  static void add_group(const hid_t hdf5_file_handle,
                        const std::string dataset_path);
  static bool has_group(const hid_t hdf5_file_handle,
                        const std::string group_name);
private:
  template <typename T>
  static void add_attribute_value(const hid_t dset_id,
                                  const std::string attribute_name,
                                  const T& attribute_value);
  template <typename T>
  static hid_t hdf5_type() {
    std::cout << 
      "HDF5Interface.h no specialised function for this data type " <<
      std::endl;
    return 0;
  }
};

template <> inline hid_t HDF5Interface::hdf5_type<float>() {
  return H5T_NATIVE_FLOAT;
}

template <> inline hid_t HDF5Interface::hdf5_type<double>() {
  return H5T_NATIVE_DOUBLE;
}

template <> inline hid_t HDF5Interface::hdf5_type<int>() {
  return H5T_NATIVE_INT;
}

template <> inline hid_t HDF5Interface::hdf5_type<std::int64_t>() {
  return H5T_NATIVE_INT64;
}

template <> inline hid_t HDF5Interface::hdf5_type<std::size_t>() {
  if (sizeof(std::size_t) == sizeof(unsigned long)) {
    return H5T_NATIVE_ULONG;
  } else if (sizeof(std::size_t) == sizeof(unsigned int)) {
    return H5T_NATIVE_UINT;
  } else {
    std::cout << 
      "HDF5Interface.h size of std::size_t not the same size as long or int" <<
      std::endl;
  }
  return 0;
}

template <typename T>
inline void HDF5Interface::add_attribute(const hid_t hdf5_file_handle,
                                         const std::string dataset_path,
                                         const std::string attribute_name,
                                         const T& attribute_value) {
  // Open named dataset or group
  hid_t dset_id = H5Oopen(hdf5_file_handle, dataset_path.c_str(),
                          H5P_DEFAULT);
  // Check if attribute already exists and delete if so
  htri_t has_attr = H5Aexists(dset_id, attribute_name.c_str());
  if (has_attr > 0) {
    H5Adelete(dset_id, attribute_name.c_str());
  }
  // Add attribute of appropriate type
  add_attribute_value(dset_id, attribute_name, attribute_value);
  // Close dataset or group
 H5Oclose(dset_id);
}

template<typename T>
inline void
HDF5Interface::add_attribute_value(const hid_t dset_id,
                                   const std::string attribute_name,
                                   const T& attribute_value) {
  // Create a scalar dataspace
  hid_t dataspace_id = H5Screate(H5S_SCALAR);
  const hid_t h5type = hdf5_type<T>();
  // Create attribute of type std::size_t
  hid_t attribute_id = H5Acreate2(dset_id, attribute_name.c_str(),
                                 h5type, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT);
  // Write attribute to dataset
  H5Awrite(attribute_id, h5type, &attribute_value);
  // Close dataspace
  H5Sclose(dataspace_id);
  // Close attribute
  H5Aclose(attribute_id);
}

template <typename T>
inline void
HDF5Interface::write_dataset(const hid_t file_handle,
                             const std::string dataset_path,
                             const std::vector<T>& data,
                             const std::pair<std::int64_t, std::int64_t> range,
                             const std::vector<int64_t> global_size,
                             bool use_chunking) {
  // Data rank
  const std::size_t rank = global_size.size();

  // Get HDF5 data type
  const hid_t h5type = hdf5_type<T>();

  // Hyperslab selection parameters
  std::vector<hsize_t> count(global_size.begin(), global_size.end());
  count[0] = range.second - range.first;

  // Data offsets
  std::vector<hsize_t> offset(rank, 0);
  offset[0] = range.first;

  // Dataset dimensions
  const std::vector<hsize_t> dimsf(global_size.begin(), global_size.end());

  // Create a global data space
  const hid_t filespace0 = H5Screate_simple(rank, dimsf.data(), NULL);

  // Set chunking parameters
  hid_t chunking_properties;
  if (use_chunking)
  {
    // Set chunk size and limit to 1kB min/1MB max
    hsize_t chunk_size = dimsf[0]/2;
    if (chunk_size > 1048576)
      chunk_size = 1048576;
    if (chunk_size < 1024)
      chunk_size = 1024;

    hsize_t chunk_dims[2] = {chunk_size, dimsf[1]};
    chunking_properties = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(chunking_properties, rank, chunk_dims);
  }
  else
    chunking_properties = H5P_DEFAULT;

  // Check that group exists and recursively create if required
  const std::string group_name(dataset_path, 0, dataset_path.rfind('/'));
  add_group(file_handle, group_name);

  // Create global dataset (using dataset_path)
  const hid_t dset_id = H5Dcreate2(file_handle, dataset_path.c_str(), h5type,
                                   filespace0, H5P_DEFAULT,
                                   chunking_properties, H5P_DEFAULT);
  // Close global data space
  H5Sclose(filespace0);

  // Create a local data space
  const hid_t memspace = H5Screate_simple(rank, count.data(), NULL);

  // Create a file dataspace within the global space - a hyperslab
  const hid_t filespace1 = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace1, H5S_SELECT_SET, offset.data(),
                               NULL, count.data(), NULL);

  // Set parallel access
  const hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);

  // Write local dataset into selected hyperslab
  H5Dwrite(dset_id, h5type, memspace, filespace1, plist_id,
                    data.data());

  if (use_chunking)
  {
    // Close chunking properties
    H5Pclose(chunking_properties);
  }

  // Close dataset collectively
  H5Dclose(dset_id);

  // Close hyperslab
  H5Sclose(filespace1);

  // Close local dataset
  H5Sclose(memspace);

  // Release file-access template
  H5Pclose(plist_id);
}

#endif /* __HDF5Interface_hpp */
