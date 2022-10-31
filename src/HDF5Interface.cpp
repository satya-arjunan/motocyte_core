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

#include <HDF5Interface.hpp>


hid_t HDF5Interface::open_file(const std::string filename,
                               const std::string mode) {
  // Set parallel access with communicator
  const hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hid_t file_id = HDF5_FAIL;
  if (mode == "w") {
    // Create file for write, (overwriting existing file, if present)
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  } else {
    // Check that file exists
    if (!boost::filesystem::is_regular_file(filename)) {
      std::cout << "HDF5Interface.cpp open HDF5 file which does not exist" <<
        std::endl;
    }
    if (mode == "a") {
      // Open file existing file for append
      file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
    } else if (mode == "r") {
      file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
    } else {
      std::cout << "HDF5Interface.cpp open HDF5 file unknown file mode" <<
        std::endl;
    }
  }
  H5Pclose(plist_id);
  return file_id;
}

std::string HDF5Interface::get_filename(hid_t hdf5_file_handle) {
  // Get length of filename
  ssize_t length = H5Fget_name(hdf5_file_handle, NULL, 0);
  // Allocate memory
  std::vector<char> name(length +1);
  // Retrive filename
  length = H5Fget_name(hdf5_file_handle, name.data(), length + 1);
  return std::string(name.begin(), name.end());
}

void HDF5Interface::close_file(const hid_t hdf5_file_handle) {
  H5Fclose(hdf5_file_handle);
}

void HDF5Interface::flush_file(const hid_t hdf5_file_handle) {
  H5Fflush(hdf5_file_handle, H5F_SCOPE_GLOBAL);
}

void HDF5Interface::add_group(const hid_t hdf5_file_handle,
                              const std::string group_name)
{
  std::string _group_name(group_name);

  // Cannot create the root level group
  if (_group_name.size() == 0 || _group_name == "/")
    return;

  // Prepend a slash if missing
  if (_group_name[0] != '/')
    _group_name = "/" + _group_name;

  // Starting from the root level, check and create groups if needed
  std::size_t pos = 0;
  while (pos != std::string::npos)
  {
    pos++;
    pos = _group_name.find('/', pos);
    const std::string parent_name(_group_name, 0, pos);

    if (!has_group(hdf5_file_handle, parent_name))
    {
      hid_t group_id_vis = H5Gcreate2(hdf5_file_handle, parent_name.c_str(),
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group_id_vis);
    }
  }
}

bool HDF5Interface::has_group(const hid_t hdf5_file_handle,
                              const std::string group_name) {
  hid_t lapl_id = H5Pcreate(H5P_LINK_ACCESS);
  htri_t link_status = H5Lexists(hdf5_file_handle, group_name.c_str(), lapl_id);
  if(link_status==0) {
    // Close link access properties
    H5Pclose(lapl_id);
    return false;
  }

  H5O_info_t object_info;
  H5Oget_info_by_name(hdf5_file_handle, group_name.c_str(), &object_info,
                      lapl_id);

  // Close link access properties
  H5Pclose(lapl_id);

  return (object_info.type == H5O_TYPE_GROUP);
}
