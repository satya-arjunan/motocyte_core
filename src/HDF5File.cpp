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

#include <HDF5File.hpp>


HDF5File::HDF5File(const std::string filename, const std::string file_mode):
  hdf5_file_id_(0) {
  const boost::filesystem::path path(filename);
  if (path.has_parent_path() &&
      !boost::filesystem::is_directory(path.parent_path())) {
    boost::filesystem::create_directories(path.parent_path());
    if (!boost::filesystem::is_directory(path.parent_path())) {
      std::cout << "Could not create directory " << 
        path.parent_path().string().c_str() << std::endl;
    }
  }
  hdf5_file_id_ = HDF5Interface::open_file(filename, file_mode);
}

HDF5File::~HDF5File() {
  close();
}

void HDF5File::close() {
  if (hdf5_file_id_ > 0)
    HDF5Interface::close_file(hdf5_file_id_);
  hdf5_file_id_ = 0;
}

void HDF5File::flush() {
  HDF5Interface::flush_file(hdf5_file_id_);
}

