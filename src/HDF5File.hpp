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


#ifndef __HDF5File_hpp
#define __HDF5File_hpp

#include <boost/filesystem.hpp>
#include <HDF5Interface.hpp>

class HDF5File {
public:
  HDF5File(const std::string filename, const std::string file_mode);
  ~HDF5File();
  void close();
  void flush();
  HDF5File(const std::string filename);
  hid_t h5_id() const {
    return hdf5_file_id_;
  }
  template <typename T>
  static void write_dataset(const hid_t file_handle,
                            const std::string dataset_path,
                            const std::vector<T>& data,
                            const std::pair<std::int64_t, std::int64_t> range,
                            const std::vector<std::int64_t> global_size,
                            bool use_mpio, bool use_chunking);
private:
  hid_t hdf5_file_id_;

};

#endif /* __HDF5File_hpp */
