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


#ifndef __XDMFFile_hpp
#define __XDMFFile_hpp

#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <pugixml.hpp>
#include <hdf5.h>
#include <Common.hpp>
#include <HDF5File.hpp>

class XDMFFile
{
public:
  XDMFFile(const std::string filename);
  ~XDMFFile();
  void write(const std::vector<Vector<double> >& points);
  void write(const std::string name, const std::vector<Vector<double> >& points,
             const double time_step);
  void write(const std::vector<std::vector<double>>& mesh,
             const bool cell_centred,
             const std::string name,
             const std::vector<unsigned> factors,
             const double time_step);
  void close();
private:
  void add_geometry_data(pugi::xml_node& xml_node,
                         hid_t h5_id, const std::string path_prefix,
                         const std::vector<std::vector<double>>& mesh);
  template<typename T>
  void add_topology_data(pugi::xml_node& xml_node, hid_t h5_id,
                         const std::vector<std::vector<double>>& mesh);
  void add_mesh(pugi::xml_node& xml_node, hid_t h5_id,
                const std::vector<std::vector<double>>& mesh,
                const std::string path_prefix);
  std::string get_hdf5_filename(std::string xdmf_filename);
  // Add set of points to XDMF xml_node and write data
  static void add_points(pugi::xml_node& xml_node, hid_t h5_id,
                         const std::vector<Vector<double> >& points);
  // Add DataItem node to an XML node. If HDF5 is open (h5_id > 0)
  // the data is written to the HDFF5 file with the path
  // 'h5_path'. Otherwise, data is witten to the XML node and
  // 'h5_path' is ignored
  template<typename T>
  static void add_data_item(pugi::xml_node& xml_node, hid_t h5_id,
                            const std::string h5_path,
                            const T& x,
                            const std::vector<std::int64_t> dimensions,
                            const std::string number_type="");
private:
  const std::string filename_;
  std::size_t counter_;
  std::unique_ptr<pugi::xml_document> xml_doc_;
  std::unique_ptr<HDF5File> hdf5_file_;
};

  template<typename T>
    std::string container_to_string(const T& x, std::string delimiter,
                                    int precision, int linebreak=0)
  {
    std::stringstream s;
    s.precision(precision);
    if (!x.empty())
    {
      if (linebreak == 0)
      {
        s << *x.begin();
        for (auto it = x.begin() + 1; it != x.end(); ++it)
          s << delimiter << *it;
      }
      else
      {
        for (unsigned int i = 0 ; i != x.size(); ++i)
        {
          if ((i + 1)%linebreak == 0)
            s << x[i] << std::endl;
          else
            s << x[i] << delimiter;
        }
      }
    }
    return s.str();
  }

#endif /* __XDMFFile_hpp */
