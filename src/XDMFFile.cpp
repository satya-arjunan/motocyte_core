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

#include <XDMFFile.hpp>


XDMFFile::XDMFFile(const std::string filename):
  filename_(filename),
  counter_(0),
  xml_doc_(new pugi::xml_document) {
}

XDMFFile::~XDMFFile()
{
  close();
}

void XDMFFile::close()
{
  hdf5_file_.reset();
}

std::string XDMFFile::get_hdf5_filename(std::string xdmf_filename) {
  boost::filesystem::path p(xdmf_filename);
  p.replace_extension(".h5");
  if (p.string() == xdmf_filename) {
    std::cout << 
      "XDMFile.cpp: deduce name of HDF5 file from XDMF filename clash" <<
      std::endl;
  }
  return p.string();
}


template<typename T>
void XDMFFile::add_data_item(pugi::xml_node& xml_node,
                             hid_t h5_id,
                             const std::string h5_path,
                             const T& x,
                             const std::vector<std::int64_t> shape,
                             const std::string number_type) {
  // Add DataItem node
  pugi::xml_node data_item_node = xml_node.append_child("DataItem");

  // Add dimensions attribute
  data_item_node.append_attribute("Dimensions")
    = container_to_string(shape, " ", 4).c_str();

  // Set type for topology data (needed by XDMF to prevent default to double)
  if (!number_type.empty()) {
    data_item_node.append_attribute("NumberType") = number_type.c_str();
  }

  // Add format attribute
  if (h5_id < 0) {
    data_item_node.append_attribute("Format") = "XML";
    data_item_node.append_child(pugi::node_pcdata)
      .set_value(container_to_string(x, " ", 4, shape[1]).c_str());
  } else {
    data_item_node.append_attribute("Format") = "HDF";

    // Get name of HDF5 file
    const std::string hdf5_filename = HDF5Interface::get_filename(h5_id);
    const boost::filesystem::path p(hdf5_filename);

    // Add HDF5 filename and HDF5 internal path to XML file
    const std::string xdmf_path = p.filename().string() + ":" + h5_path;
    data_item_node.append_child(pugi::node_pcdata).set_value(xdmf_path.c_str());

    // Compute total number of items and check for consistency with shape
    std::int64_t num_items_total = 1;
    for (auto n : shape) {
      num_items_total *= n;
    }
    // Compute data offset and range of values
    std::int64_t local_shape0 = x.size();
    for (std::size_t i = 1; i < shape.size(); ++i) {
      local_shape0 /= shape[i];
    }
    const std::int64_t offset = 0;
    const std::pair<std::int64_t, std::int64_t> local_range
      = {offset, offset + local_shape0};

    HDF5Interface::write_dataset(h5_id, h5_path, x, local_range, shape, false);
    // Add partitioning attribute to dataset
    std::size_t partitions(0);
    HDF5Interface::add_attribute(h5_id, h5_path, "partition", partitions);
  }
}

/*
//Write a single time-point data of point particles
void XDMFFile::write(const std::vector<Vector<double> >& points) {
  hid_t h5_id = -1;
  std::unique_ptr<HDF5File> h5_file;
  h5_file.reset(new HDF5File(get_hdf5_filename(filename_), "w")); 
  h5_id = h5_file->h5_id();
  // Create pugi doc
  xml_doc_->reset();
  // Add XDMF node and version attribute
  xml_doc_->append_child(pugi::node_doctype).set_value(
                                         "Xdmf SYSTEM \"Xdmf.dtd\" []");
  pugi::xml_node xdmf_node = xml_doc_->append_child("Xdmf");
  add_points(xdmf_node, h5_id, points);
  xml_doc_->save_file(filename_.c_str(), "  ");
}
*/

void XDMFFile::add_points(pugi::xml_node& xdmf_node, hid_t h5_id,
                          const std::vector<Vector<double> >& points) {
  xdmf_node.append_attribute("Version") = "3.0";
  xdmf_node.append_attribute("xmlns:xi") = "http://www.w3.org/2001/XInclude";
  pugi::xml_node domain_node = xdmf_node.append_child("Domain");

  // Add a Grid to the domain
  pugi::xml_node grid_node = domain_node.append_child("Grid");
  grid_node.append_attribute("GridType") = "Uniform";
  grid_node.append_attribute("Name") = "Point cloud";

  pugi::xml_node topology_node = grid_node.append_child("Topology");
  const std::int64_t n = points.size();
  topology_node.append_attribute("NumberOfElements") =
    std::to_string(n).c_str();
  topology_node.append_attribute("TopologyType") = "PolyVertex";
  topology_node.append_attribute("NodesPerElement") = 1;
  pugi::xml_node geometry_node = grid_node.append_child("Geometry");
  geometry_node.append_attribute("GeometryType") = "XYZ";
  // Pack data
  std::vector<double> x(3*n);
  for (std::int64_t i = 0; i < n; ++i) {
    x[3*i] = points[i].x;
    x[3*i+1] = points[i].y;
    x[3*i+2] = points[i].z;
  }

  const std::vector<std::int64_t> shape = {n, 3};
  add_data_item(geometry_node, h5_id, "/Points/coordinates", x, shape);
}

//Write a single time-point data of point particles
void XDMFFile::write(const std::vector<Vector<double> >& points) {
  hid_t h5_id = -1;
  std::unique_ptr<HDF5File> h5_file;
  h5_file.reset(new HDF5File(get_hdf5_filename(filename_), "w")); 
  h5_id = h5_file->h5_id();
  // Create pugi doc
  xml_doc_->reset();
  // Add XDMF node and version attribute
  xml_doc_->append_child(pugi::node_doctype).set_value(
                                         "Xdmf SYSTEM \"Xdmf.dtd\" []");
  pugi::xml_node xdmf_node = xml_doc_->append_child("Xdmf");

  xdmf_node.append_attribute("Version") = "3.0";
  xdmf_node.append_attribute("xmlns:xi") = "http://www.w3.org/2001/XInclude";
  pugi::xml_node domain_node = xdmf_node.append_child("Domain");
  // Add a Grid to the domain
  pugi::xml_node grid_node = domain_node.append_child("Grid");
  grid_node.append_attribute("GridType") = "Uniform";
  grid_node.append_attribute("Name") = "Point cloud";

  pugi::xml_node topology_node = grid_node.append_child("Topology");
  const std::int64_t n = points.size();
  topology_node.append_attribute("NumberOfElements") =
    std::to_string(n).c_str();
  topology_node.append_attribute("TopologyType") = "PolyVertex";
  topology_node.append_attribute("NodesPerElement") = 1;
  pugi::xml_node geometry_node = grid_node.append_child("Geometry");
  geometry_node.append_attribute("GeometryType") = "XYZ";
  // Pack data
  std::vector<double> x(3*n);
  for (std::int64_t i = 0; i < n; ++i) {
    x[3*i] = points[i].x;
    x[3*i+1] = points[i].y;
    x[3*i+2] = points[i].z;
  }
  const std::vector<std::int64_t> shape = {n, 3};
  add_data_item(geometry_node, h5_id, "/Points/coordinates", x, shape);

  xml_doc_->save_file(filename_.c_str(), "  ");
}

void XDMFFile::write(const std::string name,
                     const std::vector<Vector<double> >& points,
                     double time_step) {
  // Clear the pugi doc the first time
  if (counter_ == 0) {
    xml_doc_->reset();
    // Create XDMF header
    xml_doc_->append_child(pugi::node_doctype).set_value(
                                         "Xdmf SYSTEM \"Xdmf.dtd\" []");
    pugi::xml_node xdmf_node = xml_doc_->append_child("Xdmf");
    xdmf_node.append_attribute("Version") = "3.0";
    xdmf_node.append_attribute("xmlns:xi") = "http://www.w3.org/2001/XInclude";
    xdmf_node.append_child("Domain");
  }

  hid_t h5_id = -1;
  // Truncate the file the first time
  if (counter_ == 0) {
    hdf5_file_.reset(new HDF5File(get_hdf5_filename(filename_), "w"));
  } else {
    // Append to existing HDF5 file
    hdf5_file_.reset(new HDF5File(get_hdf5_filename(filename_), "a"));
  }
  h5_id = hdf5_file_->h5_id();
  pugi::xml_node xdmf_node = xml_doc_->child("Xdmf");
  pugi::xml_node domain_node = xdmf_node.child("Domain");

  // Should functions share mesh or not? By default they do not
  std::string tg_name = "TimeSeries_" + name;
  // Look for existing time series grid node with Name == tg_name
  std::string time_step_str = boost::lexical_cast<std::string>(time_step);
  pugi::xml_node timegrid_node;
  timegrid_node = domain_node.find_child_by_attribute("Grid", "Name",
                                                      tg_name.c_str());
  // Ensure that we have a time series grid node
  if (timegrid_node) {
    // Get existing mesh grid node with the correct time step if it exist
    // (otherwise null)
    std::string xpath = std::string("Grid[Time/@Value=\"") + time_step_str +
      std::string("\"]");
  } else {
    //  Create a new time series grid node with Name = tg_name
    timegrid_node = domain_node.append_child("Grid");
    timegrid_node.append_attribute("Name") = tg_name.c_str();
    timegrid_node.append_attribute("GridType") = "Collection";
    timegrid_node.append_attribute("CollectionType") = "Temporal";
  }
  pugi::xml_node grid_node = timegrid_node.append_child("Grid");
  grid_node.append_attribute("GridType") = "Uniform";
  grid_node.append_attribute("Name") = "Point cloud";

  pugi::xml_node topology_node = grid_node.append_child("Topology");
  const std::int64_t n = points.size();
  topology_node.append_attribute("NumberOfElements") =
    std::to_string(n).c_str();
  topology_node.append_attribute("TopologyType") = "PolyVertex";
  topology_node.append_attribute("NodesPerElement") = 1;
  pugi::xml_node geometry_node = grid_node.append_child("Geometry");
  geometry_node.append_attribute("GeometryType") = "XYZ";
  // Pack data
  std::vector<double> x(3*n);
  for (std::int64_t i = 0; i < n; ++i) {
    x[3*i] = points[i].x;
    x[3*i+1] = points[i].y;
    x[3*i+2] = points[i].z;
  }

  const std::string dataset_name = "/Points/coordinates"
    + std::to_string(counter_);

  const std::vector<std::int64_t> shape = {n, 3};
  add_data_item(geometry_node, h5_id, dataset_name, x, shape);
  pugi::xml_node time_node = grid_node.append_child("Time");
  time_node.append_attribute("Value") = time_step_str.c_str();

  xml_doc_->save_file(filename_.c_str(), "  ");
  hdf5_file_.reset();
  ++counter_;
}

void XDMFFile::add_mesh(pugi::xml_node& xml_node,
                        hid_t h5_id,
                        const std::vector<std::vector<double>>& mesh,
                        const std::string path_prefix) {
  pugi::xml_node grid_node = xml_node.append_child("Grid");
  grid_node.append_attribute("Name") = "Mesoscopic_Space";
  grid_node.append_attribute("GridType") = "Uniform";

  // Add topology node and attributes (including writing data)
  add_topology_data<std::int32_t>(grid_node, h5_id, mesh);

  // Add geometry node and attributes (including writing data)
  add_geometry_data(grid_node, h5_id, path_prefix, mesh);
}

void XDMFFile::write(const std::vector<std::vector<double>>& mesh,
                     const bool cell_centred,
                     const std::string name,
                     const std::vector<unsigned> factors, double time_step) {
  // Clear the pugi doc the first time
  if (counter_ == 0) {
    xml_doc_->reset();
    // Create XDMF header
    xml_doc_->append_child(pugi::node_doctype).set_value(
                                         "Xdmf SYSTEM \"Xdmf.dtd\" []");
    pugi::xml_node xdmf_node = xml_doc_->append_child("Xdmf");
    xdmf_node.append_attribute("Version") = "3.0";
    xdmf_node.append_attribute("xmlns:xi") = "http://www.w3.org/2001/XInclude";
    xdmf_node.append_child("Domain");
  }

  hid_t h5_id = -1;
  // Truncate the file the first time
  if (counter_ == 0) {
    hdf5_file_.reset(new HDF5File(get_hdf5_filename(filename_), "w"));
  } else {
    // Append to existing HDF5 file
    hdf5_file_.reset(new HDF5File(get_hdf5_filename(filename_), "a"));
  }
  h5_id = hdf5_file_->h5_id();
  pugi::xml_node xdmf_node = xml_doc_->child("Xdmf");
  pugi::xml_node domain_node = xdmf_node.child("Domain");

  // Should functions share mesh or not? By default they do not
  //std::string tg_name = "TimeSeries_" + name;
  std::string tg_name = "TimeSeries";
  // Look for existing time series grid node with Name == tg_name
  std::string time_step_str = boost::lexical_cast<std::string>(time_step);
  pugi::xml_node timegrid_node, mesh_node;
  timegrid_node = domain_node.find_child_by_attribute("Grid", "Name",
                                                      tg_name.c_str());
  // Ensure that we have a time series grid node
  bool new_timegrid = false;
  if (timegrid_node) {
    // Get existing mesh grid node with the correct time step if it exist
    // (otherwise null)
    std::string xpath = std::string("Grid[Time/@Value=\"") + time_step_str +
      std::string("\"]");
    mesh_node = timegrid_node.select_node(xpath.c_str()).node();
  } else {
    //  Create a new time series grid node with Name = tg_name
    timegrid_node = domain_node.append_child("Grid");
    timegrid_node.append_attribute("Name") = tg_name.c_str();
    timegrid_node.append_attribute("GridType") = "Collection";
    timegrid_node.append_attribute("CollectionType") = "Temporal";
    new_timegrid = true;
  }

  // Only add mesh grid node at this time step if no other function has
  // previously added it (and parameters["functions_share_mesh"] == true)
  if (!mesh_node) {
    // Add the mesh grid node to to the time series grid node
    if (new_timegrid) {
      add_mesh(timegrid_node, h5_id, mesh, "/Mesh/" + std::to_string(counter_));
    } else {
      // Make a grid node that references back to first mesh grid node of the
      // time series
      pugi::xml_node grid_node = timegrid_node.append_child("Grid");

      // Reference to previous topology and geometry document nodes via XInclude
      std::string xpointer = std::string("xpointer(//Grid[@Name=\"") + tg_name +
        std::string("\"]/Grid[1]/*[self::Topology or self::Geometry])");
      pugi::xml_node reference = grid_node.append_child("xi:include");
      reference.append_attribute("xpointer") = xpointer.c_str();
    }

    // Get the newly created mesh grid node
    mesh_node = timegrid_node.last_child();

    // Add time value to mesh grid node
    pugi::xml_node time_node = mesh_node.append_child("Time");
    time_node.append_attribute("Value") = time_step_str.c_str();
  }

  // Add attribute node
  pugi::xml_node attribute_node = mesh_node.append_child("Attribute");
  attribute_node.append_attribute("Name") = name.c_str();
  attribute_node.append_attribute("AttributeType") = "Scalar";
  attribute_node.append_attribute("Center") = cell_centred ? "Cell" : "Node";

  const std::vector<std::int64_t> shape = {int64_t(mesh[0].size()),
    int64_t(mesh[1].size()), int64_t(mesh[2].size())};

  const std::string dataset_name = "/VisualisationVector/"
                                   + std::to_string(counter_);

  std::vector<double> data;
  for (unsigned i(0); i != factors.size(); ++i) {
    data.push_back(factors[i]);
  }
  add_data_item(attribute_node, h5_id, dataset_name, data, shape);

  xml_doc_->save_file(filename_.c_str(), "  ");
  hdf5_file_.reset();
  ++counter_;
}

template<typename T>
void XDMFFile::add_topology_data(pugi::xml_node& xml_node,
                                 hid_t h5_id,
                                 const std::vector<std::vector<double>>& mesh) {
  const std::vector<std::int64_t> shape = {int64_t(mesh[0].size()),
    int64_t(mesh[1].size()), int64_t(mesh[2].size())};
  pugi::xml_node topology_node = xml_node.append_child("Topology");
  topology_node.append_attribute("TopologyType") = "3DRECTMesh";
  topology_node.append_attribute("NumberOfElements") = 
    container_to_string(shape, " ", 4).c_str();
}

void XDMFFile::add_geometry_data(pugi::xml_node& xml_node,
                                 hid_t h5_id, const std::string path_prefix,
                                 const std::vector<std::vector<double>>& mesh) {
  // Add geometry node and attributes
  pugi::xml_node geometry_node = xml_node.append_child("Geometry");
  geometry_node.append_attribute("GeometryType") = "VXVYVZ";

  std::string group_name = path_prefix + "/" + "xcoord";
  std::string h5_path = group_name + "/geometry";
  std::vector<std::int64_t> shape = {int64_t(mesh[0].size())};
  add_data_item(geometry_node, h5_id, h5_path, mesh[0], shape);

  group_name = path_prefix + "/" + "ycoord";
  h5_path = group_name + "/geometry";
  shape = {int64_t(mesh[1].size())};
  add_data_item(geometry_node, h5_id, h5_path, mesh[1], shape);

  group_name = path_prefix + "/" + "zcoord";
  h5_path = group_name + "/geometry";
  shape = {int64_t(mesh[2].size())};
  add_data_item(geometry_node, h5_id, h5_path, mesh[2], shape);
}





