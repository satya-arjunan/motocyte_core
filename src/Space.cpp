//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Motocyte
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
//

#include <Space.hpp>
#include <SpaceCompartment.hpp>

Space::Space(SpaceCompartment& space_compartment,
             const Vector<double>& dimensions, const double agent_radius,
             const int Mult, const int Div):
  //see Bansal and Ardell, Metallography 5 (1972)
  voxel_length_(std::max(double(agent_radius*2), Mult*0.554/
         cbrt(MAX_AGENTS*3.0/(dimensions.x*dimensions.y*dimensions.z)))/Div),
  dim_x_(uint16_t(std::round((dimensions.x)/voxel_length_))),
  dim_y_(uint16_t(std::round((dimensions.y)/voxel_length_))),
  dim_z_(uint16_t(std::round((dimensions.z)/voxel_length_))),
  num_voxels_(dim_x_*dim_y_*dim_z_), 
  dimensions_(Vector<double>(voxel_length_*dim_x_, voxel_length_*dim_y_,
                            voxel_length_*dim_z_)),
  dimensions_mod_(dimensions_.x-voxel_length_*0.01,
                  dimensions_.y-voxel_length_*0.01,
                  dimensions_.z-voxel_length_*0.01), 
  space_compartment_(space_compartment),
  model_(space_compartment.get_model()) {} 

const size_t Space::xyz_to_voxel_id(const Vector<double>& xyz) const {
  return coord_to_voxel_id(xyz_to_coord(xyz));
}

const double Space::get_voxel_length() const {
  return voxel_length_;
}

const Vector<double> Space::voxel_id_to_xyz(const size_t voxel_id) const {
  const Vector<unsigned> coord(voxel_id_to_coord(voxel_id));
  const Vector<double> xyz(coord.x*voxel_length_+voxel_length_/2.0,
                          coord.y*voxel_length_+voxel_length_/2.0,
                          coord.z*voxel_length_+voxel_length_/2.0);
  return xyz;
}

const Vector<double>& Space::get_dimensions_mod() const {
  return dimensions_mod_;
}

const Vector<unsigned> Space::xyz_to_coord(const Vector<double>& xyz) const {
  Vector<double> fcoord(xyz/Vector<double>(voxel_length_, voxel_length_,
                                         voxel_length_));
  Vector<unsigned> coord(unsigned(fcoord.x), unsigned(fcoord.y),
                         unsigned(fcoord.z));
  return coord;
}

const size_t Space::coord_to_voxel_id(const Vector<unsigned>& coord) const {
  return coord.x + coord.y*dim_x_+ coord.z*dim_x_*dim_y_;
}

const Vector<unsigned> Space::voxel_id_to_coord(const size_t voxel_id) const {
  Vector<unsigned> coord;
  coord.z = voxel_id/(dim_x_*dim_y_);
  coord.y = voxel_id%(dim_x_*dim_y_)/(dim_x_);
  coord.x = voxel_id%(dim_x_*dim_y_)%(dim_x_);
  return coord;
}

const Vector<double>& Space::get_dimensions() const {
  return dimensions_;
}

void Space::push_walker(Walker& walker) {
  walkers_.push_back(&walker);
}

Model& Space::get_model() {
  return model_;
}

SpaceCompartment& Space::get_space_compartment() {
  return space_compartment_;
}

Species& Space::get_species(const unsigned species_id) {
  return get_space_compartment().get_meso_species_list()[species_id].get();
}
