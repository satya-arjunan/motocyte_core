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

#include <MesoSpecies.hpp>
#include <VolumeSpecies.hpp>
#include <Model.hpp>

MesoSpecies::MesoSpecies(const std::string name, Model& model,
                         SpaceCompartment& space_compartment,
                         const unsigned init_size, const double D):
  //set interval from, voxel_length = 6Dt (volume diffusion):
  MotileSpecies(name, model, dynamic_cast<VolumeSpecies&>(space_compartment),
                init_size, D, (D > 0)?
                space_compartment.get_meso_space().get_voxel_length()/(6*D):
                std::numeric_limits<double>::infinity()),
  space_compartment_(space_compartment) {
    get_compartment().push_species(dynamic_cast<Species&>(*this));
  }

MesoSpecies::MesoSpecies(const std::string name, Model& model,
                         const unsigned init_size, const double D):
  MesoSpecies(name, model, model.get_space_compartment(), init_size, D) {}

SpaceCompartment& MesoSpecies::get_space_compartment() {
  return space_compartment_;
}
 
std::vector<unsigned>& MesoSpecies::get_sizes() {
  return sizes_;
}

MesoSpace& MesoSpecies::get_meso_space() {
  return get_space_compartment().get_meso_space();
}

const std::vector<std::vector<double>>& MesoSpecies::get_rect_mesh() {
  return get_meso_space().get_rect_mesh();
}

void MesoSpecies::walk() {
  get_space_compartment().get_meso_space().walk(get_species_cid());
}

std::vector<unsigned>& MesoSpecies::get_occupied_voxels() {
  return get_space_compartment().get_meso_space().
    get_species_occupied_voxels(get_species_cid());
}
