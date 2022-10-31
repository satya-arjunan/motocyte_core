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

#include <SpaceCompartment.hpp>
#include <Model.hpp>


SpaceCompartment::SpaceCompartment(const std::string name, Model& model,
                                   const Vector<double>& dimensions,
                                   const double agent_radius,
                                   SpaceCompartment& space_compartment,
                                   const unsigned init_size,
                                   const double D):
  VolumeSpecies(name, model, space_compartment, init_size, D), 
  micro_space_(*this, dimensions, agent_radius),
  meso_space_(*this, dimensions, agent_radius) {}

SpaceCompartment::SpaceCompartment(const std::string name, Model& model, 
                                   const Vector<double>& dimensions,
                                   const double agent_radius,
                                   const unsigned init_size,
                                   const double D):
  SpaceCompartment(name, model, dimensions, agent_radius,
                   model.get_space_compartment(), init_size, D) {}

void SpaceCompartment::push_species(Species& species) {
  VolumeSpecies::push_species(species);
  MesoSpecies* meso_species(dynamic_cast<MesoSpecies*>(&species));
  if (meso_species) {
    push_species(*meso_species);
  } else {
    MicroSpecies* micro_species(dynamic_cast<MicroSpecies*>(&species));
    if (micro_species) {
      push_species(*micro_species);
    } else {
      std::cout << "Unknown species class in SpaceCompartment" << std::endl;
    }
  }
}

void SpaceCompartment::push_species(MesoSpecies& meso_species) {
  meso_species.set_species_cid(meso_species_list_.size());
  meso_species_list_.push_back(meso_species);
  get_meso_space().push_factors(meso_species.get_sizes());
  get_meso_space().push_walker(meso_species.get_walker());
}

void SpaceCompartment::push_species(MicroSpecies& micro_species) {
  micro_species.set_species_cid(micro_species_list_.size());
  micro_species_list_.push_back(micro_species);
  get_micro_space().push_agents_xyz(micro_species.get_init_positions());
  get_micro_space().push_walker(micro_species.get_walker());
}

void SpaceCompartment::populate_all() {
  get_micro_space().populate_all();
  get_meso_space().populate_all();
  VolumeSpecies::populate_all();
}

void SpaceCompartment::update_species() {
  micro_space_.update_species();
}

const Vector<double>& SpaceCompartment::get_dimensions() const {
  return micro_space_.get_dimensions();
}

MicroSpace& SpaceCompartment::get_micro_space() {
  return micro_space_;
}

MesoSpace& SpaceCompartment::get_meso_space() {
  return meso_space_;
}

std::vector<std::reference_wrapper<MicroSpecies>>& 
SpaceCompartment::get_micro_species_list() {
  return micro_species_list_;
}

std::vector<std::reference_wrapper<MesoSpecies>>& 
SpaceCompartment::get_meso_species_list() {
  return meso_species_list_;
}






