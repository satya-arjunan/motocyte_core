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

#include <VolumeSpecies.hpp>
#include <Model.hpp>
#include <PointSpecies.hpp>


VolumeSpecies::VolumeSpecies(const std::string name, Model& model,
                         SpaceCompartment& space_compartment,
                         const unsigned init_size,
                         const double D):
  MicroSpecies(name, model,  dynamic_cast<VolumeSpecies&>(space_compartment),
               init_size, D),
  space_compartment_(space_compartment) {
    if (&get_space_compartment() != this) {
      get_compartment().push_species(dynamic_cast<Species&>(*this));
    }
  }

VolumeSpecies::VolumeSpecies(const std::string name, Model& model, 
                         const unsigned init_size,
                         const double D):
  VolumeSpecies(name, model, model.get_space_compartment(), init_size, D) {}

void VolumeSpecies::push_species(Species& species) {
  species_list_.push_back(species);
  VolumeSpecies* compartment(dynamic_cast<VolumeSpecies*>(&species));
  if (compartment) {
    push_species(*compartment);
  } else {
    PointSpecies* point_species(dynamic_cast<PointSpecies*>(&species));
    if (point_species) {
      push_species(*point_species);
    }
  }
}

void VolumeSpecies::set_check_collision(const bool check_collision) {
  get_space_compartment().get_micro_space().set_check_collision(
                                                    check_collision);
}



void VolumeSpecies::push_species(PointSpecies& point) {
  point_species_list_.push_back(&point);
}

void VolumeSpecies::push_species(VolumeSpecies& compartment) {
  compartment_species_list_.push_back(compartment);
}

SpaceCompartment& VolumeSpecies::get_space_compartment() {
  return space_compartment_;
}

double VolumeSpecies::get_agent_radius() const {
  return space_compartment_.get_micro_space().get_agent_radius();
}

std::vector<std::reference_wrapper<Species>>& VolumeSpecies::get_species_list() {
  return species_list_;
}

void VolumeSpecies::populate_all() {
  //std::cout << "pouplate all my name:" << get_name() << std::endl;
  populate_compartment();
  populate_child_compartments();
}

void VolumeSpecies::populate_child_compartments() {
  for (unsigned i(0); i != compartment_species_list_.size(); ++i) {
    //std::cout << "my name:" << get_name() << " child:" <<
   //   compartment_species_list_[i].get().get_name() << std::endl;
    compartment_species_list_[i].get().populate_all();
  }
}

void VolumeSpecies::populate_compartment() {
  for (unsigned i(0); i != point_species_list_.size(); ++i) {
    point_species_list_[i]->populate();
  }
}

void VolumeSpecies::walk() {
  get_space_compartment().get_micro_space().walk();
  MicroSpecies::walk();
}

std::vector<unsigned>& VolumeSpecies::get_agent_cids() {
  return agent_cids_;
}

void VolumeSpecies::remove_agent(const unsigned agent_idx) {
  remove_xyz(agent_idx);
  get_populator().remove_agent(agent_idx, size());
  get_walker().remove_agent(agent_idx);
  agent_cids_[agent_idx] = agent_cids_.back();
  agent_cids_.pop_back();
  for (unsigned i(0); i < point_species_list_.size(); ++i) {
    point_species_list_[i]->remove_points(agent_idx);
  }
}

void VolumeSpecies::add_agent(const Vector<double>& xyz) {
  add_xyz(xyz);
  get_populator().add_agent(xyz);
  get_walker().add_agent();
  for (unsigned i(0); i < point_species_list_.size(); ++i) {
    point_species_list_[i]->add_points();
  }
}







