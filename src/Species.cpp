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

#include <Species.hpp>
#include <Model.hpp>

Species::Species(const std::string name, Model& model,
                 VolumeSpecies& compartment, const unsigned init_size):
  init_size_(init_size),
  compartment_(compartment),
  name_(get_init_name(name)),
  model_(model),
  species_id_(model.push_species(*this)) {
    if (&compartment != this) {
    }
  }

Species::Species(const std::string name, Model& model,
                 const unsigned init_size):
  Species(name, model, model.get_space_compartment(), init_size) {}

unsigned Species::get_species_id() const {
  return species_id_;
}

Model& Species::get_model() {
  return model_;
}

VolumeSpecies& Species::get_compartment() {
  return compartment_;
}

const std::string& Species::get_name() const {
  return name_;
}

void Species::set_species_cid(const unsigned species_cid) {
  species_cid_ = species_cid;
}

unsigned Species::get_species_cid() const {
  return species_cid_;
}

std::string Species::get_name_id() const {
  std::stringstream sid;
  sid << get_species_id();
  return std::string(get_name()+":"+sid.str());
}

unsigned Species::get_init_size() const {
  return init_size_;
}

unsigned& Species::get_size() {
  return size_;
}

unsigned Species::size() {
  return size_;
}

void Species::set_size(const unsigned size) {
  size_ = size;
}

std::string Species::get_init_name(const std::string name) const {
  if (&compartment_ == this) {
    return name;
  }
  return std::string(compartment_.get_name()+"/"+name);
}

