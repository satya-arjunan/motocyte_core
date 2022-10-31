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

#include <iostream>
#include <string>
#include <time.h>
#include <MicroSpecies.hpp>
#include <VolumeSpecies.hpp>
#include <Model.hpp>

MicroSpecies::MicroSpecies(const std::string name, Model& model,
                           VolumeSpecies& compartment,
                           const unsigned init_size,
                           const double D):
  MotileSpecies(name, model, compartment, init_size, D) {
    xyz_list_.resize(init_size);
  }

MicroSpecies::MicroSpecies(const std::string name, Model& model,
                           const unsigned init_size, const double D):
  MicroSpecies(name, model, model.get_space_compartment(), init_size, D) {}

std::vector<Vector<double>>& MicroSpecies::get_init_positions() {
  return xyz_list_;
}

void MicroSpecies::remove_xyz(const unsigned index) {
  xyz_list_[index] = xyz_list_.back();
  xyz_list_.pop_back();
}

void MicroSpecies::add_xyz(const Vector<double>& xyz) {
  xyz_list_.push_back(xyz);
}

std::vector<Vector<double>>& MicroSpecies::get_relative_positions() {
  return get_xyz_list();
}

std::vector<Vector<double>>& MicroSpecies::get_xyz_list() {
  get_compartment().update_species();
  return xyz_list_;
}

unsigned MicroSpecies::size() {
  return xyz_list_.size();
}

void MicroSpecies::set_size_changed(bool state) {
  is_size_changed = state;
}

void MicroSpecies::set_relative_positions_changed(bool state) {
  is_relative_positions_changed = state;
}

bool MicroSpecies::size_changed() {
  return is_size_changed;
}

bool MicroSpecies::relative_positions_changed() {
  return is_relative_positions_changed;
}

void MicroSpecies::walk() {
  set_relative_positions_changed(true);
}

