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

#include <PointSpecies.hpp>
#include <Model.hpp>


PointSpecies::PointSpecies(const std::string name, VolumeSpecies& compartment,
             const unsigned init_size, const double D):
  MicroSpecies(name, compartment.get_model(), compartment, init_size, D),
  uni_z_(std::uniform_real_distribution<>(-1, 1)),
  uni_t_(std::uniform_real_distribution<>(0, 2*M_PI)) {
    get_compartment().push_species(dynamic_cast<Species&>(*this));
    set_size(init_size);
  }

unsigned PointSpecies::size() {
  return Species::size();
}

void PointSpecies::populate() {
  //don't use get_size, it is inaccurate:
  const unsigned compartment_init_size(get_compartment().size());
  unsigned cnt(0);
  const double agent_radius(get_compartment().get_space_compartment().
                           get_agent_radius());
  for (unsigned i(0); i != compartment_init_size; ++i) {
    compartment_points_sizes_.push_back(get_init_size());
    for (unsigned j(0); j != get_init_size(); ++j) {
      double z(uni_z_(get_model().get_rng()));
      double t(uni_t_(get_model().get_rng()));
      Vector<double> position(sqrt(1-pow(z,2))*cos(t)*agent_radius,
                             sqrt(1-pow(z,2))*sin(t)*agent_radius,
                             z*agent_radius);
      relative_positions_.push_back(position); 
      direction_vectors_.push_back(position.norm());
      ++cnt;
    }
  }
  get_populator().set_is_populated();
}

std::vector<Vector<double>>& PointSpecies::get_relative_positions() {
  return relative_positions_;
}

std::vector<Vector<double>>& PointSpecies::get_positions() {
  std::vector<Vector<double>>& comp_positions(get_compartment().get_xyz_list());
  std::vector<Vector<double>>& positions(MicroSpecies::get_xyz_list());
  unsigned cnt(0);
  positions.resize(0);
  for (unsigned i(0); i != comp_positions.size(); ++i) {
    for (unsigned j(0); j != get_init_size(); ++j) { 
      Vector<double> position(relative_positions_[cnt]+comp_positions[i]);
      //periodic boundary:
      position.mod(get_model().get_space_compartment().get_micro_space().
                   get_dimensions_mod());
      positions.push_back(position);
      ++cnt;
    }
  }
  return positions;
}


const std::vector<Vector<double>>& PointSpecies::get_direction_vectors() const {
  return direction_vectors_;
}

const std::vector<unsigned>& PointSpecies::get_compartment_ids() {
  return get_compartment().get_agent_cids();
}

const std::vector<unsigned>& PointSpecies::get_compartment_points_sizes()
  const {
  return compartment_points_sizes_;
}

void PointSpecies::add_points() {
  compartment_points_sizes_.push_back(get_init_size());
  const double agent_radius(get_compartment().get_space_compartment().
                           get_agent_radius());
  for (unsigned i(0); i != get_init_size(); ++i) {
    double z(uni_z_(get_model().get_rng()));
    double t(uni_t_(get_model().get_rng()));
    Vector<double> position(sqrt(1-pow(z,2))*cos(t)*agent_radius,
                           sqrt(1-pow(z,2))*sin(t)*agent_radius,
                           z*agent_radius);
    relative_positions_.push_back(position); 
    direction_vectors_.push_back(position.norm());
  }
  set_size_changed(true);
}

void PointSpecies::remove_points(const unsigned index) {
  const unsigned point_size(compartment_points_sizes_[index]);
  const unsigned back_start(relative_positions_.size()-point_size);
  unsigned front_start(0);
  for (unsigned i(0); i < index; ++i) {
    front_start += compartment_points_sizes_[i];
  }
  for (unsigned i(0); i < point_size; ++i) {
    relative_positions_[front_start+i] = relative_positions_[back_start+i];
    direction_vectors_[front_start+i] = direction_vectors_[back_start+i];
  }
  relative_positions_.resize(back_start);
  direction_vectors_.resize(back_start);
  compartment_points_sizes_[index] = compartment_points_sizes_.back();
  compartment_points_sizes_.pop_back();
  set_size_changed(true);
}
