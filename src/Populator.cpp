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

#include <map>
#include <Species.hpp>
#include <Populator.hpp>
#include <Model.hpp>
#include <TrackSpecies.hpp>

Populator::Populator(Model& model, MotileSpecies& motile_species,
                     Vector<double> origin, Vector<double> cuboid):
  Process(std::string("Populator "+motile_species.get_name()), model),
  motile_species_(motile_species),
  is_populated_(false),
  origin_(origin),
  cuboid_(cuboid) {}

void Populator::set_origin(const double x, const double y, const double z) {
  origin_.x = x;
  origin_.y = y;
  origin_.z = z;
}

void Populator::set_cuboid(const double x, const double y, const double z) {
  cuboid_.x = x;
  cuboid_.y = y;
  cuboid_.z = z;
}

void Populator::set_radius(const double radius) {
  radius_ = radius;
}

bool Populator::get_is_populated() const {
  return is_populated_;
}

void Populator::populate() {
  volume_species_ = dynamic_cast<VolumeSpecies*>(&motile_species_);
  track_species_ = dynamic_cast<TrackSpecies*>(&motile_species_);
  meso_species_ = dynamic_cast<MesoSpecies*>(&motile_species_);
  if (track_species_) {
    populate_track_species();
    set_is_populated();
  }
  else if (cuboid_.x != 1 || cuboid_.y != 1 || cuboid_.z != 1 || radius_ > 0) {
    if (volume_species_) {
      populate_volume_species();
      set_is_populated();
    }
    else if (meso_species_) {
      populate_meso_species();
      set_is_populated();
    }
  }
}

void Populator::set_is_populated() {
  is_populated_ = true;
}

void Populator::populate_meso_species() {
  MesoSpecies& meso_species(*meso_species_);
  MesoSpace& meso_space(meso_species.get_space_compartment().get_meso_space());
  const unsigned species_cid(meso_species.get_species_cid());
  const Vector<double> dimensions(meso_space.get_dimensions());
  const Vector<double> new_dimensions(dimensions*cuboid_);
  const Vector<double> new_origin(dimensions*origin_+dimensions/2);
  const Vector<double> min(new_origin-new_dimensions/2);
  const Vector<double> max(new_origin+new_dimensions/2);
  const double radius(radius_*dimensions.x);
  std::uniform_real_distribution<> uni_x(min.x, max.x);
  std::uniform_real_distribution<> uni_y(min.y, max.y);
  std::uniform_real_distribution<> uni_z(min.z, max.z);
  unsigned init_size(meso_species.get_init_size());
  Vector<double> xyz;
  for (unsigned agent_idx(0); agent_idx < init_size; ++agent_idx) {
    do {
      xyz = Vector<double>(uni_x(get_model().get_rng()),
                          uni_y(get_model().get_rng()),
                          uni_z(get_model().get_rng()));
    } while (radius > 0 && xyz.distance(new_origin) > radius);
    meso_space.add_factor(species_cid, meso_species, xyz);
  }
  meso_species.set_size(init_size);
}

void Populator::populate_volume_species() {
  VolumeSpecies& volume_species(*volume_species_);
  MicroSpace& micro_space(
                    volume_species.get_space_compartment().get_micro_space());
  const unsigned species_cid(volume_species.get_species_cid());
  const Vector<double> dimensions(micro_space.get_dimensions());
  const Vector<double> new_dimensions(dimensions*cuboid_);
  const Vector<double> new_origin(dimensions*origin_+dimensions/2);
  const Vector<double> min(new_origin-new_dimensions/2);
  const Vector<double> max(new_origin+new_dimensions/2);
  const double radius(radius_*dimensions.x);
  std::uniform_real_distribution<> uni_x(min.x, max.x);
  std::uniform_real_distribution<> uni_y(min.y, max.y);
  std::uniform_real_distribution<> uni_z(min.z, max.z);
  unsigned init_size(volume_species.get_init_size());
  for (unsigned agent_idx(0); agent_idx < init_size; ++agent_idx) {
    Vector<double> xyz;
    unsigned attempts(0);
    do {
      do {
        xyz = Vector<double>(uni_x(get_model().get_rng()),
                            uni_y(get_model().get_rng()),
                            uni_z(get_model().get_rng()));
      } while (radius > 0 && xyz.distance(new_origin) > radius);
      ++attempts;
    } while (!micro_space.check_collision_and_populate_agent(xyz,
               agent_idx, species_cid) && attempts < 1000);
    if (attempts >= 1000) {
      std::cout << "Not enough space to populate " << volume_species.get_name()
        << ". Only populated " << agent_idx << " agents out of " << 
        volume_species.get_init_size() << std::endl; 
      volume_species.get_init_positions().resize(agent_idx);
      volume_species.set_size(agent_idx);
      init_size = volume_species.get_size();
      break;
    }
  }
}

void Populator::remove_agent(const unsigned agent_idx,
                             const unsigned back_agent_idx) {
  VolumeSpecies& volume_species(*volume_species_);
  MicroSpace& micro_space(
                    volume_species.get_space_compartment().get_micro_space());
  const unsigned species_cid(volume_species.get_species_cid());
  micro_space.remove_agent(agent_idx, back_agent_idx, species_cid); 
}

void Populator::add_agent(const Vector<double>& xyz) {
  VolumeSpecies& volume_species(*volume_species_);
  MicroSpace& micro_space(
                    volume_species.get_space_compartment().get_micro_space());
  const unsigned agent_idx(volume_species.get_agent_cids().size());
  const unsigned species_cid(volume_species.get_species_cid());
  micro_space.populate_agent(xyz, agent_idx, species_cid);
}

//Three types of ids refering to a single agent:
//agent_idx: agent idx that is unique within the agent species. It serves
//           as the index to the agent in the agent_list of the species. For
//           fast access from lattice.
//agent_cid: agent id that is unique within the entire space compartment
//                 required for fast VisualWriter
//track id: agent id that is represented by a track
void Populator::populate_track_species() {
  TrackSpecies& track_species(*track_species_);
  MicroSpace& micro_space(
                      track_species.get_space_compartment().get_micro_space());
  const unsigned species_cid(track_species.get_species_cid());
  std::vector<Vector<double>> xyz_list(track_species.get_initial_xyz_list());
  for (unsigned agent_idx(0); agent_idx < xyz_list.size(); ++agent_idx) {
    Vector<double>& xyz(xyz_list[agent_idx]);
    if (!micro_space.check_collision_and_populate_agent(xyz, agent_idx,
                                                        species_cid)) { 
      std::cout << "Unable to populate track species because of collision " 
        << track_species.get_name() << ". Only populated " << agent_idx <<
        " agents out of " << track_species.get_init_size() << std::endl; 
      break;
    }
  }
}

