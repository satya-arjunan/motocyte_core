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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURxyzE.
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

#include <cmath>
#include <cstring>
#include <random>
#include <algorithm>
#include <Motocyte.hpp>
#include <VolumeSpecies.hpp>
#include <Model.hpp>
#include <Walker.hpp>
#include <MicroSpace.hpp>

MicroSpace::MicroSpace(SpaceCompartment& space_compartment,
                       const Vector<double>& dimensions,
                       const double agent_radius): 
  Space(space_compartment, dimensions, agent_radius, 1, 1),
  species_list_(space_compartment.get_micro_species_list()),
  agent_radius_(agent_radius) {
    occupied_sizes_.resize(num_voxels_, 0);
    pending_occupied_sizes_.resize(num_voxels_, 0);
  }

void MicroSpace::push_agents_xyz(std::vector<Vector<double> >& agents_xyz) {
  species_agents_xyz_.push_back(agents_xyz);
}

double MicroSpace::get_agent_radius() const {
  return agent_radius_;
}

void MicroSpace::set_check_collision(const bool check_collision) {
  is_check_collision_ = check_collision;
}

void MicroSpace::populate_all() {
  std::vector<std::reference_wrapper<MicroSpecies>> species_list;
  for (unsigned i(0); i != species_list_.size(); ++i) {
    //populate using Populator with bounding box (MotileSpecies::populate):
    species_list_[i].get().populate();
    //if it is still not populated:
    if (!species_list_[i].get().get_is_populated()) {
      species_list.push_back(species_list_[i].get());
    }
  }
  unsigned total_agents(0);
  for (unsigned i(0); i != species_list.size(); ++i) {
    total_agents += species_list[i].get().get_init_size();
  }
  if (total_agents > num_voxels_) {
    std::cout << "Too many agents to populate, insufficient number of voxels."
      << std::endl << "Number of voxel:" << num_voxels_ << std::endl <<
      "Number of agents:" << total_agents << std::endl << 
      "MAX_AGENTS:" << MAX_AGENTS << std::endl << 
      "Change MAX_AGENTS in Motocyte.hpp to the max number of agents." <<
      std::endl;
    exit(0);
  }

  const Vector<unsigned> dim(dimensions_.x/(agent_radius_*2),
                             dimensions_.y/(agent_radius_*2),
                             dimensions_.z/(agent_radius_*2));
  const unsigned size(dim.x*dim.y*dim.z);
  std::vector<unsigned> voxel_ids;
  voxel_ids.resize(size);
  for (unsigned i(0); i != size; ++i) {
    voxel_ids[i] = i;
  }

  //extra spots for collision avoidance
  total_agents = std::min(total_agents*3, size);

  std::vector<unsigned> sampled_voxel_ids;
  std::sample(voxel_ids.begin(), voxel_ids.end(), 
      std::back_inserter(sampled_voxel_ids), total_agents,
      space_compartment_.get_model().get_rng());
  std::shuffle(sampled_voxel_ids.begin(), sampled_voxel_ids.end(), 
               space_compartment_.get_model().get_rng());

  unsigned cnt(0);
  for (unsigned species_cid(0); species_cid != species_list_.size();
       ++species_cid) {
    if (!species_list_[species_cid].get().get_is_populated()) {
      for(unsigned agent_idx(0);
          agent_idx != species_list_[species_cid].get().get_init_size();
          ++agent_idx) { 
        Vector<double> xyz;
        do {
          unsigned voxel_id(sampled_voxel_ids[cnt]);
          Vector<unsigned> coord;
          coord.z = voxel_id/(dim.x*dim.y);
          coord.y = voxel_id%(dim.x*dim.y)/(dim.x);
          coord.x = voxel_id%(dim.x*dim.y)%(dim.x);
          xyz = Vector<double>(coord.x*agent_radius_*2+agent_radius_,
                              coord.y*agent_radius_*2+agent_radius_,
                              coord.z*agent_radius_*2+agent_radius_);
          xyz.mod(get_dimensions_mod());
          ++cnt;
        } while(!check_collision_and_populate_agent(xyz, agent_idx, species_cid)
                && cnt < total_agents);
      }
    }
  }
  pending_occupied_sizes_ = occupied_sizes_;
  std::cout << "Populated agents uniformly" << std::endl;
}

unsigned MicroSpace::agent_idx_to_global_id(const unsigned agent_idx,
                                           const unsigned species_cid) const {
  return id_stride_*species_cid+agent_idx;
}

unsigned MicroSpace::global_id_to_agent_idx(const unsigned global_id,
                                           unsigned& species_cid) const {
  species_cid = global_id/id_stride_;
  return global_id%id_stride_;
}

void MicroSpace::set_agent_cid(const unsigned species_cid) {
  MicroSpecies& micro_species(species_list_[species_cid]);
  VolumeSpecies* volume_species(dynamic_cast<VolumeSpecies*>(&micro_species));
  volume_species->get_agent_cids().push_back(
          space_compartment_.get_model().get_next_agent_cid());
}

bool MicroSpace::populate_agent(const Vector<double>& xyz,
                                const unsigned global_id) {
  const unsigned voxel_id(xyz_to_voxel_id(xyz));
  ++occupied_sizes_[voxel_id];
  pending_occupied_sizes_[voxel_id] = occupied_sizes_[voxel_id];
  const uint16_t size(occupied_sizes_[voxel_id]);
  if (xyz_.size() < size) {
    xyz_.resize(size);
    xyz_.back().resize(num_voxels_);
    global_ids_.resize(size);
    global_ids_.back().resize(num_voxels_);
  }
  xyz_[size-1][voxel_id] = xyz;
  global_ids_[size-1][voxel_id] = global_id;
  set_agent_cid(global_id/id_stride_);
  return true;
}

bool MicroSpace::check_collision_and_populate_agent(const Vector<double>& xyz,
                                                    const size_t agent_idx,
                                                    const size_t species_cid) {
  if (is_colliding(xyz)) {
    return false;
  }
  populate_agent(xyz, agent_idx_to_global_id(agent_idx, species_cid));
  return true;
}


//populate an agent without checking for overlap (used by TrackSpecies):
void MicroSpace::populate_agent(const Vector<double>& xyz,
                                const size_t agent_idx,
                                const size_t species_cid) {
  populate_agent(xyz, agent_idx_to_global_id(agent_idx, species_cid));
}

void MicroSpace::check_overlap() const {
  std::cout << "in overlap" << std::endl;
  for (size_t i(0); i < num_voxels_; ++i) {
    for (size_t j(0); j < occupied_sizes_[i]; ++j) {
      for (size_t k(0); k < num_voxels_; ++k) {
        for (size_t l(0); l < occupied_sizes_[k]; ++l) {
          Vector<double> agent1(xyz_[j][i]);
          Vector<double> agent2(xyz_[l][k]);
          if(i != k || j != l) {
            if (agent1.distance(agent2) < agent_radius_*1.99) {
              Vector<unsigned> c1(xyz_to_coord(agent1));
              Vector<unsigned> c2(xyz_to_coord(agent2));
              std::cout << "Overlap distance:" << agent1.distance(agent2) << 
                std::endl;
              return;
            }
          }
        }
      }
    }
  }
  std::cout << "out overlap" << std::endl;
}

bool MicroSpace::is_colliding(const Vector<double>& new_xyz) const {
  const size_t voxel_id(xyz_to_voxel_id(new_xyz));
  return is_colliding(new_xyz, voxel_id);
}

bool MicroSpace::is_colliding(const Vector<double>& new_xyz,
                              const size_t voxel_id) const {
  if (!is_check_collision_) {
    return false;
  }
  const Vector<unsigned> new_coord(xyz_to_coord(new_xyz));
  Vector<unsigned> neighbor_coord;
  Vector<double> stride(0,0,0);
  for (int z(int(new_coord.z)-1); z <= int(new_coord.z)+1; ++z) {
    if(z < 0) {
      neighbor_coord.z = z + dim_z_;
      stride.z = -dim_z_*voxel_length_;
    } else if (z >= dim_z_) {
      neighbor_coord.z = z - dim_z_;
      stride.z = dim_z_*voxel_length_;
    } else {
      neighbor_coord.z = z;
      stride.z = 0;
    }
    for (int y(int(new_coord.y)-1); y <= int(new_coord.y)+1; ++y) {
      if(y < 0) {
        neighbor_coord.y = y + dim_y_;
        stride.y = -dim_y_*voxel_length_;
      } else if (y >= dim_y_) {
        neighbor_coord.y = y - dim_y_;
        stride.y = dim_y_*voxel_length_;
      } else {
        neighbor_coord.y = y;
        stride.y = 0;
      }
      for (int x(int(new_coord.x)-1); x <= int(new_coord.x)+1; ++x) {
        if(x < 0) {
          neighbor_coord.x = x + dim_x_;
          stride.x = -dim_x_*voxel_length_;
        } else if (x >= dim_x_) {
          neighbor_coord.x = x - dim_x_;
          stride.x = dim_x_*voxel_length_;
        } else {
          neighbor_coord.x = x;
          stride.x = 0;
        }
        const size_t neighbor_voxel_id(coord_to_voxel_id(neighbor_coord));
        const uint16_t agent_size(occupied_sizes_[neighbor_voxel_id]);
        if (agent_size != 0) {
          for (uint16_t j(0); j < agent_size; ++j) {
            if(new_xyz.distance(xyz_[j][neighbor_voxel_id]+stride) <
               agent_radius_*2) {
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}

unsigned MicroSpace::is_colliding(const Vector<double>& new_xyz,
                                  const size_t voxel_id,
                                  const uint16_t current_idx) const {
  if (!is_check_collision_) {
    return false;
  }
  unsigned cnt(0);
  const Vector<unsigned> new_coord(xyz_to_coord(new_xyz));
  Vector<unsigned> neighbor_coord;
  Vector<double> stride(0,0,0);
  for (int z(int(new_coord.z)-1); z <= int(new_coord.z)+1; ++z) {
    if(z < 0) {
      neighbor_coord.z = z + dim_z_;
      stride.z = -dim_z_*voxel_length_;
    } else if (z >= dim_z_) {
      neighbor_coord.z = z - dim_z_;
      stride.z = dim_z_*voxel_length_;
    } else {
      neighbor_coord.z = z;
      stride.z = 0;
    }
    for (int y(int(new_coord.y)-1); y <= int(new_coord.y)+1; ++y) {
      if(y < 0) {
        neighbor_coord.y = y + dim_y_;
        stride.y = -dim_y_*voxel_length_;
      } else if (y >= dim_y_) {
        neighbor_coord.y = y - dim_y_;
        stride.y = dim_y_*voxel_length_;
      } else {
        neighbor_coord.y = y;
        stride.y = 0;
      }
      for (int x(int(new_coord.x)-1); x <= int(new_coord.x)+1; ++x) {
        if(x < 0) {
          neighbor_coord.x = x + dim_x_;
          stride.x = -dim_x_*voxel_length_;
        } else if (x >= dim_x_) {
          neighbor_coord.x = x - dim_x_;
          stride.x = dim_x_*voxel_length_;
        } else {
          neighbor_coord.x = x;
          stride.x = 0;
        }
        const size_t neighbor_voxel_id(coord_to_voxel_id(neighbor_coord));
        const uint16_t agent_size(occupied_sizes_[neighbor_voxel_id]);
        if (agent_size != 0) {
          for (uint16_t j(0); j < agent_size; ++j) {
            if(!(voxel_id == neighbor_voxel_id && current_idx == j) &&
               new_xyz.distance(xyz_[j][neighbor_voxel_id]+stride) <
               agent_radius_*2) {
              cnt++;
            }
          }
        }
      }
    }
  }
  return cnt;
}

/*
//Although faster, cannot use the method below because 
//volume_species::positions_ and volume_species::walker::orientations_ will no 
//longer be consistent
void MicroSpace::update_species() {
  for (unsigned i(0); i != species_agents_xyz_.size(); ++i) {
    species_agents_xyz_[i].get().resize(0);
  }
  for (unsigned i(0); i < num_voxels_; ++i) {
    const uint16_t agent_size(occupied_sizes_[i]);
    if (agent_size != 0) {
      for (uint16_t j(0); j < agent_size; ++j) {
        species_agents_xyz_[global_ids_[j][i]/id_stride_].get().push_back(
                                                xyz_[j][i]);
      }
    }
  }
}
*/

void MicroSpace::update_species() {
  for (unsigned i(0); i < num_voxels_; ++i) {
    const uint16_t agent_size(occupied_sizes_[i]);
    if (agent_size != 0) {
      for (uint16_t j(0); j < agent_size; ++j) {
        unsigned species_cid;
        const unsigned agent_idx(global_id_to_agent_idx(global_ids_[j][i],
                                                      species_cid)); 
        species_agents_xyz_[species_cid].get()[agent_idx] = xyz_[j][i];
      }
    }
  }
}

unsigned MicroSpace::get_occupancy(const unsigned _species_cid) const {
  int cnt(0);
  for (size_t i(0); i < num_voxels_; ++i) {
    for (size_t j(0); j < occupied_sizes_[i]; ++j) {
      const unsigned global_id(global_ids_[j][i]);
      unsigned species_cid;
      const unsigned agent_idx(global_id_to_agent_idx(global_id, species_cid));
      if (species_cid == _species_cid) {
        ++cnt;
      }
      else {
        std::cout << "some other species:" << species_cid << " agent:" <<
          agent_idx << std::endl;
      }
    }
  }
  return cnt;
}

void MicroSpace::check_occupancy() const {
  std::cout << "in check occupancy" << std::endl;
  int cnt(0);
  for (size_t i(0); i < num_voxels_; ++i) {
    for (size_t j(0); j < pending_occupied_sizes_[i]; ++j) {
      cnt++;
      const unsigned global_id(global_ids_[j][i]);
      unsigned species_cid;
      const unsigned agent_idx(global_id_to_agent_idx(global_id, species_cid));
      if (agent_idx > 299) {
        std::cout << "agent_idx:" << agent_idx << " cnt:" << cnt << std::endl;
      }
    }
  }
  std::cout << "out check occupancy" << std::endl;
}

//can be optimized further by combining both loops below:
void MicroSpace::remove_agent(const unsigned agent_idx,
                              const unsigned back_agent_idx,
                              const unsigned species_cid) {
  const unsigned global_id(agent_idx_to_global_id(agent_idx, species_cid));
  const unsigned back_global_id(
                        agent_idx_to_global_id(back_agent_idx, species_cid));
  for (unsigned i(0); i < num_voxels_; ++i) {
    for (uint16_t j(0); j < occupied_sizes_[i]; ++j) {
      if (global_ids_[j][i] == global_id) {
        //move xyz.back() to current agent_id:
        xyz_[j][i] = xyz_[occupied_sizes_[i]-1][i];
        global_ids_[j][i] = global_ids_[occupied_sizes_[i]-1][i];
        --occupied_sizes_[i];
        pending_occupied_sizes_[i] = occupied_sizes_[i];
        break;
      }
    }
  }
  for (unsigned i(0); i < num_voxels_; ++i) {
    for (uint16_t j(0); j < occupied_sizes_[i]; ++j) {
      if (global_ids_[j][i] == back_global_id) {
        global_ids_[j][i] = global_id;
        return;
      }
    }
  }
}


void MicroSpace::walk_to_xyz(const unsigned _species_cid, 
                             const std::vector<Vector<double>>& xyz_list) {
  for (unsigned i(0); i < num_voxels_; ++i) {
    if (pending_occupied_sizes_[i] != 0) {
      for (uint16_t j(0); j < pending_occupied_sizes_[i]; ++j) {
        const unsigned global_id(global_ids_[j][i]);
        unsigned species_cid;
        const unsigned agent_idx(
                           global_id_to_agent_idx(global_id, species_cid));
        if (species_cid == _species_cid) {
          Vector<double> agent(xyz_list[agent_idx]);
          //agent.mod(dimensions_mod_);
          const size_t voxel_id(xyz_to_voxel_id(agent));
          //If we need to move the agent to another voxel:
          if (voxel_id != i) {
            ++occupied_sizes_[voxel_id];
            //get the occupied size of the target voxel id:
            const uint16_t size(occupied_sizes_[voxel_id]);
            if (xyz_.size() < size) {
              xyz_.resize(size);
              xyz_.back().resize(num_voxels_);
              global_ids_.resize(size);
              global_ids_.back().resize(num_voxels_);
            }
            xyz_[size-1][voxel_id] = agent;
            xyz_[j][i] = xyz_[occupied_sizes_[i]-1][i];
            global_ids_[size-1][voxel_id] = global_id;
            global_ids_[j][i] = global_ids_[occupied_sizes_[i]-1][i];

            if(occupied_sizes_[i] == pending_occupied_sizes_[i]) {
              --pending_occupied_sizes_[i];
              --j;
            }
            --occupied_sizes_[i]; //voxels_[i].pop_back();
            //keep already executed pending_occupied_sizes updated for the
            //next time step
            if (voxel_id < i) {
              ++pending_occupied_sizes_[voxel_id];
            }
          }
          else { //moved within the same voxel: 
            xyz_[j][i] = agent;
          }
        }
      }
    }
    //keep pending_occupied_sizes updated for the next time step
    pending_occupied_sizes_[i] = occupied_sizes_[i];
  }
}

//pending_occupied_sizes_: sizes of agents that are yet to have walked
void MicroSpace::walk() {
  for (unsigned i(0); i < num_voxels_; ++i) {
    if (pending_occupied_sizes_[i] != 0) {
      for (uint16_t j(0); j < pending_occupied_sizes_[i]; ++j) {
        Vector<double> agent(xyz_[j][i]);
        const unsigned global_id(global_ids_[j][i]);
        unsigned species_cid;
        const unsigned agent_idx(
                           global_id_to_agent_idx(global_id, species_cid));
        agent += walkers_[species_cid]->get_displacement(agent, agent_idx);
        agent.mod(dimensions_mod_);
        if(is_colliding(agent, i, j)) {
          continue;
        }
        const size_t voxel_id(xyz_to_voxel_id(agent));
        //If we need to move the agent to another voxel:
        if (voxel_id != i) {
          ++occupied_sizes_[voxel_id];
          //get the occupied size of the target voxel id:
          const uint16_t size(occupied_sizes_[voxel_id]);
          if (xyz_.size() < size) {
            xyz_.resize(size);
            xyz_.back().resize(num_voxels_);
            global_ids_.resize(size);
            global_ids_.back().resize(num_voxels_);
          }
          xyz_[size-1][voxel_id] = agent;
          xyz_[j][i] = xyz_[occupied_sizes_[i]-1][i];
          global_ids_[size-1][voxel_id] = global_id;
          global_ids_[j][i] = global_ids_[occupied_sizes_[i]-1][i];

          if(occupied_sizes_[i] == pending_occupied_sizes_[i]) {
            --pending_occupied_sizes_[i];
            --j;
          }
          --occupied_sizes_[i]; //voxels_[i].pop_back();
          //keep already executed pending_occupied_sizes updated for the
          //next time step
          if (voxel_id < i) {
            ++pending_occupied_sizes_[voxel_id];
          }
        }
        else { //moved within the same voxel: 
          xyz_[j][i] = agent;
        }
      }
    }
    //keep pending_occupied_sizes updated for the next time step
    pending_occupied_sizes_[i] = occupied_sizes_[i];
  }
  //check_overlap();
}


