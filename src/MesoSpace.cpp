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

#include <cmath>
#include <cstring>
#include <random>
#include <algorithm>
#include <Motocyte.hpp>
#include <SpaceCompartment.hpp>
#include <Model.hpp>
#include <Walker.hpp>
#include <MesoSpace.hpp>

MesoSpace::MesoSpace(SpaceCompartment& space_compartment,
                     const Vector<double>& dimensions, const double agent_radius):
  Space(space_compartment, dimensions, agent_radius, 0, 4),
  coord_offsets_({Vector<int>(-1,0,0),
                  Vector<int>(0,-1,0),
                  Vector<int>(0,0,-1),
                  Vector<int>(1,0,0),
                  Vector<int>(0,1,0),
                  Vector<int>(0,0,1)}),
  walk_dist_(0,5),
  drift_walk_dist_(0,2) {
    rect_mesh_.resize(3);
    for (unsigned i(0); i < dim_x_; ++i) {
      rect_mesh_[0].push_back(i*voxel_length_);
    }
    for (unsigned i(0); i < dim_y_; ++i) {
      rect_mesh_[1].push_back(i*voxel_length_);
    }
    for (unsigned i(0); i < dim_z_; ++i) {
      rect_mesh_[2].push_back(i*voxel_length_);
    }
  }


void MesoSpace::push_factors(std::vector<unsigned>& factors) {
  factors.resize(num_voxels_, 0);
  species_factors_.push_back(factors);
  curr_occupied_voxels_.push_back(new std::vector<unsigned>());
  next_occupied_voxels_.push_back(new std::vector<unsigned>());
  removed_voxel_ids_.resize(removed_voxel_ids_.size()+1);
  all_voxels_.resize(species_factors_.size());
  all_voxels_.back().resize(num_voxels_, 0);
}

/*
void MesoSpace::populate_all() {
  const std::vector<std::reference_wrapper<MesoSpecies>>& 
    species_list(space_compartment_.get_meso_species_list());

  std::list<unsigned> voxel_ids(num_voxels_);
  std::iota(voxel_ids.begin(), voxel_ids.end(), 0);
  std::vector<unsigned> sampled_voxel_ids;
  std::sample(voxel_ids.begin(), voxel_ids.end(), 
      std::back_inserter(sampled_voxel_ids), species_list.size(),
      get_model().get_rng());
  std::shuffle(sampled_voxel_ids.begin(), sampled_voxel_ids.end(), 
               get_model().get_rng());

  for (unsigned species_cid(0); species_cid != species_list.size();
       ++species_cid) {
    const unsigned init_size(species_list[species_cid].get().get_init_size());
    if (init_size) {
      const unsigned voxel_id(sampled_voxel_ids[species_cid]);
      species_factors_[species_cid].get()[voxel_id] = init_size;
      all_voxels_[species_cid][voxel_id] = 1;
      (*curr_occupied_voxels_)[species_cid].push_back(voxel_id);
    }
  }
  std::cout << "Populated factors" << std::endl;
}
*/


void MesoSpace::populate_all() {
  const std::vector<std::reference_wrapper<MesoSpecies>>& 
    species_list(space_compartment_.get_meso_species_list());
  for (unsigned species_cid(0); species_cid != species_list.size();
       ++species_cid) {
    species_list[species_cid].get().populate();
    if (!species_list[species_cid].get().get_is_populated()) {
      populate_species(species_list[species_cid].get());
    }
  }
}

void MesoSpace::populate_species(MesoSpecies& species) {
  const unsigned init_size(species.get_init_size());
  const unsigned species_cid(species.get_species_cid());
  if (init_size) {
    for (unsigned voxel_id(0); voxel_id < num_voxels_; ++voxel_id) {
      species_factors_[species_cid].get()[voxel_id] = init_size;
      all_voxels_[species_cid][voxel_id] = 1;
      (*curr_occupied_voxels_[species_cid]).push_back(voxel_id);
    }
    species.set_size(num_voxels_*init_size);
  }
}

void MesoSpace::add_factor(const unsigned species_cid, MesoSpecies& species,
                           Vector<double>& xyz) {
  const unsigned voxel_id(xyz_to_voxel_id(xyz));
  if (!species_factors_[species_cid].get()[voxel_id]) {
    all_voxels_[species_cid][voxel_id] = 1;
    (*curr_occupied_voxels_[species_cid]).push_back(voxel_id);
  }
  ++species_factors_[species_cid].get()[voxel_id];
  ++species.get_size();
}

void MesoSpace::add_factor(const unsigned species_cid, MesoSpecies& species,
                           Vector<double>& xyz, const unsigned max) {
  const unsigned voxel_id(xyz_to_voxel_id(xyz));
  if (!species_factors_[species_cid].get()[voxel_id] && max > 0) {
    all_voxels_[species_cid][voxel_id] = 1;
    (*curr_occupied_voxels_[species_cid]).push_back(voxel_id);
  }
  if (species_factors_[species_cid].get()[voxel_id] < max) {
    ++species_factors_[species_cid].get()[voxel_id];
    ++species.get_size();
  }
}

void MesoSpace::replace_factors(std::vector<unsigned>& src_species_cids,
                                std::vector<MesoSpecies*>& src_species_list,
                                const unsigned tar_species_cid,
                                MesoSpecies& tar_species, Vector<double>& xyz,
                                const unsigned ignore_index) {
  const unsigned voxel_id(xyz_to_voxel_id(xyz));
  for (unsigned i(0); i < src_species_cids.size() && i != ignore_index; ++i) {
    unsigned src_species_cid(src_species_cids[i]);
    MesoSpecies& src_species(*src_species_list[i]);
    std::vector<unsigned>& factors(species_factors_[src_species_cid].get()); 
    unsigned remove_size(factors[voxel_id]);
    if (remove_size) {
      factors[voxel_id] = 0;
      src_species.set_size(src_species.get_size()-remove_size);
      all_voxels_[src_species_cid][voxel_id] = 0;
      removed_voxel_ids_[src_species_cid].push_back(voxel_id);
      if (!species_factors_[tar_species_cid].get()[voxel_id]) {
        all_voxels_[tar_species_cid][voxel_id] = 1;
        (*curr_occupied_voxels_[tar_species_cid]).push_back(voxel_id);
      }
      species_factors_[tar_species_cid].get()[voxel_id] += remove_size;
      tar_species.set_size(tar_species.get_size()+remove_size);
    }
  }
}

//must call update_curr_occupied_voxels after calling this function
void MesoSpace::remove_factors(const unsigned species_cid, MesoSpecies& species,
                              Vector<double>& xyz, unsigned remove_size) {
  const unsigned voxel_id(xyz_to_voxel_id(xyz));
  std::vector<unsigned>& factors(species_factors_[species_cid].get()); 
  if (factors[voxel_id] >= remove_size) {
    factors[voxel_id] -= remove_size;
    species.set_size(species.get_size()-remove_size);
  }
  if (!factors[voxel_id]) {
    all_voxels_[species_cid][voxel_id] = 0;
    removed_voxel_ids_[species_cid].push_back(voxel_id);
  }
}

void MesoSpace::update_curr_occupied_voxels(const unsigned species_cid) {
  std::vector<unsigned>& removed_voxel_ids(removed_voxel_ids_[species_cid]);
  if (removed_voxel_ids.size()) {
    std::vector<unsigned>& 
      curr_occupied_voxels((*curr_occupied_voxels_[species_cid]));
    std::vector<uint8_t>& all_voxels(all_voxels_[species_cid]);
    if (curr_occupied_voxels.size()*removed_voxel_ids.size() <
        all_voxels.size()) {
      for (unsigned i(0); i < curr_occupied_voxels.size() &&
           removed_voxel_ids.size(); ++i) {
        const unsigned voxel_id(curr_occupied_voxels[i]);
        auto iter(std::find(removed_voxel_ids.begin(), removed_voxel_ids.end(),
                            voxel_id));
        if (iter != removed_voxel_ids.end()) {
          *iter = removed_voxel_ids.back();
          removed_voxel_ids.pop_back();
          curr_occupied_voxels[i] = curr_occupied_voxels.back();
          curr_occupied_voxels.pop_back();
          --i;
        }
      }
    }
    else {
      removed_voxel_ids.resize(0);
      curr_occupied_voxels.resize(0);
      std::vector<unsigned>& factors(species_factors_[species_cid].get());
      for (unsigned voxel_id(0); voxel_id < all_voxels.size(); ++voxel_id) {
        if (all_voxels[voxel_id]) {
          curr_occupied_voxels.push_back(voxel_id);
          if (!factors[voxel_id]) {
            std::cout << "error in factor size in update curr occupied " <<
              "voxels:" << voxel_id << std::endl;
          }
        }
      }
    }
  }
}

void MesoSpace::remove_factors(const unsigned species_cid,
                               MesoSpecies& species,
                               unsigned remove_size) {
  const unsigned current_size(species.get_size());
  if (current_size < remove_size) {
    remove_size = current_size;
  }
  std::list<unsigned> factor_ids(current_size);
  std::iota(factor_ids.begin(), factor_ids.end(), 0);
  std::vector<unsigned> sampled_ids;
  std::sample(factor_ids.begin(), factor_ids.end(), 
      std::back_inserter(sampled_ids), remove_size, get_model().get_rng());

  unsigned global_cnt(0);
  unsigned cnt(0);
  unsigned factor_id(sampled_ids[cnt]);
  std::vector<unsigned>& factors(species_factors_[species_cid].get()); 
  std::vector<unsigned>& 
    curr_occupied_voxels((*curr_occupied_voxels_[species_cid]));
  std::vector<uint8_t>& all_voxels(all_voxels_[species_cid]);
  for (unsigned i(0);
       i < curr_occupied_voxels.size() && cnt < sampled_ids.size(); ++i) {
    const unsigned voxel_id(curr_occupied_voxels[i]);
    const unsigned factor_size(factors[voxel_id]);
    while (global_cnt+factor_size > factor_id && cnt < sampled_ids.size()) {
      --factors[voxel_id];
      ++cnt;
      factor_id = sampled_ids[cnt];
    }
    global_cnt += factor_size;
    if (!factors[voxel_id]) {
      all_voxels[voxel_id] = 0;
      curr_occupied_voxels[i] = curr_occupied_voxels.back();
      curr_occupied_voxels.pop_back();
      --i;
    }
  }
  species.set_size(current_size-remove_size);
}

unsigned MesoSpace::get_position_factor(const unsigned species_cid,
                                        MesoSpecies& species,
                                        Vector<double>& xyz) {
  const unsigned voxel_id(xyz_to_voxel_id(xyz));
  return species_factors_[species_cid].get()[voxel_id];
}



std::vector<unsigned> MesoSpace::get_position_factors(const unsigned species_cid,
                                        std::vector<Vector<double>> positions,
                                        std::vector<unsigned>& agent_factors,
                                        const unsigned point_size) {
  std::vector<unsigned>& factors(species_factors_[species_cid].get());
  std::vector<unsigned> position_factors;
  for (unsigned i(0); i < agent_factors.size(); ++i) {
    for (unsigned j(0); j < point_size; ++j) {
      const unsigned k(i*point_size+j);
      const unsigned factor(factors[xyz_to_voxel_id(positions[k])]);
      agent_factors[i] += factor;
      position_factors.push_back(factor);
    }
  }
  return position_factors;
}

const std::vector<std::vector<double>>& MesoSpace::get_rect_mesh() const {
  return rect_mesh_;
}

std::vector<unsigned>& MesoSpace::get_species_occupied_voxels(
                                                  const unsigned species_cid) {

  return (*curr_occupied_voxels_[species_cid]);
}


unsigned MesoSpace::direction_to_voxel_id(const unsigned direction,
                                          const unsigned voxel_id) const {
  Vector<unsigned> coord(voxel_id_to_coord(voxel_id));
  Vector<int> offset(coord.x+coord_offsets_[direction].x,
                     coord.y+coord_offsets_[direction].y,
                     coord.z+coord_offsets_[direction].z);
  Vector<unsigned> dest;
  if(offset.z < 0) {
    dest.z = offset.z + dim_z_;
  } else if (offset.z >= dim_z_) {
    dest.z = offset.z - dim_z_;
  } else {
    dest.z = offset.z;
  }
  if(offset.y < 0) {
    dest.y = offset.y + dim_y_;
  } else if (offset.y >= dim_y_) {
    dest.y = offset.y - dim_y_;
  } else {
    dest.y = offset.y;
  }
  if(offset.x < 0) {
    dest.x = offset.x + dim_x_;
  } else if (offset.x >= dim_x_) {
    dest.x = offset.x - dim_x_;
  } else {
    dest.x = offset.x;
  }
  return coord_to_voxel_id(dest);
}

void MesoSpace::check_occupied_voxels() {
  const unsigned species_size(species_factors_.size());
  for (unsigned species_cid(0); species_cid < species_size; ++species_cid) { 
    std::vector<uint8_t>& all_voxels(all_voxels_[species_cid]);
    std::vector<unsigned>& 
      curr_occupied_voxels((*curr_occupied_voxels_[species_cid]));
    unsigned cnt1(0);
    for (unsigned i(0); i != all_voxels.size(); ++i) {
      if (all_voxels[i]) {
        ++cnt1;
      }
    }
    unsigned cnt2(0);
    for (unsigned i(0); i < curr_occupied_voxels.size(); ++i) {
      if (!all_voxels[curr_occupied_voxels[i]]) {
        std::cout << "error in all voxel" << std::endl;
      }
      ++cnt2;
    }
    if (cnt1 != cnt2) {
      std::cout << "not in sync:" << cnt1 << " " << cnt2 << " " <<
        get_species(species_cid).get_name() << std::endl;
    }
  }
}

//walk probability is 1
//for each factor, need to only decide the random direction, and move it
//should reduce the walk probability to < 1 for better accuracy when step
//intervals are large
void MesoSpace::walk(const unsigned species_cid) {
  //check_occupied_voxels();
  std::vector<unsigned>& 
    curr_occupied_voxels((*curr_occupied_voxels_[species_cid]));
  std::vector<unsigned>& 
    next_occupied_voxels((*next_occupied_voxels_[species_cid]));
  std::vector<uint8_t>& all_voxels(all_voxels_[species_cid]);
  //set all voxels occupied by current species to be vacant. This
  //will be used as a flag to push a voxel in next_occupied_voxels only once.
  for (unsigned i(0); i < curr_occupied_voxels.size(); ++i) {
    all_voxels[curr_occupied_voxels[i]] = 0;
  }
  next_occupied_voxels.resize(0);
  //species_factors_ contains the size of the species factor in each voxel
  std::vector<unsigned>& factors(species_factors_[species_cid].get());
  for (unsigned i(0); i < curr_occupied_voxels.size(); ++i) {
    const unsigned voxel_id(curr_occupied_voxels[i]);
    const unsigned factor_size(factors[voxel_id]);
    if (!factor_size) {
      std::cout << "error in factor size" << std::endl;
    }
    //walk each factor of current species in current voxel to its neighbor
    const unsigned quotient(factor_size/6);
    if (quotient) {
      const unsigned remainder(factor_size%6);
      std::vector<unsigned> dest_sizes(6, quotient);
      for (unsigned i(0); i < remainder; ++i) {
        dest_sizes[walk_dist_(get_model().get_rng())]++;
      }
      for (unsigned direction(0); direction < 6; ++direction) {
        const unsigned dest_voxel_id(
                               direction_to_voxel_id(direction, voxel_id)); 
        if (!factors[dest_voxel_id]) {
          //indicate that the dest_voxel of the species is occupied
          //(1 is occupied, 0 is vacant)
          all_voxels[dest_voxel_id] = 1;
          next_occupied_voxels.push_back(dest_voxel_id);
        }
        factors[dest_voxel_id] += dest_sizes[direction];
      }
    } else {
      for (unsigned i(0); i < factor_size; ++i) {
        const unsigned direction(walk_dist_(get_model().get_rng()));
        const unsigned dest_voxel_id(
                                 direction_to_voxel_id(direction, voxel_id));
        //if at present the dest_voxel does not have any factor of the
        //current species:
        if (!factors[dest_voxel_id]) {
          //indicate that the dest_voxel of the species is occupied
          //(1 is occupied, 0 is vacant)
          all_voxels[dest_voxel_id] = 1;
          next_occupied_voxels.push_back(dest_voxel_id);
        }
        ++factors[dest_voxel_id];
      }
    }
    factors[voxel_id] = 0;
  }
  std::vector<unsigned>* tmp(curr_occupied_voxels_[species_cid]);
  curr_occupied_voxels_[species_cid] = next_occupied_voxels_[species_cid];
  next_occupied_voxels_[species_cid] = tmp;
  //check_occupied_voxels();
}

//walk probability is 1
//for each factor, need to only decide the random direction, and move it
//should reduce the walk probability to < 1 for better accuracy when step
//intervals are large
void MesoSpace::walk() {
  //check_occupied_voxels();
  const unsigned species_size(species_factors_.size());
  for (unsigned species_cid(0); species_cid < species_size; ++species_cid) { 
    walk(species_cid);
  }
  //check_occupied_voxels();
}


//walk probability is 1
//for each factor, need to only decide the random direction, and move it
void MesoSpace::walk_individually() {
  //check_occupied_voxels();
  const unsigned species_size(species_factors_.size());
  for (unsigned species_cid(0); species_cid < species_size; ++species_cid) { 
    std::vector<unsigned>& 
      curr_occupied_voxels((*curr_occupied_voxels_[species_cid]));
    std::vector<unsigned>& 
      next_occupied_voxels((*next_occupied_voxels_[species_cid]));
    std::vector<uint8_t>& all_voxels(all_voxels_[species_cid]);
    //set all voxels occupied by current species to be vacant. This
    //will be used as a flag to push a voxel in next_occupied_voxels only once.
    for (unsigned i(0); i < curr_occupied_voxels.size(); ++i) {
      all_voxels[curr_occupied_voxels[i]] = 0;
    }
    next_occupied_voxels.resize(0);
    //species_factors_ contains the size of the species factor in each voxel
    std::vector<unsigned>& factors(species_factors_[species_cid].get());
    for (unsigned i(0); i < curr_occupied_voxels.size(); ++i) {
      const unsigned voxel_id(curr_occupied_voxels[i]);
      const unsigned factor_size(factors[voxel_id]);
      if (!factor_size) {
        std::cout << "error in factor size" << std::endl;
      }
      //walk each factor of current species in current voxel to its neighbor
      for (unsigned i(0); i < factor_size; ++i) {
        const unsigned direction(walk_dist_(get_model().get_rng()));
        const unsigned dest_voxel_id(
                                 direction_to_voxel_id(direction, voxel_id));
        --factors[voxel_id];
        //if at present the dest_voxel does not have any factor of the
        //current species:
        if (!factors[dest_voxel_id]) {
          //indicate that the dest_voxel of the species is occupied
          //(1 is occupied, 0 is vacant)
          all_voxels[dest_voxel_id] = 1;
          next_occupied_voxels.push_back(dest_voxel_id);
        }
        ++factors[dest_voxel_id];
      }
      //If current voxel is not empty and it has not been pushed into 
      //next_occupied_voxels, push it:
      if (factors[voxel_id] && !all_voxels[voxel_id]) {
        all_voxels[voxel_id] = 1;
        next_occupied_voxels.push_back(voxel_id);
      }
    }
  std::vector<unsigned>* tmp(curr_occupied_voxels_[species_cid]);
  curr_occupied_voxels_[species_cid] = next_occupied_voxels_[species_cid];
  next_occupied_voxels_[species_cid] = tmp;
  }
  //check_occupied_voxels();
}

/*
void MesoSpace::walk() {
  const unsigned species_size(species_factors_.size());
  for (unsigned species_cid(0); species_cid < species_size; ++species_cid) { 
    std::vector<unsigned>& factors(species_factors_[species_cid].get());
    for (unsigned voxel_id(0); voxel_id < num_voxels_; ++voxel_id) {
      const unsigned factor_size(factors[voxel_id]);
      if (factor_size) {
        for (unsigned i(0); i < factor_size; ++i) {
          unsigned direction(walk_dist_(get_model().get_rng()));
          unsigned dest_voxel_id(direction_to_voxel_id(direction, voxel_id));
          --factors[voxel_id];
          ++factors[dest_voxel_id];
        }
      }
    }
  }
}
*/

void MesoSpace::walk_drift() {
  const unsigned species_size(species_factors_.size());
  for (unsigned species_cid(0); species_cid < species_size; ++species_cid) { 
    std::vector<unsigned>& factors(species_factors_[species_cid].get());
    for (unsigned voxel_id(0); voxel_id < num_voxels_; ++voxel_id) {
      const unsigned factor_size(factors[voxel_id]);
      if (factor_size) {
        for (unsigned i(0); i < factor_size; ++i) {
          unsigned direction(drift_walk_dist_(get_model().get_rng()));
          unsigned dest_voxel_id(direction_to_voxel_id(direction, voxel_id));
          --factors[voxel_id];
          ++factors[dest_voxel_id];
        }
      }
    }
  }
}
