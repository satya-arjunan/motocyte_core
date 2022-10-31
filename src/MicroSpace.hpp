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


#ifndef __MicroSpace_hpp
#define __MicroSpace_hpp

#include <functional>
#include <limits>
#include <Common.hpp>
#include <Space.hpp>

//minimum size for voxel edge length = agent radius*2
//if smaller, will get collisions from agents that are in 2nd 
//level neighbor voxels

class MicroSpace: public Space {
 public: 
  MicroSpace(SpaceCompartment& space_compartment,
             const Vector<double>& dimensions, const double agent_radius);
  ~MicroSpace() {}
  void populate_all();
  void walk();
  void update_species();
  std::vector<std::size_t>& get_unmoved_sizes();
  void update_unmoved_sizes();
  void set_check_collision(const bool check_collision);
  unsigned is_colliding(const Vector<double>& agent_xyz,
                        const size_t voxel_id,
                        const uint16_t current_idx) const;
  bool is_colliding(const Vector<double>& agent_xyz,
                    const size_t voxel_id) const;
  bool is_colliding(const Vector<double>& agent_xyz) const;
  void push_agents_xyz(std::vector<Vector<double> >& agents_xyz);
  bool check_collision_and_populate_agent(const Vector<double>& agent,
                                          const size_t voxel_id,
                                          const size_t agent_idx);
  void populate_agent(const Vector<double>& xyz, const size_t agent_idx,
                      const size_t species_cid);
  double get_agent_radius() const;
  void remove_agent(const unsigned agent_idx, const unsigned back_agent_idx,
                    const unsigned species_cid);
  void walk_to_xyz(const unsigned _species_cid,
                   const std::vector<Vector<double>>& xyz_list);
  unsigned get_occupancy(const unsigned _species_cid) const;
 private:
  void check_occupancy() const;
  void set_agent_cid(const unsigned species_cid);
  void check_overlap() const;
  void check_agent_overlap(const Vector<double>& agent_xyz,
                           const size_t voxel_id,
                           const size_t agent_idx) const;
  bool populate_agent(const Vector<double>&, const unsigned);
  unsigned agent_idx_to_global_id(const unsigned agent_idx,
                                  const unsigned species_cid) const; 
  unsigned global_id_to_agent_idx(const unsigned global_id,
                                  unsigned& species_cid) const;
 private:
  bool is_check_collision_ = false;
  std::vector<std::reference_wrapper<MicroSpecies>>& species_list_;
  constexpr static unsigned id_stride_ = 
    std::numeric_limits<unsigned>::max()/MAX_SPECIES;
  const double agent_radius_;
  std::vector<uint16_t> occupied_sizes_; //num_voxels_
  std::vector<uint16_t> pending_occupied_sizes_; //num_voxels_
  std::vector<std::vector<unsigned>> 
    global_ids_; //max(occupied_sizes_)*num_voxels_
  std::vector<std::vector<Vector<double>>> 
    xyz_;//max(occupied_sizes_)*num_voxels_
  std::vector<std::reference_wrapper<std::vector<Vector<double>>>>
    species_agents_xyz_; //species_size*agent_size
};

#endif /* __MicroSpace_hpp */

