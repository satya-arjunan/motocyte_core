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

#include <TrackIDSpecies.hpp>
#include <Model.hpp>


TrackIDSpecies::TrackIDSpecies(const std::string name, Model& model, 
                               TrackReader& track_reader,
                               const unsigned track_id):
  TrackSpecies(name, model, track_reader,
               track_reader.get_initial_track_id_size(track_id)),
  track_id_(track_id) {}

std::vector<Vector<double>> TrackIDSpecies::get_initial_xyz_list() {
  std::vector<unsigned>& track_ids(track_reader_.get_track_ids_at_step(0));
  std::vector<Vector<double>>& initial_xyz_list(
                         get_track_reader().get_xyz_list_at_step(0));
  track_ids_.resize(0);
  std::vector<Vector<double>> xyz_list;
  for (unsigned i(0); i < track_ids.size(); ++i) {
    if (track_ids[i] == track_id_) {
      xyz_list.push_back(initial_xyz_list[i]);
      track_ids_.push_back(track_id_);
      return xyz_list;
    }
  }
  return xyz_list;
}

void TrackIDSpecies::walk() {
  const unsigned step(get_walker().get_time()/track_reader_.get_interval());
  const std::vector<Vector<double>>& 
    next_xyz_list(track_reader_.get_xyz_list_at_step(step));
  const std::vector<unsigned>& 
    next_track_ids(track_reader_.get_track_ids_at_step(step));
  std::vector<Vector<double>> xyz_list;

  auto iter(std::find(next_track_ids.begin(), next_track_ids.end(), track_id_));
  if (iter == next_track_ids.end()) {
    if (track_ids_.size()) {
      remove_agent(0);
      track_ids_.resize(0);
    }
  }
  else {
    if (!track_ids_.size()) {
      add_agent(next_xyz_list[iter-next_track_ids.begin()]);
      track_ids_.push_back(track_id_);
    }
    xyz_list.push_back(next_xyz_list[iter-next_track_ids.begin()]);
  }

  /*
  std::cout << "\nstep:" << step << std::endl;
  std::cout << "size:" << xyz_list.size() << " occupancy:" <<
    get_space_compartment().get_micro_space().get_occupancy(get_species_cid())
    << std::endl;
  for (unsigned i(0); i < track_ids_.size(); ++i) {
    for (unsigned j(0); j < next_track_ids.size(); ++j) {
      if (track_ids_[i] == next_track_ids[j]) {
        std::cout << "i:" << i << std::endl;
        std::cout << xyz_list[i].x << " " << xyz_list[i].y << " " << 
          xyz_list[i].z << std::endl;
        std::cout << next_xyz_list[j].x << " " << next_xyz_list[j].y << " " << 
          next_xyz_list[j].z << std::endl;
      }
    }
  }
  */
  get_space_compartment().get_micro_space().walk_to_xyz(get_species_cid(),
                                                        xyz_list);
  /*
  xyz_list = get_xyz_list();
  std::cout << "\nwalked:" << step << std::endl;
  std::cout << "size:" << xyz_list.size() << " occupancy:" <<
    get_space_compartment().get_micro_space().get_occupancy(get_species_cid())
    << std::endl;
  for (unsigned i(0); i < track_ids_.size(); ++i) {
    for (unsigned j(0); j < next_track_ids.size(); ++j) {
      if (track_ids_[i] == next_track_ids[j]) {
        std::cout << "i:" << i << std::endl;
        std::cout << xyz_list[i].x << " " << xyz_list[i].y << " " << 
          xyz_list[i].z << std::endl;
        std::cout << next_xyz_list[j].x << " " << next_xyz_list[j].y << " " << 
          next_xyz_list[j].z << std::endl;
      }
    }
  }
  //exit(0);
  */
  MicroSpecies::walk();
}





