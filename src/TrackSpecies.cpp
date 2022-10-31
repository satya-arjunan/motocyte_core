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

#include <TrackSpecies.hpp>
#include <Model.hpp>


TrackSpecies::TrackSpecies(const std::string name, Model& model, 
                           TrackReader& track_reader, const int size):
  VolumeSpecies(name, model, size),
  track_reader_(track_reader) {
    //set higher priority than other walkers to reduce volume
    //exclusion/collision:
    get_walker().set_priority(10);
    get_walker().set_interval(track_reader_.get_interval());
    track_ids_ = track_reader_.get_track_ids_at_step(0);
  }

TrackSpecies::TrackSpecies(const std::string name, Model& model, 
                           TrackReader& track_reader):
  TrackSpecies(name, model, track_reader,
               track_reader.get_initial_track_ids_size()) {}

std::vector<Vector<double>> TrackSpecies::get_initial_xyz_list() {
  return track_reader_.get_xyz_list_at_step(0);
}

TrackReader& TrackSpecies::get_track_reader() {
  return track_reader_;
}

std::vector<unsigned>& TrackSpecies::get_track_ids() {
  return track_ids_;
}

void TrackSpecies::walk() {
  const unsigned step(get_walker().get_time()/track_reader_.get_interval());
  const std::vector<Vector<double>>& 
    next_xyz_list(track_reader_.get_xyz_list_at_step(step));
  const std::vector<unsigned>& 
    next_track_ids(track_reader_.get_track_ids_at_step(step));
  std::vector<Vector<double>> xyz_list(next_xyz_list.size());
  for (unsigned i(0); i < track_ids_.size(); ++i) { 
    const unsigned track_id(track_ids_[i]);
    auto iter(std::find(next_track_ids.begin(), next_track_ids.end(),
                        track_id));
    if (iter == next_track_ids.end()) {
      remove_agent(i);
      track_ids_[i] = track_ids_.back();
      track_ids_.pop_back();
      --i;
    }
    else {
      xyz_list[i] = next_xyz_list[iter-next_track_ids.begin()];
    }
  }
  for (unsigned i(0); i < next_track_ids.size(); ++i) { 
    const unsigned track_id(next_track_ids[i]);
    auto iter(std::find(track_ids_.begin(), track_ids_.end(), track_id));
    if (iter == track_ids_.end()) {
      add_agent(next_xyz_list[i]);
      track_ids_.push_back(track_id);
      xyz_list[size()-1] = next_xyz_list[i];
    }
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





