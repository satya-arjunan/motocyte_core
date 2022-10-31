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

#include <sstream>
#include <climits>
#include <algorithm>
#include <TrackReader.hpp>

TrackReader::TrackReader(std::string file) {
  std::stringstream filename(file);
  std::vector<std::string> split_strings;
  std::string str; 
  while (std::getline(filename, str, '_')) {
    split_strings.push_back(str);
  }
  interval_ = std::stof(split_strings[split_strings.size()-2]);
  track_file_.open(filename.str().c_str(), std::ios::in);
  parse_file();
}

void TrackReader::parse_file() {
  const unsigned skiprows(1);
  unsigned cnt(0);
  std::string line;
  while (cnt < skiprows && std::getline(track_file_, line)) { 
    ++cnt;
  }

  while (std::getline(track_file_, line)) { 
    std::stringstream str(line); 
    std::string word;
    std::vector<std::string> row;
    while (std::getline(str, word, ',')) { 
      row.push_back(word); 
    } 
    const Vector<double> xyz(
          std::stof(row[0]), std::stof(row[1]), std::stof(row[2]));
    const unsigned step(std::stoi(row[6]));
    const unsigned agent_id(std::stoi(row[7]));
    max_.x = std::max(max_.x, xyz.x);
    max_.y = std::max(max_.y, xyz.y);
    max_.z = std::max(max_.z, xyz.z);
    min_.x = std::min(min_.x, xyz.x);
    min_.y = std::min(min_.y, xyz.y);
    min_.z = std::min(min_.z, xyz.z);
    if (step_track_ids_.size() < step) { 
      step_track_ids_.resize(step);
      step_xyz_list_.resize(step);
    }
    //step starts at 1, not 0:
    step_track_ids_[step-1].push_back(agent_id);
    step_xyz_list_[step-1].push_back(xyz);
  } 
  dimensions_ = Vector<double>(max_.x-min_.x+10, max_.y-min_.y+10,
                              max_.z-min_.z+10); 

  for (unsigned i(0); i < step_track_ids_.size(); ++i) {
    for (unsigned j(0); j < step_track_ids_[i].size(); ++j) {
      //make chronological agent_ids
      const unsigned id(step_track_ids_[i][j]);
      if (agent_id_map_.find(id) == agent_id_map_.end()) {
        agent_id_map_[id] = agent_id_cnt_;
        ++agent_id_cnt_;
      }
      step_track_ids_[i][j] = agent_id_map_[id];
      //translate the points to set min_ as origin:
      //step_xyz_list_[i][j] -= min_;
      step_xyz_list_[i][j].x += 2.5;
      step_xyz_list_[i][j].y += 2.5;
      step_xyz_list_[i][j].z += 2.5;
    }
  }

  std::cout << step_xyz_list_[0][0].x << " " << step_xyz_list_[0][0].y << " " 
    << step_xyz_list_[0][0].z << std::endl;
  std::cout << min_.x << " " << min_.y << " " << min_.z << std::endl;
}

Vector<double>& TrackReader::get_dimensions() {
  std::cout << "dimensions:" << dimensions_.x << " " << dimensions_.y << " " 
    << dimensions_.z << std::endl;
  return dimensions_;
}

unsigned TrackReader::get_size() {
  return agent_id_cnt_;
}

double TrackReader::get_interval() const {
  return interval_;
} 

double TrackReader::get_end_time() const {
  return (step_track_ids_.size()-1)*interval_;
} 

unsigned TrackReader::get_mean_agents_per_step() const {
  unsigned size(0);
  for (unsigned i(0); i < step_track_ids_.size(); ++i) {
    size += step_track_ids_[i].size();
  }
  return std::round(double(size)/step_track_ids_.size());
}

unsigned TrackReader::get_initial_track_id_size(const int track_id) {
  std::vector<unsigned>& track_ids(get_track_ids_at_step(0));
  if (std::find(track_ids.begin(), track_ids.end(), track_id) ==
      track_ids.end()) { 
    return 0;
  }
  return 1;
}

unsigned TrackReader::get_initial_track_ids_size() {
  return get_track_ids_at_step(0).size();
}

std::vector<unsigned>& TrackReader::get_track_ids_at_step(const unsigned step) {
  return step_track_ids_[step];
}

std::vector<Vector<double>>& TrackReader::get_xyz_list_at_step(
                                                        const unsigned step) {
  return step_xyz_list_[step];
}

