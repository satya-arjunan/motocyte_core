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
#include <iomanip>
#include <iostream>
#include <Motocyte.hpp>
#include <TrackWriter.hpp>
#include <Model.hpp>
#include <Stepper.hpp>


TrackWriter::TrackWriter(Model& model):
  NumberWriter(model),
  volume_species_(NULL) {
    set_file_name("_Position.csv");
    set_priority(-5); //lowest priority in process queue
  }

void TrackWriter::add(MesoSpecies& meso_species) {
  std::cout << "TrackWriter: Only VolumeSpecies is supported." <<
    " MesoSpecies is not supported:" << meso_species.get_name() << std::endl;
  exit(1);
}

void TrackWriter::add(MicroSpecies& micro_species) {
  VolumeSpecies* volume_species(dynamic_cast<VolumeSpecies*>(&micro_species));
  if (!volume_species) {
    std::cout << "TrackWriter: Only VolumeSpecies is supported." <<
      " This MicroSpecies is not supported:" << micro_species.get_name() <<
      std::endl;
    exit(1);
  }
  volume_species_ = volume_species;
  const double interval(volume_species_->get_walker().get_interval());
  std::ostringstream filename;
  filename << interval << "_Position.csv";
  set_file_name(filename.str());
}

void TrackWriter::set_padding(const double padding) {
  padding_ = padding;
}
void TrackWriter::set_skip_frames(const unsigned skip_frames) {
  skip_frames_ = skip_frames;
}

void TrackWriter::initialize() {
  if (!volume_species_) {
    return;
  }
  NumberWriter::initialize();
}

void TrackWriter::initialize_log() { 
  rect_max_ = 
    volume_species_->get_space_compartment().get_micro_space().get_dimensions();
  rect_max_ -= padding_;
  rect_min_ = Vector<double>(padding_, padding_, padding_);

  const std::vector<unsigned>& agent_cids(volume_species_->get_agent_cids());
  agent_cids_.resize(agent_cids.size());
  active_cids_.resize(agent_cids.size(), 0);

  get_log_file() << "Position X,Position Y,Position Z,Unit,Category," <<
    "Collection,Time,TrackID,ID," << std::endl;
}

void TrackWriter::log_species() {
  ++frame_cnt_;
  if (frame_cnt_ < skip_frames_) {
    return;
  }
  const std::vector<Vector<double>>&
    agents(volume_species_->get_relative_positions());
  const unsigned size(agents.size());
  if (size > 0) {
    for(unsigned i(0); i < size; ++i) {
      const Vector<double>& agent(agents[i]);
      if (agent.x >= rect_min_.x && agent.x <= rect_max_.x &&
          agent.y >= rect_min_.y && agent.y <= rect_max_.y &&
          agent.z >= rect_min_.z && agent.z <= rect_max_.z) {
        if (!active_cids_[i]) {
          agent_cids_[i] = ++id_cnt_; 
          active_cids_[i] = 1;
        }
        Vector<double> u;
        if (prev_agents_.size()) {
          const Vector<double>& p(prev_agents_[i]);
          u = agent-p;
          u = u.norm();
        }
        get_log_file() << std::setprecision(20) <<  std::scientific <<
          agent.x << "," <<
          agent.y << "," <<
          agent.z << "," << std::resetiosflags(std::ios::fixed |
                                               std::ios::scientific ) <<
          "um,Surface,Position," <<
          int(get_stepper().get_time()/get_interval()) << "," <<
          1000000000+agent_cids_[i] << "," << agent_cids_[i] << "," <<
          std::endl;
      }
      else {
        active_cids_[i] = 0;
      }
    }
  }
  prev_agents_ = agents;
}


