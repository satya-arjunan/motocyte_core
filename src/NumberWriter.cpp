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
#include <Motocyte.hpp>
#include <NumberWriter.hpp>
#include <Model.hpp>
#include <Stepper.hpp>


NumberWriter::NumberWriter(Model& model):
  Process("NumberWriter", model, std::numeric_limits<double>::infinity()),
  file_name_("NumberLog.csv"),
  space_compartment_(model.get_space_compartment()) {
    set_queue_id(get_stepper().get_process_queue().push(this));
    set_priority(-5); //lowest priority in process queue
  }

void NumberWriter::set_file_name(std::string file_name) {
  file_name_ = file_name;
}

double NumberWriter::step() {
  log_species();
  log_file_.flush();
  return Process::step();
}

void NumberWriter::initialize() {
  std::ostringstream fileName;
  fileName << file_name_ << std::ends;
  log_file_.open(fileName.str().c_str(), std::ios::out | std::ios::trunc);
  initialize_log();
  log_species();
  log_file_.flush();
}

std::ofstream& NumberWriter::get_log_file() {
  return log_file_;
}

void NumberWriter::add(MicroSpecies& micro_species) {
  species_.push_back(reinterpret_cast<Species*>(&micro_species));
}

void NumberWriter::add(MesoSpecies& meso_species) {
  species_.push_back(reinterpret_cast<Species*>(&meso_species));
}

void NumberWriter::initialize_log() { 
  log_file_ << "Time";
  for(unsigned i(0); i != species_.size(); ++i) {
    log_file_ << "," << species_[i]->get_name();
  }
  log_file_ << std::endl;
}

void NumberWriter::log_species() {
  const double currentTime(get_stepper().get_time());
  log_file_ << currentTime;
  for(unsigned i(0); i < species_.size(); ++i) {
    MicroSpecies* micro_species(dynamic_cast<MicroSpecies*>(species_[i]));
    if (micro_species) {
      VolumeSpecies* volume_species(
                        dynamic_cast<VolumeSpecies*>(micro_species));
      if (volume_species) {
        log(*volume_species, i);
      }
      else {
        PointSpecies* point_species(dynamic_cast<PointSpecies*>(micro_species));
        if (point_species) {
          log(*point_species, i);
        }
      }
    }
    else {
      MesoSpecies* meso_species(dynamic_cast<MesoSpecies*>(species_[i]));
      if (meso_species) {
        log(*meso_species, i);
      }
      else {
        std::cout << "Unknown species type in number writer" << std::endl;
      }
    }
  }
  log_file_ << std::endl;
}

void NumberWriter::log(VolumeSpecies& volume_species, const unsigned index) {
  const std::vector<Vector<double>>& 
    agents(volume_species.get_relative_positions());
  log_file_ << "," << agents.size();
}  

void NumberWriter::log(PointSpecies& point_species, const unsigned index) {
  const std::vector<unsigned>& 
    agents_sizes(point_species.get_compartment_points_sizes());
  unsigned size(0);
  for (unsigned i(0); i < agents_sizes.size(); ++i) {
    size += agents_sizes[i];
  }
  log_file_ << "," << size;
}

void NumberWriter::log(MesoSpecies& meso_species, const unsigned index) {
  const std::vector<unsigned>& sizes(meso_species.get_sizes());
  const std::vector<unsigned>& 
    occupied_voxels(meso_species.get_occupied_voxels());
  int size(occupied_voxels.size());
  unsigned factor_size(0);
  for (int i(0); i < size; ++i) {
    const unsigned voxel_id(occupied_voxels[i]);
    factor_size += sizes[voxel_id];
  }
  log_file_ << "," << size;
}

