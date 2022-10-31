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

#include <stdexcept>
#include <math.h>
#include <sys/time.h>
#include <Motocyte.hpp>
#include <Model.hpp>
#include <MotileSpecies.hpp>

Model::Model(const Vector<double>& dimensions, const double agent_radius,
             unsigned long seed):
  rng_(rd_()),
  gsl_rng_(gsl_rng_alloc (gsl_rng_mt19937)),
  visual_writer_(*this),
  number_writer_(*this),
  track_writer_(*this),
  space_compartment_("r", *this, dimensions, agent_radius) {
    if (!seed) {
      struct timeval tv;
      gettimeofday(&tv, 0);
      seed = tv.tv_sec + tv.tv_usec;
    }
    gsl_rng_set(gsl_rng_, seed);
  }
      
void Model::initialize() {
  get_space_compartment().populate_all();
  initialize_processes();
  get_stepper().initialize();
}

void Model::initialize_processes() {
  for (unsigned i(0); i < process_list_.size(); ++i) {
    process_list_[i]->initialize();
  }
}

void Model::run(const double interval) {
  const double prev_time(stepper_.get_time());
  const double end_time(prev_time+interval);
  while (stepper_.step() <= end_time) {};
}

void Model::step(const unsigned steps) {
  for (unsigned i(0); i != steps; ++i) {
      stepper_.step();
    }
}

unsigned Model::push_species(Species& species) {
  species_list_.push_back(&species);
  if (species_list_.size() > MAX_SPECIES) {
      std::stringstream error_message;
      error_message << "\nMAX_SPECIES is set to " << MAX_SPECIES <<
        ", which is less than the actual size " << species_list_.size() << 
        ".\nSet MAX_SPECIES to " << species_list_.size() << ".";
      throw std::runtime_error(error_message.str());
  }
  return species_list_.size()-1;
}

void Model::push_process(Process& process) {
  process_list_.push_back(&process);
}

SpaceCompartment& Model::get_space_compartment() {
  return space_compartment_;
}

VisualWriter& Model::get_visual_writer() {
  return visual_writer_;
}

NumberWriter& Model::get_number_writer() {
  return number_writer_;
}

TrackWriter& Model::get_track_writer() {
  return track_writer_;
}

Stepper& Model::get_stepper() {
  return stepper_;
}

std::vector<Species*>& Model::get_species_list() {
  return species_list_;
}

std::mt19937_64& Model::get_rng() {
  return rng_;
}

gsl_rng* Model::get_gsl_rng() {
  return gsl_rng_;
}

unsigned Model::get_next_agent_cid() {
  return ++agent_cid_-1;
}
