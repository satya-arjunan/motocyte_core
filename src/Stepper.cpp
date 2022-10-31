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

#include <Stepper.hpp>
#include <Walker.hpp>
#include <Reaction.hpp>
#include <VisualWriter.hpp>
#include <MicroSpace.hpp>

Stepper::Stepper():
  time_(0) {
    get_process_queue().clear();
}

double Stepper::step() {
  //std::cout << "\n stepping" << std::endl;
  do {
    /*
    std::cout << "curr time:" << time_ << std::endl;
    std::cout << "executing:" << get_process_queue().get_top()->get_name() <<
      " interval:" << get_process_queue().get_top()->get_interval() <<
      std::endl;
      */
    get_process_queue().get_top()->step();

  } while (get_process_queue().get_top()->get_time() == get_time());
  set_time(get_process_queue().get_top()->get_time());
  return get_time();
}

void Stepper::initialize() {
  set_time(get_process_queue().get_top()->get_time());
}

ProcessQueue& Stepper::get_process_queue() {
  return process_queue_;
}

double Stepper::get_time() const {
  return time_;
}

void Stepper::set_time(const double time) {
  time_ = time;
}

