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

#include <Process.hpp>
#include <Model.hpp>

Process::Process(const std::string name, Model& model, const double interval):
  name_(name),
  model_(model),
  interval_(interval),
  priority_(0),
  time_(interval) {
    model_.push_process(*this);
  }

double Process::step() {
  time_ += interval_;
  get_stepper().get_process_queue().move_top();
  return time_;
}

Stepper& Process::get_stepper() {
  return get_model().get_stepper();
}

Model& Process::get_model() {
  return model_;
}

void Process::set_priority(const int priority) {
  priority_ = priority;
}

void Process::set_time(const double time) {
  time_ = time;
}

const std::string Process::get_name() const {
  return name_;
}

void Process::set_interval(const double interval) {
  const double old_time(time_);
  if (time_ == std::numeric_limits<double>::infinity()) {
    time_ = get_stepper().get_time()+interval;
  }
  else {
    time_ = time_-interval_+interval;
  }
  interval_ = interval;
  if (time_ >= old_time) {
    get_stepper().get_process_queue().move_down(queue_id_);
  } else {
    get_stepper().get_process_queue().move_up(queue_id_);
  }          
}

void Process::set_queue_id(QueueID queue_id) {
  queue_id_ = queue_id;
}

int Process::get_priority() const {
  return priority_;
}

double Process::get_time() const {
  return time_;
}

double Process::get_interval() const {
  return interval_;
}

QueueID Process::get_queue_id() const {
  return queue_id_;
}


