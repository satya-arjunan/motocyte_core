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


#ifndef __Process_hpp
#define __Process_hpp

#include <random>
#include <limits>
#include <Common.hpp>
#include <PriorityQueue.hpp>

typedef PriorityQueue<Process*> ProcessQueue;
typedef ProcessQueue::ID QueueID;

class Process
{ 
public: 
  Process(const std::string, Model&,
          const double interval = std::numeric_limits<double>::infinity());
  ~Process() {}
  void set_queue_id(QueueID);
  void set_priority(int);
  void set_time(double);
  virtual void set_interval(double);
  QueueID get_queue_id() const;
  int get_priority() const;
  double get_time() const;
  double get_interval() const;
  virtual double step();
  virtual void initialize() {};
  virtual const std::string get_name() const;
  Stepper& get_stepper();
  Model& get_model();
private:
  const std::string name_;
  Model& model_;
  double interval_;
  QueueID queue_id_;
  int priority_;
  double time_;
};

#endif /* __Process_hpp */

