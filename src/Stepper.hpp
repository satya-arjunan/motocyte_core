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


#ifndef __Stepper_hpp
#define __Stepper_hpp

#include <Common.hpp>
#include <SpaceCompartment.hpp>
#include <PriorityQueue.hpp>
#include <Process.hpp>

class Stepper { 
public: 
  Stepper();
   ~Stepper() {}
  double step();
  double get_time() const;
  void set_time(const double);
  ProcessQueue& get_process_queue();
  void initialize();
private:
  double time_;
  ProcessQueue process_queue_; 
};

#endif /* __Stepper_hpp */

