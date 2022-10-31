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

#include <iostream> 
#include <boost/date_time/posix_time/posix_time.hpp>
#include <Motocyte.hpp>
#include <Model.hpp>
#include <VisualWriter.hpp>
#include <NumberWriter.hpp>
#include <TrackReader.hpp>
#include <TrackSpecies.hpp>
#include <PointSpecies.hpp>
#include <Reaction.hpp>
#include <Populator.hpp>

int main() {
  const Vector<float> dims(383, 383, 60);
  const float agent_radius(5);
  Model model(dims, agent_radius);
  MesoSpecies Am("Am", model, 100000, 1);
  Am.get_populator().set_cuboid(0.01, 0.01, 0.01);
  model.get_visual_writer().set_interval(20);
  model.get_visual_writer().add(Am);

  model.initialize();
  //model.step(50);
  model.run(700);
  return 0;
}


