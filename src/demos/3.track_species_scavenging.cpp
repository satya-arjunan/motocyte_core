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
  std::string filename("20.628_Position.csv");
  TrackReader tracks(filename);
  const Vector<float> dims(tracks.get_dimensions());
  const float agent_radius(5);
  Model model(dims, agent_radius);
  TrackSpecies A("A", model, tracks);
  PointSpecies Ap("Ap", A, 200);

  MesoSpecies Ac("Ac", model, 1, 0);
  Reaction scavenging(Ap, Ac, 50, Ap); //Ap + Ac -> Ap
  scavenging.set_interval(5);

  model.get_visual_writer().set_interval(tracks.get_interval());
  model.get_visual_writer().add(A);
  model.get_visual_writer().add(Ap);
  model.get_visual_writer().add(Ac);

  model.initialize();
  //model.step(50);
  model.run(tracks.get_end_time());
  return 0;
}


