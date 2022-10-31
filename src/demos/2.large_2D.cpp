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
#include <PointSpecies.hpp>
#include <Reaction.hpp>
#include <Populator.hpp>

int main() {
  Model model;
  VolumeSpecies A("A", model, 300);
  A.get_walker().set_interval(5);
  A.get_walker().set_homogeneous_crw();
  PointSpecies Ap("Ap", A, 10);

  MesoSpecies Ac("Ac", model, 1, 0);

  MesoSpecies Am("Am", model, 0);
  Am.get_walker().set_interval(5);


  Reaction exocytosis(Ap, 1, Am); //Ap -> Am
  exocytosis.set_interval(5);
  Reaction degradation(Am, 0.005); //Am -> 0
  degradation.set_interval(4);


  PointSpecies As("As", A, 20);
  Reaction chemotaxis(As, Am, 1000); //As + Am -> change direction
  chemotaxis.set_interval(4);

  /*

  //model.get_visual_writer().add(Ac);

  model.get_number_writer().set_interval(100);
  model.get_number_writer().add(Ac);

  model.get_track_writer().set_interval(20);
  model.get_track_writer().set_padding(40);
  model.get_track_writer().add(A);
  */

  model.get_visual_writer().set_interval(100);
  model.get_visual_writer().add(A);
  model.get_visual_writer().add(Ap);
  model.get_visual_writer().add(Am);

  model.initialize();
  //model.step(50);
  model.run(100000);
  return 0;
}

