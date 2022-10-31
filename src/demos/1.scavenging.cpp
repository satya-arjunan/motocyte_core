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
#include <VisualLogger.hpp>
#include <PointSpecies.hpp>
#include <Reaction.hpp>
#include <Populator.hpp>

int main() {
  Model model;
  VolumeSpecies A("A", model, 44);
  A.get_walker().set_interval(5);
  A.get_walker().set_homogeneous_crw();
  PointSpecies Ap("Ap", A, 200);

  /*
  MesoSpecies Am("Am", model, 0, 0);
  Am.get_walker().set_interval(2);
  */

  MesoSpecies Ac("Ac", model, 1, 0);

  Reaction scavenging(Ap, Ac, 50, Ap); //Ap + Ac -> Ap
  scavenging.set_interval(5);

  //Reaction degradation2(Ac, 0.005); //Ac -> 0
  //degradation2.set_interval(5);


  /*
  Reaction exocytosis(Ap, 50, Am); //Ap -> Am
  exocytosis.set_interval(4);
  Reaction degradation(Am, 0.02); //Am -> 0
  degradation.set_interval(4);
  */

  /*
  PointSpecies As("As", A, 500);
  Reaction chemotaxis(As, Am, -1); //As + Am -> change direction
  chemotaxis.set_interval(4);

  model.get_writer().add(As);
  model.get_writer().add(Am);
  */

  model.get_writer().set_interval(100);
  model.get_writer().add(A);
  model.get_writer().add(Ap);
  model.get_writer().add(Ac);
  model.initialize();
  //model.step(50);
  model.run(100000);
  return 0;
}

