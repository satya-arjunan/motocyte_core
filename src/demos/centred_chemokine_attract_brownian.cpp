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

  //Chemokine source
  VolumeSpecies A("A", model, 100, 0);
  A.get_populator().set_cuboid(0.14, 0.14, 0.14);
  PointSpecies Ap("Ap", A, 10);
  MesoSpecies Am1("Am1", model, 0, 1);
  Am1.get_walker().set_interval(0.1);
  Reaction exocytosis(Ap, 500, Am1); //Ap -> Am1
  exocytosis.set_interval(0.2);
  //Adjust degradation rate to control radius of chemokine spread:
  Reaction degradation(Am1, 0.02); //Am1 -> 0
  degradation.set_interval(0.2);

  //T cells
  VolumeSpecies B("B", model, 100, 15);
  B.get_walker().set_interval(0.2);
  //B.get_walker().set_homogeneous_crw();
  PointSpecies Bp("Bp", B, 1000);
  Reaction chemotaxis(Bp, Am1, 1); //Bp + Am1 -> change direction
  chemotaxis.set_interval(0.2);

  model.get_writer().set_interval(1);
  model.get_writer().add(Bp);
  model.get_writer().add(Ap);
  model.get_writer().add(Am1);

  model.initialize();

  model.step(10000);
  /*

  exocytosis1.set_rate(0);
  exocytosis2.set_rate(5);

  model.step(1000);
  */

  return 0;
}

