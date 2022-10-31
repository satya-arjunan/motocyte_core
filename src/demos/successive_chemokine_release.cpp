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

  VolumeSpecies A("A", model, 100, 0);
  A.get_populator().set_cuboid(0.8, 0.1, 0.1);

  VolumeSpecies B("B", model, 1000, 1e-12);
  B.get_walker().set_interval(0.1e-6);
  VolumeSpecies C("C", model, 1000, 1e-12);
  C.get_walker().set_interval(0.1e-6);
  PointSpecies Ap("Ap", A, 10);
  PointSpecies Bp("Bp", B, 100);
  PointSpecies Cp("Cp", C, 100);
  MesoSpecies Am1("Am1", model, 0, 1e-12);
  Am1.get_walker().set_interval(0.1e-6);
  MesoSpecies Am2("Am2", model, 0, 1e-12);
  Am2.get_walker().set_interval(0.1e-6);

  Reaction exocytosis1(Ap, 5, Am1); //Ap -> Am1
  exocytosis1.set_interval(0.1e-6);
  Reaction exocytosis2(Ap, 0, Am2); //Ap -> Am2
  exocytosis2.set_interval(0.1e-6);
  Reaction chemotaxis1(Bp, Am1, 5); //Bp + Am1 -> change direction
  chemotaxis1.set_interval(0.1e-6);
  Reaction chemotaxis2(Cp, Am2, 5); //Cp + Am2 -> change direction
  chemotaxis2.set_interval(0.1e-6);
  Reaction degradation1(Am1, 30000); //Am1 -> 0
  degradation1.set_interval(0.1e-6);
  Reaction degradation2(Am2, 30000); //Am2 -> 0
  degradation2.set_interval(0.1e-6);

  model.get_writer().set_interval(1e-6);
  model.get_writer().add(A);
  model.get_writer().add(B);
  model.get_writer().add(C);
  model.get_writer().add(Ap);
  model.get_writer().add(Bp);
  model.get_writer().add(Cp);
  model.get_writer().add(Am1);
  model.get_writer().add(Am2);

  model.initialize();

  model.step(6000);

  exocytosis1.set_rate(0);
  exocytosis2.set_rate(5);

  model.step(6000);

  return 0;
}

