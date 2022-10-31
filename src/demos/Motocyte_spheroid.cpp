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

  VolumeSpecies AA("AA", model, 0, 0);
  VolumeSpecies A("A", model, 550, 0);
  A.get_populator().set_cuboid(0.9, 1.0, 1.0);
  A.get_populator().set_origin(0.0, 0.0, 0.0);
  A.get_populator().set_radius(0.12);
  //VolumeSpecies A("A", model, 1550, 0);
  //A.get_populator().set_cuboid(0.1, 1.0, 1.0);
  //A.get_populator().set_origin(0.25, 0, 0);

  VolumeSpecies B("B", model, 100, 1e-12);
  VolumeSpecies C("C", model, 100, 1e-12);
  //PointSpecies Dp("Dp", D, 500);
  PointSpecies AAp("Ap", AA, 500);
  PointSpecies Ap("Ap", A, 500);
  PointSpecies App("App", A, 1000);
  PointSpecies Bpp("Bpp", B, 1000);
  PointSpecies Cpp("Cpp", C, 1000);
  PointSpecies Bp("Bp", B, 500);
  PointSpecies Bpr("Bpr", B, 6000);
  PointSpecies Cp("Cp", C, 500);
  MesoSpecies Am1("Am1", model, 0, 1e-12);
  MesoSpecies Bm("Bm", model, 0, 1e-12);

  Reaction exocytosis1(Ap, 5, Am1); //Ap -> Am1
  exocytosis1.set_interval(1e-5);
  Reaction chemotaxis1(Bp, Am1, 1); //Bp + Am1 -> change direction
  Reaction syntheticligand(Bpr, Am1, 15, Bm); //Bp + Am1 -> release Bm
  Reaction chemotaxis2(Cp, Bm, 7); //Cp + Bm -> change direction

  Reaction degradation1(Am1, 50000); //Am1 -> 0
  degradation1.set_interval(1e-5);
  Reaction degradation2(Bm, 5000); //Bm-> 0
  degradation2.set_interval(1e-5);

  model.get_writer().set_interval(1e-6);
  //model.get_writer().add(Dp);
  model.get_writer().add(AAp);
  //model.get_writer().add(Am1);
  model.get_writer().add(Bm);
  model.get_writer().add(App);
  model.get_writer().add(Bpp);
  model.get_writer().add(Cpp);

  model.initialize();

  model.step(4000);
  /*

  exocytosis1.set_rate(0);
  exocytosis2.set_rate(5);

  model.step(1000);
  */

  return 0;
}

