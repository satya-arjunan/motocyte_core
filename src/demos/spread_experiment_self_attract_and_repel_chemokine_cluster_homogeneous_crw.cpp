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

  //T cells
  VolumeSpecies B("B", model, 50, 15);
  //B.get_populator().set_cuboid(0.14, 0.14, 1);
  B.get_walker().set_interval(0.2);
  B.get_walker().set_homogeneous_crw();

  MesoSpecies Bm("Bm", model, 0, 1);
  Bm.get_walker().set_interval(0.5);

  PointSpecies Bp("Bp", B, 100);
  Reaction exocytosis(Bp, 20, Bm); //Bp -> Bm
  exocytosis.set_interval(0.2);
  Reaction degradation(Bm, 0.1); //Bm -> 0
  degradation.set_interval(0.2);

  PointSpecies Bs("Bs", B, 500);
  //higher k, lower chemokine potency
  Reaction chemotaxis(Bs, Bm, 20); //Bs + Bm -> change direction
  chemotaxis.set_interval(0.2);

  MesoSpecies Bm2("Bm2", model, 0, 1);
  Bm2.get_walker().set_interval(0.5);

  PointSpecies Bd("Bd", B, 100);
  Reaction exocytosis2(Bd, 50, Bm2); //Bp -> Bm
  exocytosis2.set_interval(0.2);
  Reaction degradation2(Bm2, 0.4); //Bm -> 0
  degradation2.set_interval(0.2);

  PointSpecies Bt("Bt", B, 500);
  Reaction chemotaxis2(Bt, Bm2, -1); //Bs + Bm -> change direction
  chemotaxis2.set_interval(0.2);

  model.get_writer().set_interval(0.2);
  model.get_writer().add(Bm);
  model.get_writer().add(Bm2);
  model.get_writer().add(Bp);
  model.get_writer().add(Bs);
  model.get_writer().add(Bd);
  model.get_writer().add(Bt);

  model.initialize();

  model.step(10000);

  return 0;
}

