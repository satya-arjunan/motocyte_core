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


#ifndef __Populator_hpp
#define __Populator_hpp

#include <Common.hpp>
#include <Process.hpp>

class Populator: public Process { 
public: 
  Populator(Model& model, MotileSpecies&,
              Vector<double> origin = Vector<double>(0, 0, 0),
              Vector<double> cuboid = Vector<double>(1, 1, 1));
  ~Populator() {}
  void set_origin(const double, const double, const double);
  void set_cuboid(const double, const double, const double);
  void set_radius(const double radius);
  void populate();
  bool get_is_populated() const;
  void set_is_populated();
  void add_agent(const Vector<double>& xyz);
  void remove_agent(const unsigned agent_idx, const unsigned back_agent_idx);
protected:
  void populate_volume_species();
  void populate_track_species();
  void populate_meso_species();
private:
  MotileSpecies& motile_species_;
  VolumeSpecies* volume_species_;
  TrackSpecies* track_species_;
  MesoSpecies* meso_species_;
  bool is_populated_;
  double radius_ = -1;
  Vector<double> origin_;
  Vector<double> cuboid_;
};

#endif /* __Populator_hpp */

