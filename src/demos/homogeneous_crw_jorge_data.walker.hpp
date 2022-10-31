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


#ifndef __Walker_hpp
#define __Walker_hpp

#include <random>
#include <Common.hpp>
#include <Process.hpp>
#include <Quaternion.hpp>

typedef const Vector<float> (Walker::*DisplaceMethod)(const unsigned);

class Walker: public Process { 
public: 
  Walker(const double, Model& model, MotileSpecies&);
  ~Walker() {}
  const Vector<float> get_displacement(const unsigned);
  void set_directed_displacement();
  void set_homogeneous_crw();
  float step();
  const Vector<float> get_random_displacement(const unsigned);
  const Vector<float> get_directed_displacement(const unsigned);
  const Vector<float> get_homogeneous_crw_displacement(const unsigned agent_id);
  std::vector<Vector<float>>& get_directions();
  std::vector<unsigned>& get_reorientation_magnitudes();
  virtual void set_interval(float);
protected:
  void init_orientations();
  const Quaternion get_homogeneous_crw_orientation(const unsigned agent_id);
  std::vector<Vector<double>> get_random_directions(const unsigned size);
private:
  const float speed_mean_ = 0.115;
  const float speed_std_ = 0.3;
  const double max_speed_ = 0.68;
  const float roll_rate_mean_ = -1.0; 
  const float roll_rate_std_ = 0.0; 
  const float pitch_rate_mean_ = 0.0;
  const float pitch_rate_std_ = 0.1;
  double D_;
  MotileSpecies& motile_species_;
  std::normal_distribution<> norm_dist_;
  std::normal_distribution<> speed_norm_dist_;
  std::normal_distribution<> roll_norm_dist_;
  std::normal_distribution<> pitch_norm_dist_;
  std::uniform_real_distribution<> uni_dist_;
  DisplaceMethod displace_method_;
  std::vector<Vector<float>> directions_;
  std::vector<unsigned> reorientation_magnitudes_;
  std::vector<Quaternion> orientations_;
};

#endif /* __Walker_hpp */

