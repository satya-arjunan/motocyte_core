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


#ifndef __PointSpecies_hpp
#define __PointSpecies_hpp

#include <MicroSpecies.hpp>

class PointSpecies: public MicroSpecies { 
 public: 
  PointSpecies(const std::string name, VolumeSpecies& compartment,
        const unsigned init_size = 0, const double D = 0);
  ~PointSpecies() {}
  virtual void initialize() {}
  void populate();
  virtual std::vector<Vector<double>>& get_positions(); 
  virtual std::vector<Vector<double>>& get_relative_positions();
  const std::vector<Vector<double>>& get_direction_vectors() const;
  const std::vector<unsigned>& get_compartment_ids();
  const std::vector<unsigned>& get_compartment_points_sizes() const;
  unsigned size();
  void remove_points(const unsigned index);
  void add_points();
 private:
  std::uniform_real_distribution<> uni_z_;
  std::uniform_real_distribution<> uni_t_;
  std::vector<unsigned> compartment_points_sizes_;
  std::vector<Vector<double>> relative_positions_;
  std::vector<Vector<double>> direction_vectors_;
};



#endif /* __PointSpecies_hpp */

