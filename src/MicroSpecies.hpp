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


#ifndef __MicroSpecies_hpp
#define __MicroSpecies_hpp

#include <Common.hpp>
#include <MotileSpecies.hpp>

class MicroSpecies: public MotileSpecies { 
public: 
  MicroSpecies(const std::string, Model&, VolumeSpecies&,
               const unsigned init_size = 1, const double D = 0);
  MicroSpecies(const std::string, Model&, const unsigned init_size = 1,
               const double D = 0);
  ~MicroSpecies() {}
  std::vector<Vector<double>>& get_init_positions(); 
  virtual std::vector<Vector<double>>& get_xyz_list(); 
  virtual std::vector<Vector<double>>& get_relative_positions(); 
  virtual void walk();
  virtual unsigned size();
  void set_size_changed(bool);
  void set_relative_positions_changed(bool);
  bool size_changed();
  bool relative_positions_changed();
  void remove_xyz(const unsigned index);
  void add_xyz(const Vector<double>& xyz);
private:
  bool is_size_changed = true;
  bool is_relative_positions_changed = true;
  std::vector<Vector<double>> xyz_list_;
};

#endif /* __MicroSpecies_hpp */

