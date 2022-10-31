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


#ifndef __MesoSpecies_hpp
#define __MesoSpecies_hpp

#include <Common.hpp>
#include <MotileSpecies.hpp>

class MesoSpecies: public MotileSpecies { 
public: 
  MesoSpecies(const std::string, Model&, SpaceCompartment&,
              const unsigned init_size = 1, const double D = 0);
  MesoSpecies(const std::string, Model&, const unsigned init_size = 1,
              const double D = 0);
  ~MesoSpecies() {}
  std::vector<unsigned>& get_sizes();
  MesoSpace& get_meso_space();
  const std::vector<std::vector<double>>& get_rect_mesh();
  SpaceCompartment& get_space_compartment();
  std::vector<unsigned>& get_occupied_voxels();
  void walk();
 private:
  std::vector<unsigned> sizes_;
  SpaceCompartment& space_compartment_;
};

#endif /* __MesoSpecies_hpp */

