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


#ifndef __MacroSpecies_hpp
#define __MacroSpecies_hpp

#include <Common.hpp>
#include <Species.hpp>

class MacroSpecies: public Species { 
public: 
  MacroSpecies(const std::string, Model&, VolumeSpecies&,
               const unsigned init_size = 1);
  MacroSpecies(const std::string, Model&, const unsigned init_size = 1);
  ~MacroSpecies() {}
  unsigned get_size() const;
 private:
  unsigned size_;
};

#endif /* __MacroSpecies_hpp */

