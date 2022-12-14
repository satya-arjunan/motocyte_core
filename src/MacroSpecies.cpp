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

#include <Model.hpp>
#include <MacroSpecies.hpp>

MacroSpecies::MacroSpecies(const std::string name, Model& model,
                           VolumeSpecies& compartment,
                           const unsigned init_size):
  Species(name, model, compartment, init_size),
  size_(init_size) {}

MacroSpecies::MacroSpecies(const std::string name, Model& model,
                           const unsigned init_size):
  MacroSpecies(name, model, model.get_space_compartment(), init_size) {}

unsigned MacroSpecies::get_size() const {
  return size_;
}
