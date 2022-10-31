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


#ifndef __SpaceCompartment_hpp
#define __SpaceCompartment_hpp

#include <VolumeSpecies.hpp>
#include <MesoSpecies.hpp>
#include <MicroSpace.hpp>
#include <MesoSpace.hpp>

class SpaceCompartment: public VolumeSpecies { 
public: 
  SpaceCompartment(const std::string name, Model& model,
                   const Vector<double>& dimensions, const double agent_radius,
                   SpaceCompartment& space_compartment,
                   const unsigned init_size = 1, const double D = 0);
  SpaceCompartment(const std::string name, Model& model, 
                   const Vector<double>& dimensions, const double agent_radius,
                   const unsigned init_size = 1, const double D = 0);

  ~SpaceCompartment() {}
  void populate_all();
  void update_species();
  const Vector<double>& get_dimensions() const;
  Vector<double> get_center() const;
  MicroSpace& get_micro_space();
  MesoSpace& get_meso_space();
  std::vector<std::reference_wrapper<MicroSpecies>>& get_micro_species_list();
  std::vector<std::reference_wrapper<MesoSpecies>>& get_meso_species_list();
  void push_species(Species&);
  void push_species(MesoSpecies&);
  void push_species(MicroSpecies&);
 private:
  std::vector<std::reference_wrapper<MicroSpecies>> micro_species_list_;
  std::vector<std::reference_wrapper<MesoSpecies>> meso_species_list_;
  MicroSpace micro_space_;
  MesoSpace meso_space_;
};

#endif /* __SpaceCompartment_hpp */

