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


#ifndef __VolumeSpecies_hpp
#define __VolumeSpecies_hpp

#include <MicroSpecies.hpp>
#include <MesoSpecies.hpp>

class VolumeSpecies: public MicroSpecies { 
 public: 
  VolumeSpecies(const std::string name, Model& model, SpaceCompartment&,
              const unsigned init_size = 1, const double D = 0);
  VolumeSpecies(const std::string name, Model& model, 
              const unsigned init_size = 1, const double D = 0);
  ~VolumeSpecies() {}
  virtual void initialize() {}
  virtual void populate_all();
  virtual void update_species() {}
  virtual void push_species(Species&);
  virtual void push_species(PointSpecies&);
  virtual void push_species(VolumeSpecies&);
  double get_agent_radius() const;
  void walk();
  Vector<double> get_center() const;
  SpaceCompartment& get_space_compartment();
  std::vector<std::reference_wrapper<Species>>& get_species_list();
  std::vector<unsigned>& get_agent_cids();
  void remove_agent(const unsigned agent_idx);
  void add_agent(const Vector<double>& xyz);
  void set_check_collision(const bool check_collision);
 private:
  void populate_compartment();
  void populate_child_compartments();
 private:
  std::vector<unsigned> agent_cids_; //agent ids unique across whole space
  SpaceCompartment& space_compartment_;
  std::vector<std::reference_wrapper<Species>> species_list_;
  std::vector<std::reference_wrapper<VolumeSpecies>> compartment_species_list_;
  std::vector<PointSpecies*> point_species_list_;
};

#endif /* __VolumeSpecies_hpp */

