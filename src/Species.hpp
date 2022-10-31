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


#ifndef __Species_hpp
#define __Species_hpp

#include <Common.hpp>

class Species { 
public: 
  Species(const std::string name, Model& model,
          VolumeSpecies& compartment, const unsigned init_size = 1);
  Species(const std::string name, Model& model, const unsigned init_size = 1);
  ~Species() {}
  virtual void initialize () {}
  VolumeSpecies& get_compartment();
  Model& get_model();
  const std::string& get_name() const;
  std::string get_name_id() const;
  unsigned get_species_id() const;
  unsigned get_init_size() const;
  virtual unsigned& get_size();
  virtual unsigned size();
  void set_size(unsigned);
  void set_species_cid(const unsigned);
  unsigned get_species_cid() const;
private:
  std::string get_init_name(const std::string) const;
private:
  const unsigned init_size_;
  unsigned size_;
  VolumeSpecies& compartment_;
  const std::string name_;
  Model& model_;
  unsigned species_id_; //species_id unique to the whole model
  unsigned species_cid_; //species_id unique to the compartment
};

#endif /* __Species_hpp */

