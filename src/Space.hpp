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


#ifndef __Space_hpp
#define __Space_hpp

#include <Common.hpp>
#include <Motocyte.hpp>

class Space {
public: 
  Space(SpaceCompartment& space_compartment, const Vector<double>& dimensions,
        const double agent_radius, const int Mult, const int Div);
  ~Space() {}
  virtual void populate_all() = 0;
  virtual void walk() = 0;
  virtual const Vector<double>& get_dimensions() const;
  virtual const Vector<double>& get_dimensions_mod() const;
  virtual void push_walker(Walker& walker);
  const size_t xyz_to_voxel_id(const Vector<double>& xyz) const;
  const double get_voxel_length() const;
  const Vector<double> voxel_id_to_xyz(const size_t voxel_id) const;
  const Vector<unsigned> xyz_to_coord(const Vector<double>& xyz) const;
  const size_t coord_to_voxel_id(const Vector<unsigned>& coord) const;
  const Vector<unsigned> voxel_id_to_coord(const size_t voxel_id) const;
  Model& get_model();
  SpaceCompartment& get_space_compartment();
  Species& get_species(const unsigned species_id);
protected:
  const double voxel_length_;
  const uint16_t dim_x_;
  const uint16_t dim_y_;
  const uint16_t dim_z_;
  const unsigned num_voxels_;
  const Vector<double> dimensions_;
  const Vector<double> dimensions_mod_;
  SpaceCompartment& space_compartment_;
  std::vector<Walker*> walkers_;
  Model& model_;
};

#endif /* __Space_hpp */

