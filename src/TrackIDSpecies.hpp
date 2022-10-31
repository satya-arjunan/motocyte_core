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


#ifndef __TrackIDSpecies_hpp
#define __TrackIDSpecies_hpp

#include <VolumeSpecies.hpp>
#include <TrackSpecies.hpp>

class TrackIDSpecies: public TrackSpecies { 
 public: 
  TrackIDSpecies(const std::string name, Model& model, TrackReader& tracks,
                 const unsigned track_id); 
  ~TrackIDSpecies() {}
  virtual std::vector<Vector<double>> get_initial_xyz_list();
  virtual void walk();
 private:
  const unsigned track_id_;
};

#endif /* __TrackIDSpecies_hpp */

