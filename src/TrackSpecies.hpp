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


#ifndef __TrackSpecies_hpp
#define __TrackSpecies_hpp

#include <VolumeSpecies.hpp>
#include <TrackReader.hpp>

class TrackSpecies: public VolumeSpecies { 
 public: 
  TrackSpecies(const std::string name, Model& model, TrackReader& tracks,
               const int size);
  TrackSpecies(const std::string name, Model& model, TrackReader& tracks);
  ~TrackSpecies() {}
  virtual std::vector<Vector<double>> get_initial_xyz_list();
  virtual void walk();
  TrackReader& get_track_reader();
  std::vector<unsigned>& get_track_ids();
 protected:
  TrackReader& track_reader_;
  std::vector<unsigned> track_ids_;
};

#endif /* __TrackSpecies_hpp */

