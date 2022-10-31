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


#ifndef __TrackWriter_hpp
#define __TrackWriter_hpp

#include <fstream>
#include <Common.hpp>
#include <MicroSpecies.hpp>
#include <MesoSpecies.hpp>
#include <VolumeSpecies.hpp>
#include <PointSpecies.hpp>
#include <XDMFFile.hpp>
#include <NumberWriter.hpp>

class TrackWriter: public NumberWriter {
public:
  TrackWriter(Model& model);
  ~TrackWriter() {}
  virtual void add(MesoSpecies&);
  virtual void add(MicroSpecies&);
  virtual void initialize();
  void set_padding(const double padding);
  void set_skip_frames(const unsigned skip_frames);
protected:
  virtual void initialize_log();
  virtual void log_species();
private:
  unsigned skip_frames_ = 0;
  unsigned frame_cnt_ = 0;
  std::vector<Vector<double>> prev_agents_;
  VolumeSpecies* volume_species_;
  double padding_ = 0;
  Vector<double> rect_min_;
  Vector<double> rect_max_;
  std::vector<unsigned> agent_cids_;
  std::vector<unsigned> active_cids_;
  unsigned id_cnt_ = 0;
};

#endif /* __TrackWriter_hpp */
