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


#ifndef __TrackReader_hpp
#define __TrackReader_hpp

#include <fstream>
#include <map>
#include <Common.hpp>

class TrackReader {
public:
  TrackReader(std::string filename);
  ~TrackReader() {}
  void parse_file();
  Vector<double>& get_dimensions();
  double get_interval() const;
  double get_end_time() const;
  unsigned get_size();
  unsigned get_mean_agents_per_step() const;
  unsigned get_initial_track_id_size(const int id);
  unsigned get_initial_track_ids_size();
  std::vector<unsigned>& get_track_ids_at_step(const unsigned step);
  std::vector<Vector<double>>& get_xyz_list_at_step(const unsigned step);

private:
  std::string filename_;
  std::ifstream track_file_;
  double interval_;
  Vector<double> dimensions_;
  Vector<double> max_ = Vector<double>(-1, -1, -1);
  Vector<double> min_ = Vector<double>(std::numeric_limits<double>::infinity(), 
                                     std::numeric_limits<double>::infinity(), 
                                     std::numeric_limits<double>::infinity());
  unsigned agent_id_cnt_ = 0;
  std::map<unsigned, unsigned> agent_id_map_;
  std::vector<std::vector<unsigned>> step_track_ids_;
  std::vector<std::vector<Vector<double>>> step_xyz_list_;
};

#endif /* __TrackReader_hpp */
