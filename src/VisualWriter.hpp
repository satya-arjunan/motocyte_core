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


#ifndef __VisualWriter_hpp
#define __VisualWriter_hpp

#include <fstream>
#include <Common.hpp>
#include <MicroSpecies.hpp>
#include <MesoSpecies.hpp>
#include <VolumeSpecies.hpp>
#include <PointSpecies.hpp>
#include <XDMFFile.hpp>
#include <Process.hpp>

class VisualWriter: public Process {
public:
  VisualWriter(Model& model, const bool write_xdmf = false);
  ~VisualWriter() {}
  double step();
  virtual void initialize();
  void add(MicroSpecies&);
  void add(MesoSpecies&);
private:
  void initialize_log();
  void log_structure_species();
  void log_species();
  void log_xdmf(MicroSpecies&);
  void log_xdmf(MesoSpecies&);
  void log_points(MesoSpecies&, const unsigned);
  void log_points(VolumeSpecies&, const unsigned);
  void log_points(PointSpecies&, const unsigned);
private:
  const bool write_xdmf_;
  unsigned marker_;
  std::string filename_;
  std::ofstream logfile_;
  SpaceCompartment& space_compartment_;
  std::vector<Species*> species_;
  XDMFFile micro_xdmf_file_;
  XDMFFile meso_xdmf_file_;
};

#endif /* __VisualWriter_hpp */
