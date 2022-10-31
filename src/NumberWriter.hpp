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


#ifndef __NumberWriter_hpp
#define __NumberWriter_hpp

#include <fstream>
#include <Common.hpp>
#include <MicroSpecies.hpp>
#include <MesoSpecies.hpp>
#include <VolumeSpecies.hpp>
#include <PointSpecies.hpp>
#include <XDMFFile.hpp>
#include <Process.hpp>

class NumberWriter: public Process {
public:
  NumberWriter(Model& model);
  ~NumberWriter() {}
  virtual double step();
  virtual void initialize();
  virtual void add(MicroSpecies&);
  virtual void add(MesoSpecies&);
  std::ofstream& get_log_file();
  void set_file_name(std::string file_name);
protected:
  virtual void initialize_log();
  virtual void log_species();
  virtual void log(VolumeSpecies&, const unsigned);
  void log(MesoSpecies&, const unsigned);
  void log(PointSpecies&, const unsigned);
private:
  std::string file_name_;
  std::ofstream log_file_;
  SpaceCompartment& space_compartment_;
  std::vector<Species*> species_;
};

#endif /* __NumberWriter_hpp */
