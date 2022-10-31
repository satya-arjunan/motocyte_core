//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Motocyte
//
//        Copyright (C) 2019-2020 Satya N.V. Arjunan
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

#include <XMLModelWriter.hpp>

XMLModelWriter::XMLModelWriter(std::string file) {
  auto declarationNode = doc_.append_child(pugi::node_declaration);
  declarationNode.append_attribute("version")    = "1.0";
  declarationNode.append_attribute("encoding")   = "ISO-8859-1";
  declarationNode.append_attribute("standalone") = "yes";
  //auto root = doc_.append_child("MyRoot");
  doc_.append_child("MyRoot");
  //bool saveSucceeded = doc_.save_file(file.c_str(), PUGIXML_TEXT("  "));
  doc_.save_file(file.c_str(), PUGIXML_TEXT("  "));
}

void XMLModelWriter::set_dimensions(const Vector<double>& dimensions) {
  dimensions_ = dimensions;
}

void XMLModelWriter::set_size(unsigned size) {
  size_ = size;
}

void XMLModelWriter::set_interval(double interval) {
  interval_ = interval;
} 

void XMLModelWriter::set_end_time(double end_time) {
  steps_ = end_time/interval_;
} 
