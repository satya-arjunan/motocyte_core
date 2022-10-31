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

#include <sstream>
#include <climits>
#include <Motocyte.hpp>
#include <VisualWriter.hpp>
#include <Model.hpp>
#include <Stepper.hpp>


VisualWriter::VisualWriter(Model& model, const bool write_xdmf):
  Process("VisualWriter", model, std::numeric_limits<double>::infinity()),
  write_xdmf_(write_xdmf),
  marker_(std::numeric_limits<unsigned>::max()),
  filename_("VisualLog.dat"),
  space_compartment_(model.get_space_compartment()),
  micro_xdmf_file_(XDMFFile("motocyte0.xdmf")),
  meso_xdmf_file_(XDMFFile("motocyte1.xdmf")) {
    set_queue_id(get_stepper().get_process_queue().push(this));
    set_priority(-5); //lowest priority in process queue
  }

double VisualWriter::step() {
  log_species();
  logfile_.flush();
  return Process::step();
}

void VisualWriter::initialize() {
  std::ostringstream fileName;
  fileName << filename_ << std::ends;
  logfile_.open(fileName.str().c_str(), std::ios::binary | std::ios::trunc);
  initialize_log();
  log_structure_species();
  log_species();
  logfile_.flush();
}

void VisualWriter::add(MicroSpecies& micro_species) {
  species_.push_back(reinterpret_cast<Species*>(&micro_species));
}

void VisualWriter::add(MesoSpecies& meso_species) {
  species_.push_back(reinterpret_cast<Species*>(&meso_species));
}

void VisualWriter::initialize_log() {
  const unsigned latticeType(1); //CUBIC
  logfile_.write((char*)(&latticeType), sizeof(latticeType));
  const unsigned meanCount(0);
  logfile_.write((char*)(&meanCount), sizeof(meanCount));
  const unsigned startCoord(0);
  logfile_.write((char*)(&startCoord), sizeof(startCoord));
  const Vector<double>& real_dimensions(space_compartment_.get_dimensions());
  const double agent_radius(space_compartment_.get_agent_radius());
  const Vector<unsigned> dimensions(
                              unsigned(real_dimensions.x/(agent_radius*2)+0.5),
                              unsigned(real_dimensions.y/(agent_radius*2)+0.5),
                              unsigned(real_dimensions.z/(agent_radius*2)+0.5));
  logfile_.write((char*)(&dimensions.x), sizeof(dimensions.x));
  logfile_.write((char*)(&dimensions.y), sizeof(dimensions.y));
  logfile_.write((char*)(&dimensions.z), sizeof(dimensions.z));
  Vector<double> mod(get_model().get_space_compartment().get_micro_space().
                    get_dimensions_mod()/
                    Vector<double>(agent_radius*2, agent_radius*2, agent_radius*2));
  logfile_.write((char*)(&mod), sizeof(mod));
  const double voxel_radius(agent_radius);
  const Vector<double> min_point(0,0,0);
  logfile_.write((char*)(&min_point), sizeof(min_point));
  const Vector<double> max_point(
    space_compartment_.get_dimensions()/
    Vector<double>(agent_radius*2, agent_radius*2, agent_radius*2));
  logfile_.write((char*)(&max_point), sizeof(max_point));
  const unsigned latticeSpSize(0);
  logfile_.write((char*)(&latticeSpSize), sizeof(latticeSpSize));
  const unsigned polymerSize(0);
  logfile_.write((char*)(&polymerSize), sizeof(polymerSize));
  const unsigned reservedSize(0);
  logfile_.write((char*)(&reservedSize), sizeof(reservedSize));
  const unsigned offLatticeSpSize(species_.size());
  logfile_.write((char*)(&offLatticeSpSize), sizeof(offLatticeSpSize));
  logfile_.write((char*)(&marker_), sizeof(marker_));
  logfile_.write((char*)(&voxel_radius), sizeof(voxel_radius));
  for(unsigned i(0); i != species_.size(); ++i) {
    const unsigned stringSize(species_[i]->get_name_id().size());
    logfile_.write((char*)(&stringSize), sizeof(stringSize));
    logfile_.write(species_[i]->get_name_id().c_str(), stringSize);
    MesoSpecies* meso_species(dynamic_cast<MesoSpecies*>(species_[i]));
    double molecule_radius(voxel_radius);
    unsigned type(0);
    unsigned parent_species_id(species_.size());
    if (meso_species) { 
      molecule_radius = meso_species->get_meso_space().get_voxel_length()/2;
      type = 1; //meso species
    } else { 
      VolumeSpecies* volume_species(dynamic_cast<VolumeSpecies*>(species_[i]));
      if (volume_species) {
        type = 2; //volume_species
      } else {
        PointSpecies* point_species(dynamic_cast<PointSpecies*>(species_[i]));
        if (point_species) {
          std::vector<Species*>::iterator iterator_idx(
                 std::find(species_.begin(), species_.end(),
                           &(point_species->get_compartment())));
          if (iterator_idx == species_.end()) {
            std::cout << "In VisualWriter, must add the parent species of " <<
              point_species->get_name_id() << ", which is " <<
              point_species->get_compartment().get_name_id() << std::endl;
            exit(0);
          } 
          parent_species_id = iterator_idx-species_.begin();
          type = 3;
        }
      }
    }
    logfile_.write((char*)(&molecule_radius), sizeof(molecule_radius));
    logfile_.write((char*)(&type), sizeof(type));
    logfile_.write((char*)(&parent_species_id), sizeof(parent_species_id));
  }
}

void VisualWriter::log_structure_species() {
  const double currentTime(get_stepper().get_time());
  logfile_.write((char*)(&currentTime), sizeof(currentTime));
  logfile_.write((char*)(&marker_), sizeof(marker_));
  logfile_.write((char*)(&marker_), sizeof(marker_));
}

void VisualWriter::log_species() {
  const double currentTime(get_stepper().get_time());
  logfile_.write((char*)(&currentTime), sizeof(currentTime));
  logfile_.write((char*)(&marker_), sizeof(marker_));
  for(unsigned i(0); i != species_.size(); ++i) {
    MicroSpecies* micro_species(dynamic_cast<MicroSpecies*>(species_[i]));
    if (micro_species) {
      VolumeSpecies* volume_species(
                        dynamic_cast<VolumeSpecies*>(micro_species));
      if (volume_species) {
        log_points(*volume_species, i);
      } else {
        PointSpecies* point_species(dynamic_cast<PointSpecies*>(micro_species));
        if (point_species) {
          log_points(*point_species, i);
        }
      }
      if (write_xdmf_) {
        log_xdmf(*micro_species);
      }
    } else {
      MesoSpecies* meso_species(dynamic_cast<MesoSpecies*>(species_[i]));
      if (meso_species) {
        log_points(*meso_species, i);
        if (write_xdmf_) {
          log_xdmf(*meso_species);
        }
      } else {
        std::cout << "Unknown species type in visual logger" << std::endl;
      }
    }
  }
  logfile_.write((char*)(&marker_), sizeof(marker_));
}

void VisualWriter::log_points(VolumeSpecies& volume_species,
                              const unsigned index) {
  logfile_.write((char*)(&index), sizeof(index));
  int size(-1);
  const std::vector<Vector<double>>&
    agents(volume_species.get_relative_positions());
  if (volume_species.relative_positions_changed() || 
      volume_species.size_changed()) {
    size = agents.size();
    volume_species.set_relative_positions_changed(false);
    volume_species.set_size_changed(false);
  }
  logfile_.write((char*)(&size), sizeof(size)); 
  if (size >= 0) {
    const double agent_radius(volume_species.get_agent_radius());
    for(int i(0); i < size; ++i) {
      const Vector<double>& agent(agents[i]/
                                 Vector<double>(agent_radius*2,
                                               agent_radius*2,
                                               agent_radius*2));
      logfile_.write((char*)(&agent), sizeof(agent));
    }
    const std::vector<unsigned>& agent_cids(volume_species.get_agent_cids());
    logfile_.write(reinterpret_cast<const char*>(&agent_cids[0]),
                   sizeof(unsigned)*size);
  }
}  

void VisualWriter::log_points(PointSpecies& point_species,
                              const unsigned index) {
  logfile_.write((char*)(&index), sizeof(index));
  int size(-1);
  const std::vector<Vector<double>>& 
    agents(point_species.get_relative_positions());
  if (point_species.relative_positions_changed() || 
      point_species.size_changed()) {
    size = agents.size();
    point_species.set_relative_positions_changed(false);
    point_species.set_size_changed(false);
  }
  logfile_.write((char*)(&size), sizeof(size)); 
  if (size >= 0) { 
    const double agent_radius(
                     point_species.get_compartment().get_agent_radius());
    for(int i(0); i < size; ++i) {
      const Vector<double>& agent(agents[i]/
                                 Vector<double>(agent_radius*2,
                                               agent_radius*2,
                                               agent_radius*2));
      logfile_.write((char*)(&agent), sizeof(agent));
    }
    const std::vector<unsigned>& compartment_ids(
                                         point_species.get_compartment_ids());
    size = compartment_ids.size();
    logfile_.write((char*)(&size), sizeof(size)); 
    logfile_.write(reinterpret_cast<const char*>(&compartment_ids[0]),
                   sizeof(unsigned)*size);
    const std::vector<unsigned>& 
      agents_sizes(point_species.get_compartment_points_sizes());
    logfile_.write(reinterpret_cast<const char*>(&agents_sizes[0]),
                   sizeof(unsigned)*size);
  }
}

void VisualWriter::log_points(MesoSpecies& meso_species, const unsigned index) {
  logfile_.write((char*)(&index), sizeof(index));
  const std::vector<unsigned>& sizes(meso_species.get_sizes());
  const std::vector<unsigned>& 
    occupied_voxels(meso_species.get_occupied_voxels());
  int size(occupied_voxels.size());
  logfile_.write((char*)(&size), sizeof(size)); 
  std::vector<unsigned> factor_sizes;
  const double agent_radius(space_compartment_.get_agent_radius());
  for (int i(0); i < size; ++i) {
    const unsigned voxel_id(occupied_voxels[i]);
    factor_sizes.push_back(sizes[voxel_id]);
    Vector<double> center(
         meso_species.get_meso_space().voxel_id_to_xyz(voxel_id)/
         Vector<double>(agent_radius*2, agent_radius*2, agent_radius*2));
    logfile_.write((char*)(&center), sizeof(center));
  }
  logfile_.write(reinterpret_cast<char*>(&factor_sizes[0]),
                 sizeof(unsigned)*size);
}

/*
void VisualWriter::log_points(MesoSpecies& meso_species, const unsigned index) {
  logfile_.write((char*)(&index), sizeof(index));
  const std::vector<unsigned>& sizes(meso_species.get_sizes());
  const std::vector<unsigned>& 
    occupied_voxels(meso_species.get_occupied_voxels());
  unsigned size(0);
  for (unsigned i(0); i < occupied_voxels.size(); ++i) {
    const unsigned voxel_id(occupied_voxels[i]);
    size += sizes[voxel_id];
  }
  logfile_.write((char*)(&size), sizeof(size)); 
  for (unsigned i(0); i < occupied_voxels.size(); ++i) {
    const unsigned voxel_id(occupied_voxels[i]);
    const unsigned factor_size(sizes[voxel_id]);
    const double half(meso_species.get_meso_space().get_voxel_length()/2);
    Vector<double> center(
                       meso_species.get_meso_space().voxel_id_to_xyz(voxel_id));
    Vector<double> min(center-half);
    Vector<double> max(center+half);
    std::uniform_real_distribution<> uni_x(min.x, max.x);
    std::uniform_real_distribution<> uni_y(min.y, max.y);
    std::uniform_real_distribution<> uni_z(min.z, max.z);
    for (unsigned i(0); i != factor_size; ++i) {
      Vector<double> 
        factor(uni_x(space_compartment_.get_model().get_rng())/(agent_radius*2),
               uni_y(space_compartment_.get_model().get_rng())/(agent_radius*2),
               uni_z(space_compartment_.get_model().get_rng())/(agent_radius*2));
      logfile_.write((char*)(&factor), sizeof(factor));
    }
  }
}  
*/

void VisualWriter::log_xdmf(MicroSpecies& micro_species) {
  const std::vector<Vector<double>>& agents(micro_species.get_xyz_list());
  micro_xdmf_file_.write(micro_species.get_name(), agents,
                   get_stepper().get_time());
}

void VisualWriter::log_xdmf(MesoSpecies& meso_species) {
  const std::vector<unsigned>& sizes(meso_species.get_sizes());
  const std::vector<std::vector<double>>& mesh(meso_species.get_rect_mesh());
  const bool cell_centred(false);
  meso_xdmf_file_.write(mesh, cell_centred, meso_species.get_name(), sizes,
                        get_stepper().get_time());
}

