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

#include <Species.hpp>
#include <PointSpecies.hpp>
#include <Reaction.hpp>
#include <Model.hpp>
#include <Stepper.hpp>


Reaction::Reaction(PointSpecies& Ap, const double k, MesoSpecies& Cm):
  Process("Reaction::Ap_to_Cm", Ap.get_model()),
  k_(k),
  Ap_(&Ap),
  Cm_(&Cm),
  react_method_(&Reaction::Ap_to_Cm) {
    set_queue_id(get_stepper().get_process_queue().push(this));
  }

//exocytosis maxed:
Reaction::Reaction(PointSpecies& Ap, const double k, MesoSpecies& Cm,
                   const unsigned max):
  Process("Reaction::Ap_to_Cm_max", Ap.get_model()),
  k_(k),
  Ap_(&Ap),
  Cm_(&Cm),
  max_(max),
  react_method_(&Reaction::Ap_to_Cm_max) {
    set_queue_id(get_stepper().get_process_queue().push(this));
  }

Reaction::Reaction(MesoSpecies& Cm, const double k):
  Process("Reaction::Cm_to_0", Cm.get_model(), 0.1e-6),
  k_(k),
  Cm_(&Cm),
  react_method_(&Reaction::Cm_to_0) {
    set_queue_id(get_stepper().get_process_queue().push(this));
  }

Reaction::Reaction(PointSpecies& Ap, MesoSpecies& Cm, const double k):
  Process("Reaction::Ap_Cm_to_direction", Ap.get_model(), 0.1e-6),
  k_(k),
  Ap_(&Ap),
  Cm_(&Cm),
  react_method_(&Reaction::Ap_Cm_to_direction) {
    set_queue_id(get_stepper().get_process_queue().push(this));
    //the line below is not required, for BD agents should automatically by 
    //redirected according to agent_redirections and
    //agent_redirection_magnitudes:
    //Ap_->get_compartment().get_walker().set_directed_displacement();
  }

Reaction::Reaction(PointSpecies& Ap, MesoSpecies& Bm, const double k,
                   MesoSpecies& Cm):
  Process("Reaction::Ap_Bm_to_Cm", Ap.get_model(), 0.1e-6),
  k_(k),
  Ap_(&Ap),
  Bm_(&Bm),
  Cm_(&Cm),
  react_method_(&Reaction::Ap_Bm_to_Cm) {
    set_queue_id(get_stepper().get_process_queue().push(this));
  }

//Apv + Bmv -> Cm (track overlap reaction)
Reaction::Reaction(std::vector<PointSpecies*>& Apv,
                   std::vector<MesoSpecies*>& Bmv, const double k,
                   MesoSpecies& Cm):
  Process("Reaction::Apv_Bmv_to_Cm", Cm.get_model(), 0.1e-6),
  k_(k),
  Apv_(&Apv),
  Bmv_(&Bmv),
  Cm_(&Cm),
  react_method_(&Reaction::Apv_Bmv_to_Cm) {
    for (unsigned i(0); i < Bmv.size(); ++i) {
      Bmv_species_cids_.push_back(Bmv[i]->get_species_cid());
    }
    set_queue_id(get_stepper().get_process_queue().push(this));
  }

// Ap + Bm -> Cp (scavenging)
Reaction::Reaction(PointSpecies& Ap, MesoSpecies& Bm, const double k,
                   PointSpecies& Cp):
  Process("Reaction::Ap_Bm_to_Cp", Ap.get_model()),
  k_(k),
  Ap_(&Ap),
  Bm_(&Bm),
  Cp_(&Cp),
  react_method_(&Reaction::Ap_Bm_to_Cp) {
    set_queue_id(get_stepper().get_process_queue().push(this));
  }

void Reaction::initialize() {
  if (react_method_ == &Reaction::Ap_Bm_to_Cm ||
      react_method_ == &Reaction::Ap_Cm_to_direction) {
    (this->*react_method_)();
  }
}

double Reaction::step() {
  (this->*react_method_)();
  return Process::step();
}

//degradation
bool Reaction::Cm_to_0() {
  const unsigned size(Cm_->get_size()*k_*get_interval());
  MesoSpace& meso_space(Cm_->get_meso_space());
  if (size) {
    meso_space.remove_factors(Cm_->get_species_cid(), *Cm_, size);
  }
  return true;
}

//scavenging
bool Reaction::Ap_Bm_to_Cp() {
  if (k_ > 0) { 
    std::vector<Vector<double>>& positions(Ap_->get_positions());
    MesoSpace& meso_space(Bm_->get_meso_space());
    for (unsigned i(0); i != positions.size(); ++i) {
      unsigned factor_size(meso_space.get_position_factor(
                      Bm_->get_species_cid(), *Bm_, positions[i]));
      if (factor_size) { 
        meso_space.remove_factors(Bm_->get_species_cid(), *Bm_, positions[i],
                                  factor_size);
      }
    }
    meso_space.update_curr_occupied_voxels(Bm_->get_species_cid());
  }
  return true;
}

//exocytosis
bool Reaction::Ap_to_Cm() {
  if (k_ > 0) { 
    std::vector<Vector<double>>& positions(Ap_->get_positions());
    MesoSpace& meso_space(Cm_->get_meso_space());
    for (unsigned i(0); i != positions.size(); ++i) {
      meso_space.add_factor(Cm_->get_species_cid(), *Cm_,
                            positions[i]);
    }
  }
  return true;
}

//exocytosis maxed
bool Reaction::Ap_to_Cm_max() {
  if (k_ > 0) { 
    std::vector<Vector<double>>& positions(Ap_->get_positions());
    MesoSpace& meso_space(Cm_->get_meso_space());
    for (unsigned i(0); i != positions.size(); ++i) {
      meso_space.add_factor(Cm_->get_species_cid(), *Cm_,
                            positions[i], max_);
    }
  }
  return true;
}

void Reaction::set_rate(const double k) {
  k_ = k;
}

//Apv + Bmv -> Cm (track overlap reaction)
bool Reaction::Apv_Bmv_to_Cm() {
  if (k_ > 0) { 
    MesoSpace& meso_space(Cm_->get_meso_space());
    for (unsigned i(0); i < Apv_->size(); ++i) {
      PointSpecies& Ap(*(*Apv_)[i]);
      std::vector<Vector<double>>& positions(Ap.get_positions());
      for (unsigned j(0); j < positions.size(); ++j) {
        meso_space.replace_factors(Bmv_species_cids_, *Bmv_, 
                                   Cm_->get_species_cid(), *Cm_, positions[j],
                                   i);
      }
    }
    for (unsigned i(0); i < Bmv_species_cids_.size(); ++i) {
      meso_space.update_curr_occupied_voxels(Bmv_species_cids_[i]);
    }
  }
  return true;
}

bool Reaction::Ap_Bm_to_Cm() {
  if (k_ > 0) { 
    std::vector<Vector<double>>& positions(Ap_->get_positions());
    MesoSpace& meso_space(Bm_->get_meso_space());
    for (unsigned i(0); i != positions.size(); ++i) {
      unsigned factor_size(meso_space.get_position_factor(
                      Bm_->get_species_cid(), *Bm_, positions[i]));
      if (factor_size > k_) { 
        meso_space.add_factor(Cm_->get_species_cid(), *Cm_, positions[i]);
      }
    }
  }
  return true;
}

//chemotaxis
bool Reaction::Ap_Cm_to_direction() {
  const std::vector<Vector<double>>& direction_vectors(
                           Ap_->get_direction_vectors());
  const unsigned agent_point_size(Ap_->size());
  std::vector<Vector<double>>& point_positions(Ap_->get_positions());
  MesoSpace& meso_space(Cm_->get_meso_space());
  std::vector<unsigned> agent_factors(Ap_->get_compartment().size(), 0);
  std::vector<unsigned> all_factors(meso_space.get_position_factors(
                                    Cm_->get_species_cid(),
                                    point_positions, agent_factors,
                                    agent_point_size));
  std::vector<Vector<double>>& agent_redirections(
            Ap_->get_compartment().get_walker().get_agent_redirections());
  std::vector<unsigned>& agent_redirection_magnitudes(
      Ap_->get_compartment().get_walker().get_agent_redirection_magnitudes());
  agent_redirections.resize(agent_factors.size());
  agent_redirection_magnitudes.resize(agent_factors.size());
  for (unsigned i(0); i != agent_factors.size(); ++i) {
    const unsigned k(i*agent_point_size);
    auto [min, max] = std::minmax_element(all_factors.begin()+k,
                                      all_factors.begin()+k+agent_point_size);
    const unsigned max_id(std::distance(all_factors.begin(), max));
    const unsigned max_val(all_factors[max_id]);
    agent_redirection_magnitudes[i] = 0;
    if (max_val > 0) {
      std::vector<unsigned> max_ids;
      for (unsigned j(0); j != agent_point_size; ++j) {
        if (all_factors[k+j] == max_val) {
          max_ids.push_back(k+j);
        }
      }
      std::uniform_int_distribution<> uni_dist(0, max_ids.size()-1);
      const unsigned selected_id(uni_dist(get_model().get_rng()));
      agent_redirections[i] = direction_vectors[max_ids[selected_id]];
      if (k_ < 0) {
        agent_redirections[i] *= -1;
      }
      const unsigned min_id(std::distance(all_factors.begin(), min));
      const unsigned min_val(all_factors[min_id]);
      const unsigned magnitude(fabs(max_val-min_val));
      if (magnitude >= fabs(k_)) {
        agent_redirection_magnitudes[i] = 1;
      }
    } else {
      //move in random direction
      std::uniform_int_distribution<> uni_dist(0, agent_factors.size()-1);
      const unsigned selected_id(uni_dist(get_model().get_rng()));
      agent_redirections[i] = direction_vectors[selected_id];
    }
  }
  return true;
}

/*
//chemotaxis
bool Reaction::Ap_Cm_to_direction() {
  const std::vector<Vector<double>>& direction_vectors(
                           Ap_->get_direction_vectors());
  const unsigned agent_point_size(Ap_->size());
  std::vector<Vector<double>>& point_positions(Ap_->get_positions());
  MesoSpace& meso_space(Cm_->get_meso_space());
  std::vector<unsigned> agent_factors(Ap_->get_compartment().size(), 0);
  std::vector<unsigned> all_factors(meso_space.get_position_factors(
                                    Cm_->get_species_cid(),
                                    point_positions, agent_factors,
                                    agent_point_size));
  std::vector<Vector<double>>& agent_redirections(
            Ap_->get_compartment().get_walker().get_agent_redirections());
  std::vector<unsigned>& agent_redirection_magnitudes(
      Ap_->get_compartment().get_walker().get_agent_redirection_magnitudes());
  agent_redirections.resize(agent_factors.size());
  agent_redirection_magnitudes.resize(agent_factors.size());
  for (unsigned i(0); i != agent_factors.size(); ++i) {
    const unsigned k(i*agent_point_size);
    auto [min, max] = std::minmax_element(all_factors.begin()+k,
                                      all_factors.begin()+k+agent_point_size);
    if (k_ < 0) {
      std::swap(min, max);
    }
    const unsigned max_id(std::distance(all_factors.begin(), max));
    const unsigned max_val(all_factors[max_id]);
    std::vector<unsigned> max_ids;
    for (unsigned j(0); j != agent_point_size; ++j) {
      //if (all_factors[k+j] == max_val || max_val < fabs(k_)) {
      if (all_factors[k+j] == max_val) {
        max_ids.push_back(k+j);
      }
    }
    std::uniform_int_distribution<> uni_dist(0, max_ids.size()-1);
    const unsigned selected_id(uni_dist(get_model().get_rng()));
    agent_redirections[i] = direction_vectors[max_ids[selected_id]];
    if (max_ids.size()) {
      std::cout << "size:" << max_ids.size() << " val:" << max_val << std::endl;
      const unsigned min_id(std::distance(all_factors.begin(), min));
      const unsigned min_val(all_factors[min_id]);
      const unsigned magnitude(fabs(max_val-min_val));
      agent_redirection_magnitudes[i] = magnitude;
    } else {
      agent_redirection_magnitudes[i] = 0;
    }
  }
  return true;
}
//chemotaxis
bool Reaction::Ap_Cm_to_direction() {
  const std::vector<Vector<double>>& direction_vectors(
                           Ap_->get_direction_vectors());
  const unsigned agent_point_size(Ap_->size());
  std::vector<Vector<double>>& point_positions(Ap_->get_positions());
  MesoSpace& meso_space(Cm_->get_meso_space());
  std::vector<unsigned> agent_factors(Ap_->get_compartment().size(), 0);
  std::vector<unsigned> all_factors(meso_space.get_position_factors(
                                    Cm_->get_species_cid(),
                                    point_positions, agent_factors,
                                    agent_point_size));
  std::vector<Vector<double>>& agent_redirections(
            Ap_->get_compartment().get_walker().get_agent_redirections());
  agent_redirections.resize(agent_factors.size());
  for (unsigned i(0); i != agent_factors.size(); ++i) {
    const unsigned k(i*agent_point_size);
    std::vector<unsigned>::iterator max;
    if (k_ < 0) {
      max = (std::min_element(all_factors.begin()+k,
                                   all_factors.begin()+k+agent_point_size));
    } else {
      max = (std::max_element(all_factors.begin()+k,
                                   all_factors.begin()+k+agent_point_size));
    }
    const unsigned max_id(std::distance(all_factors.begin(), max));
    const unsigned max_val(all_factors[max_id]);
    std::vector<unsigned> max_ids;
    for (unsigned j(0); j != agent_point_size; ++j) {
      if (all_factors[k+j] == max_val || max_val < fabs(k_)) {
        max_ids.push_back(k+j);
      }
    }
    std::uniform_int_distribution<> uni_dist(0, max_ids.size()-1);
    const unsigned selected_id(uni_dist(get_model().get_rng()));
    agent_redirections[i] = direction_vectors[max_ids[selected_id]];
  }
  return true;
}

bool Reaction::Ap_Cm_to_direction() {
  const std::vector<Vector<double>>& direction_vectors(
                           Ap_->get_direction_vectors());
  const unsigned agent_point_size(Ap_->size());
  std::vector<Vector<double>>& point_positions(Ap_->get_positions());
  MesoSpace& meso_space(Cm_->get_meso_space());
  std::vector<unsigned> agent_factors(Ap_->get_compartment().size(), 0);
  std::vector<unsigned> all_factors(meso_space.get_position_factors(
                                    Cm_->get_species_cid(),
                                    point_positions, agent_factors,
                                    agent_point_size));
  std::vector<Vector<double>>& agent_redirections(
            Ap_->get_compartment().get_walker().get_agent_redirections());
  agent_redirections.resize(agent_factors.size());
  for (unsigned i(0); i != agent_factors.size(); ++i) {
    const unsigned k(i*agent_point_size);
    std::vector<unsigned>::iterator max(std::max_element(all_factors.begin()+k,
                                   all_factors.begin()+k+agent_point_size));
    unsigned max_id(std::distance(all_factors.begin(), max));
    agent_redirections[i] = direction_vectors[max_id];
  }
  return true;
}
*/

/*
bool Reaction::Ap_Cm_to_direction() {
  const std::vector<Vector<double>>& direction_vectors(
                           Ap_->get_direction_vectors());
  const unsigned agent_point_size(Ap_->size());
  std::vector<Vector<double>>& point_positions(Ap_->get_positions());
  MesoSpace& meso_space(Cm_->get_meso_space());
  std::vector<unsigned> agent_factors(Ap_->get_compartment().size(), 0);
  std::vector<unsigned> all_factors(meso_space.get_position_factors(
                                    Cm_->get_species_cid(),
                                    point_positions, agent_factors,
                                    agent_point_size));
  std::vector<Vector<double>>& agent_redirections(
            Ap_->get_compartment().get_walker().get_agent_redirections());
  agent_redirections.resize(agent_factors.size());
  for (unsigned i(0); i != agent_factors.size(); ++i) {
    const unsigned k(i*agent_point_size);
    unsigned max(0);
    agent_redirections[i] = direction_vectors[k];
    for (unsigned j(0); j != agent_point_size; ++j) {
      if (all_factors[k+j] > max) {
        agent_redirections[i] = direction_vectors[k+j];
        max = all_factors[k+j];
      }
    }
  }
  return true;
}
*/

/*
bool Reaction::Ap_Cm_to_direction() {
  const std::vector<Vector<double>>& direction_vectors(
                           Ap_->get_direction_vectors());
  const unsigned agent_point_size(Ap_->size());
  std::vector<Vector<double>>& point_positions(Ap_->get_positions());
  MesoSpace& meso_space(Cm_->get_meso_space());
  std::vector<unsigned> agent_factors(Ap_->get_compartment().size(), 0);
  std::vector<unsigned> all_factors(meso_space.get_position_factors(
                                    Cm_->get_species_cid(),
                                    point_positions, agent_factors,
                                    agent_point_size));
  std::vector<Vector<double>>& agent_redirections(
            Ap_->get_compartment().get_walker().get_agent_redirections());
  agent_redirections.resize(agent_factors.size());
  for (unsigned i(0); i != agent_factors.size(); ++i) {
    std::uniform_int_distribution<> uni_dist(1,
                                             agent_factors[i]+agent_point_size);
    const unsigned k(i*agent_point_size);
    const unsigned selected_factor(uni_dist(get_model().get_rng()));
    unsigned cnt(0);
    for (unsigned j(0); j != agent_point_size; ++j) {
      if (cnt + all_factors[k+j] + 1 >= selected_factor) {
        agent_redirections[i] = direction_vectors[k+j];
        break;
      }
      cnt += all_factors[k+j] + 1;
    }
  }
  return true;
}
*/

