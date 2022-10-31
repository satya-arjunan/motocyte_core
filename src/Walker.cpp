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
#include <limits>
#include <iomanip>
#include <Species.hpp>
#include <Walker.hpp>
#include <Model.hpp>


Walker::Walker(const double D, Model& model, MotileSpecies& motile_species,
               const double interval):
  //set the interval for default discrete walk on lattice voxels:
  Process(std::string("Walker "+motile_species.get_name()), model, interval),
  D_(D),
  motile_species_(motile_species) {
    set_queue_id(get_stepper().get_process_queue().push(this));
    std::stringstream filename(
     //when changing the filename, don't forget to change the window size in
     //Walker.hpp!! set window_size = 1 if no windows used
     //"simulation_parameters_exp_model_d_vs_t.csv");
     //"simulation_parameters_exp_model_d_vs_s.csv");
     //"simulation_parameters_exp_model_c_vs_t.csv");
     //"simulation_parameters_exp_model_c_vs_s.csv");
     //"simulation_parameters_exp_model_b_vs_t.csv");
     //"simulation_parameters_exp_model_b_vs_s.csv");
     "simulation_parameters_exp_model_a.csv");
     //"simulation_parameters_exp_point3_window_5_vmf2_vonmises_intercept.csv");
     //"simulation_parameters_exp_window_20_vmf2_vonmises_intercept.csv");
     //"simulation_parameters_synthetic_window_20_walk_3D_vmf2_intercept.csv");
     //"simulation_parameters_window_5_vmf2_vonmises_intercept_exp.csv");
     //"simulation_parameters_window_20_vmf2_intercept_exp.csv");
     //"simulation_parameters_window_20_vmf2_intercept_exp.csv");
     //"simulation_parameters_window_5_vmf2_intercept_3D_walk.csv");
     //"simulation_parameters_2vmf_intercept_3D_walk.csv");
     //"simulation_parameters_directed.csv");
     //"simulation_parameters_corkscrew.csv");
     //"simulation_parameters_t_vs_s_vt_vr.csv");
     //"simulation_parameters_instant_speed_vs_s_t_s2_gamma_roll_vm.csv");
     //"simulation_parameters_instant_speed_vs_s_t_s2_gamma_roll_all_vm.csv");
           //_"simulation_parameters_instant_speed_vs_s_t_s2_gamma.csv");
    param_file_.open(filename.str().c_str(), std::ios::in);
    parse_file();
  }

void Walker::add_agent() {
  initialize_agent_parameters();
  agent_redirection_magnitudes_.push_back(0);
  agent_redirections_.push_back(Vector<double>(0, 0, 0));
  Vector<double> agent_direction(get_random_direction());
  //Vector<double> agent_direction(0.979583720515923,0.154094974946487,
  //                              0.129114186658261);
  agent_direction = agent_direction.norm();
  std::cout << "x:" << agent_direction.x << " " << agent_direction.y << " " <<
    agent_direction.z << std::endl;
  facing_.push_back(agent_direction);
  agent_orientations_.push_back(Quaternion::face_vector(agent_direction));
  agent_prev_speeds_.push_back(std::vector<double>());
  agent_prev_turns_.push_back(std::vector<double>());
  agent_prev_rolls_.push_back(std::vector<double>());
  agent_prev_xyz_.push_back(std::vector<Vector<double>>());
  agent_mean_speed_.push_back(0);
  agent_speed_cnt_.push_back(0);
}

void Walker::set_walk_type(const unsigned walk_type) {
  walk_type_ = walk_type;
}

void Walker::parse_file() {
  // read csv header
  std::string line;
  std::getline(param_file_, line);
  std::stringstream str(line); 
  std::string word;
  int col(-1);
  unsigned agent_params_size(0);
  while (std::getline(str, word, ',')) { 
    if (word == "instant_speed") {
      instant_speed_ = col;
      agent_params_size++;
    }
    else if (word == "mean_speed" || word == "windowed_mean_speed") {
      mean_speed_ = col;
      agent_params_size++;
    }
    else if (word == "std_mean_speed") {
      std_mean_speed_ = col;
      agent_params_size++;
    }
    else if (word == "mean_pitch_rate") {
      mean_pitch_rate_ = col;
      agent_params_size++;
    }
    else if (word == "std_mean_pitch_rate") {
      std_mean_pitch_rate_ = col;
      agent_params_size++;
    }
    else if (word == "mean_roll_rate") {
      mean_roll_rate_ = col;
      agent_params_size++;
    }
    else if (word == "std_mean_roll_rate") {
      std_mean_roll_rate_ = col;
      agent_params_size++;
    }
    else {
      std::string subword;
      std::stringstream substr(word); 
      std::getline(substr, subword, ' ');
      if (subword == "next_speed") {
        if (next_speed_ < 0) {
          agent_params_size++;
          next_speed_ = col;
        }
      }
      if (subword == "next_direction") {
        if (next_direction_ < 0) {
          agent_params_size++;
          next_direction_ = col;
        }
      }
      if (subword == "next_roll") {
        if (next_roll_ < 0) {
          agent_params_size++;
          next_roll_ = col;
        }
      }
      else if (subword == "instant_pitch_rate") {
        if (next_pitch_ < 0) {
          agent_params_size++;
          next_pitch_ = col;
        }
      }
      else if (subword == "instant_roll_rate") {
        if (next_roll_ < 0) {
          agent_params_size++;
          next_roll_ = col;
        }
      }
      else if (subword == "std_mean_speed") {
        if (std_at_mean_speed_ < 0) {
          agent_params_size++;
          std_at_mean_speed_ = col;
          if (std_mean_speed_ < 0) {
            std_mean_speed_ = std_at_mean_speed_;
          }
        }
        std::getline(substr, subword);
        unsigned size(std::stoi(subword));
        std_at_mean_speed_size_ = std::max(std_at_mean_speed_size_, size);
      }
    }
    col++;
  }

  const unsigned param_col_size(col);
  params_.resize(max_+1);
  unsigned row(0);
  while (std::getline(param_file_, line)) {
    str.clear();
    str.str(line);
    // skip first column
    std::getline(str, word, ',');
    col = 0;
    params_[row].resize(param_col_size);
    while (std::getline(str, word, ',')) {
      params_[row][col] = std::stof(word);
      ++col;
    }
    ++row;
  }
  if (mean_speed_ >= 0) {
    speed_max_ = params_[max_][mean_speed_];
  }
  else if (instant_speed_ >= 0) {
    speed_max_ = params_[max_][instant_speed_];
  }
  agent_params_.resize(agent_params_size);
}

double Walker::draw_stable_distribution_at_col(const unsigned col) {
  const double alpha(params_[alpha_][col]);
  const double beta(params_[beta_][col]);
  const double loc(params_[loc_][col]);
  const double scale(params_[sigma_][col]/sqrt(2));
  return loc + gsl_ran_levy_skew(get_model().get_gsl_rng(), scale, alpha, beta);
}

double Walker::draw_gamma_distribution_at_col(const int col,
                                             const unsigned agent_id) {
  if (col == std_at_mean_speed_) {
    const double mean_speed(agent_params_[mean_speed_][agent_id]);
    return get_new_std_mean_speed(agent_id, mean_speed);
  }
  else if (col == instant_speed_) {
    const double alpha(params_[alpha_][col]);
    const double loc(params_[loc_][col]);
    const double scale(params_[sigma_][col]);
    double max(params_[max_][col]);
    if (max == 0) {
      max = std::numeric_limits<double>::infinity();
    }
    double value;
    do {
      value = loc + gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale);
      /*
      std::cout << "value:" << value << " loc:" << loc << " alpha:" << alpha <<
        " scale:" << scale << std::endl;
        */
    }
    while (value < 0 || value >= max);
    return value;
  }
  return 0;
}


void Walker::initialize_agent_parameters() {
  for (unsigned i(0); i < agent_params_.size(); ++i) { 
    const unsigned agent_id(agent_params_[i].size());
    const double value(draw_gamma_distribution_at_col(i, agent_id));
    /*std::cout << "agent_id:" << agent_id << " param id:" << i << " value:" <<
      value << std::endl;
      */
    agent_params_[i].push_back(value);
  }
}

void Walker::remove_agent(const unsigned agent_id) {
  for (unsigned i(0); i < agent_params_.size(); ++i) { 
    agent_params_[i][agent_id] = agent_params_[i].back();
    agent_params_[i].pop_back();
  }
  agent_redirection_magnitudes_[agent_id] = 
    agent_redirection_magnitudes_.back();
  agent_redirection_magnitudes_.pop_back();

  agent_redirections_[agent_id] = agent_redirections_.back();
  agent_redirections_.pop_back();

  agent_orientations_[agent_id] = agent_orientations_.back();
  agent_orientations_.pop_back();

  agent_prev_speeds_[agent_id] = agent_prev_speeds_.back();
  agent_prev_speeds_.pop_back();

  agent_prev_turns_[agent_id] = agent_prev_turns_.back();
  agent_prev_turns_.pop_back();

  agent_prev_rolls_[agent_id] = agent_prev_rolls_.back();
  agent_prev_rolls_.pop_back();

  agent_mean_speed_[agent_id] = agent_mean_speed_.back();
  agent_mean_speed_.pop_back();

  agent_speed_cnt_[agent_id] = agent_speed_cnt_.back();
  agent_speed_cnt_.pop_back();
}

//called after populating agents. So we now know the size of agents.
void Walker::initialize() {
  const unsigned size(motile_species_.size());
  for (unsigned i(0); i < size; ++i) {
    add_agent();
  }
}

void Walker::set_interval(const double interval) {
  Process::set_interval(interval);
}


/*
//original
const Vector<double> Walker::get_displacement(const unsigned agent_id) {
  double agent_speed = get_new_agent_speed(agent_id);
  Quaternion orientation(get_updated_agent_orientation(agent_id, agent_speed));
  const Vector<double> x_axis(1, 0, 0);
  Vector<double> facing(orientation.transform(x_axis));		
  double disp(agent_speed*get_interval()/60.0); //convert interval in s to min
  disp = fabs(disp);
  Vector<double> nex(double(facing.x*disp), double(facing.y*disp),
                     double(facing.z*disp));			
  return nex;
}
*/

const Vector<double> Walker::get_displacement(const Vector<double> agent_xyz,
                                              const unsigned agent_id) {
  switch (walk_type_) {
  case exp_model_d_vs_t_:
    return get_exp_model_d_vs_t_displacement(agent_xyz, agent_id);
  case exp_model_d_vs_s_:
    return get_exp_model_d_vs_s_displacement(agent_xyz, agent_id);
  case exp_model_c_vs_t_:
    return get_exp_model_c_vs_t_displacement(agent_xyz, agent_id);
  case exp_model_c_vs_s_:
    return get_exp_model_c_vs_s_displacement(agent_xyz, agent_id);
  case exp_model_b_vs_t_:
    return get_exp_model_b_vs_t_displacement(agent_xyz, agent_id);
  case exp_model_b_vs_s_:
    return get_exp_model_b_vs_s_displacement(agent_xyz, agent_id);
  case exp_model_a_:
    return get_exp_model_a_displacement(agent_xyz, agent_id);
  case synthetic_window_walk_3D_fit_:
    return get_synthetic_window_walk_3D_fit_displacement(agent_xyz, agent_id);
  case synthetic_window_walk_3D_:
    return get_synthetic_window_walk_3D_displacement(agent_xyz, agent_id);
  case exp_window_vmf2_intercept_fit_:
    return get_exp_window_vmf2_intercept_fit_displacement(agent_xyz, agent_id);
  case exp_window_vmf2_von_mises_intercept_fit_:
    return get_exp_window_vmf2_von_mises_intercept_fit_displacement(agent_xyz,
                                                                agent_id);
  case window_vmf2_intercept_walk_3D_:
    return get_window_vmf2_intercept_walk_3D_displacement(agent_xyz, agent_id);
  case vmf2_intercept_walk_3D_:
    return get_vmf2_intercept_walk_3D_displacement(agent_id);
  case walk_3D_:
    return get_walk_3D_displacement(agent_id);
  case walk_2D_:
    return get_walk_2D_displacement(agent_id);
  case brownian_walk_:
    return get_brownian_walk_displacement(agent_id);
  case brownian_uniform_turn_:
    return get_brownian_uniform_turn_displacement(agent_xyz, agent_id);
  case uniform_directed_walk_:
    return get_uniform_directed_walk_displacement(agent_id);
  case random_directed_walk_:
    return get_random_directed_walk_displacement(agent_id);
    break;
  }
  return Vector<double>(0, 0, 0);
}

void Walker::get_vectors_v_and_y(const Vector<double> agent_xyz,
                                 const unsigned agent_id,
                                 Vector<double>& v,
                                 Vector<double>& y) {
  std::vector<Vector<double>>& prev_xyz(agent_prev_xyz_[agent_id]);
  if (!prev_xyz.size()) {
    for (unsigned i(0); i < 2; ++i) {
      double x(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double z(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double y(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      prev_xyz.push_back(Vector<double>(x, y, z));
    }
  }
  update_new_agent_xyz(prev_xyz, agent_xyz); // now, we have at least 3 points
  const Vector<double> pu(prev_xyz[0]); //earliest point
  const Vector<double> pv(prev_xyz[window_size_]);
  const Vector<double> pw(prev_xyz[prev_xyz.size()-1]); //latest point
  const Vector<double> u = (pv-pu).norm();
  v = (pw-pv).norm();
  y = u.cross_product(v);
  y = y.norm();
}

const Vector<double>
Walker::get_exp_model_a_displacement(const Vector<double> agent_xyz,
                                     const unsigned agent_id) {
  const double agent_speed(draw_next_speed_vs_none_gamma(agent_id));
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  Vector<double> d(draw_next_direction_vs_none_vmf());
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(draw_next_roll_vs_none_von_mises());
  Vector<double> v, y;
  get_vectors_v_and_y(agent_xyz, agent_id, v, y);
  
  const Quaternion y_nod(Quaternion(0, y).passive_rotate(roll, v));
  const Quaternion w(Quaternion(0, v).passive_rotate(turn, y_nod.v_));
  const double disp(agent_speed*get_interval()/60.0);
  return w.v_*disp;
}

const Vector<double>
Walker::get_exp_model_b_vs_s_displacement(const Vector<double> agent_xyz,
                                          const unsigned agent_id) {
  const double agent_speed(draw_next_speed_vs_none_gamma(agent_id));
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  Vector<double> d(draw_next_direction_vs_s_vmf(agent_speed));
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(draw_next_roll_vs_s_von_mises(agent_speed));
  Vector<double> v, y;
  get_vectors_v_and_y(agent_xyz, agent_id, v, y);
  
  const Quaternion y_nod(Quaternion(0, y).passive_rotate(roll, v));
  const Quaternion w(Quaternion(0, v).passive_rotate(turn, y_nod.v_));
  const double disp(agent_speed*get_interval()/60.0);
  return w.v_*disp;
}

const Vector<double>
Walker::get_exp_model_b_vs_t_displacement(const Vector<double> agent_xyz,
                                          const unsigned agent_id) {
  const double agent_speed(draw_next_speed_vs_none_gamma(agent_id));
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  Vector<double> d(draw_next_direction_vs_s_vmf(agent_speed));
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(draw_next_roll_vs_t_von_mises(turn));
  Vector<double> v, y;
  get_vectors_v_and_y(agent_xyz, agent_id, v, y);
  
  const Quaternion y_nod(Quaternion(0, y).passive_rotate(roll, v));
  const Quaternion w(Quaternion(0, v).passive_rotate(turn, y_nod.v_));
  const double disp(agent_speed*get_interval()/60.0);
  return w.v_*disp;
}

const Vector<double>
Walker::get_exp_model_c_vs_s_displacement(const Vector<double> agent_xyz,
                                          const unsigned agent_id) {
  if (agent_prev_speeds_[agent_id].size() < window_size_) {
    const double new_speed(agent_params_[instant_speed_][agent_id]);
    update_new_agent_speed(agent_prev_speeds_[agent_id], new_speed);
  }
  const double prev_speed(agent_prev_speeds_[agent_id].back());
  const double agent_speed(draw_next_speed_vs_sp_gamma(agent_id, prev_speed));
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  Vector<double> d(draw_next_direction_vs_s_vmf(agent_speed));
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(draw_next_roll_vs_s_von_mises(agent_speed));
  Vector<double> v, y;
  get_vectors_v_and_y(agent_xyz, agent_id, v, y);
  
  const Quaternion y_nod(Quaternion(0, y).passive_rotate(roll, v));
  const Quaternion w(Quaternion(0, v).passive_rotate(turn, y_nod.v_));
  const double disp(agent_speed*get_interval()/60.0);
  return w.v_*disp;
}

const Vector<double>
Walker::get_exp_model_c_vs_t_displacement(const Vector<double> agent_xyz,
                                          const unsigned agent_id) {
  if (agent_prev_speeds_[agent_id].size() < window_size_) {
    const double new_speed(agent_params_[instant_speed_][agent_id]);
    update_new_agent_speed(agent_prev_speeds_[agent_id], new_speed);
  }
  const double prev_speed(agent_prev_speeds_[agent_id].back());
  const double agent_speed(draw_next_speed_vs_sp_gamma(agent_id, prev_speed));
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  Vector<double> d(draw_next_direction_vs_s_vmf(agent_speed));
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(draw_next_roll_vs_t_von_mises(turn));
  Vector<double> v, y;
  get_vectors_v_and_y(agent_xyz, agent_id, v, y);
  
  const Quaternion y_nod(Quaternion(0, y).passive_rotate(roll, v));
  const Quaternion w(Quaternion(0, v).passive_rotate(turn, y_nod.v_));
  const double disp(agent_speed*get_interval()/60.0);
  return w.v_*disp;
}

const Vector<double>
Walker::get_exp_model_d_vs_s_displacement(const Vector<double> agent_xyz,
                                          const unsigned agent_id) {
  if (agent_prev_speeds_[agent_id].size() < window_size_) {
    const double new_speed(agent_params_[instant_speed_][agent_id]);
    update_new_agent_speed(agent_prev_speeds_[agent_id], new_speed);
  }
  const double prev_speed(agent_prev_speeds_[agent_id].back());
  if (agent_prev_turns_[agent_id].size() < window_size_) {
    const double new_turn(1/prev_speed*M_PI);
    update_new_agent_speed(agent_prev_turns_[agent_id], new_turn);
  }
  const double agent_speed(draw_next_speed_vs_sp_gamma(agent_id, prev_speed));
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  const double prev_turn(agent_prev_turns_[agent_id].back());
  Vector<double> d(draw_next_direction_vs_s_tp_vmf(agent_speed, prev_turn));
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  update_new_agent_speed(agent_prev_turns_[agent_id], turn);
  const double roll(draw_next_roll_vs_s_von_mises(agent_speed));
  Vector<double> v, y;
  get_vectors_v_and_y(agent_xyz, agent_id, v, y);
  
  const Quaternion y_nod(Quaternion(0, y).passive_rotate(roll, v));
  const Quaternion w(Quaternion(0, v).passive_rotate(turn, y_nod.v_));
  const double disp(agent_speed*get_interval()/60.0);
  return w.v_*disp;
}

const Vector<double>
Walker::get_exp_model_d_vs_t_displacement(const Vector<double> agent_xyz,
                                          const unsigned agent_id) {
  if (agent_prev_speeds_[agent_id].size() < window_size_) {
    const double new_speed(agent_params_[instant_speed_][agent_id]);
    update_new_agent_speed(agent_prev_speeds_[agent_id], new_speed);
  }
  const double prev_speed(agent_prev_speeds_[agent_id].back());
  if (agent_prev_turns_[agent_id].size() < window_size_) {
    const double new_turn(1/prev_speed*M_PI);
    update_new_agent_speed(agent_prev_turns_[agent_id], new_turn);
  }
  const double agent_speed(draw_next_speed_vs_sp_gamma(agent_id, prev_speed));
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  const double prev_turn(agent_prev_turns_[agent_id].back());
  Vector<double> d(draw_next_direction_vs_s_tp_vmf(agent_speed, prev_turn));
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  update_new_agent_speed(agent_prev_turns_[agent_id], turn);
  const double roll(draw_next_roll_vs_t_von_mises(turn));
  Vector<double> v, y;
  get_vectors_v_and_y(agent_xyz, agent_id, v, y);
  
  const Quaternion y_nod(Quaternion(0, y).passive_rotate(roll, v));
  const Quaternion w(Quaternion(0, v).passive_rotate(turn, y_nod.v_));
  const double disp(agent_speed*get_interval()/60.0);
  return w.v_*disp;
}

const Vector<double>
Walker::get_synthetic_window_walk_3D_fit_displacement(
                                                 const Vector<double> agent_xyz,
                                                 const unsigned agent_id) {
  std::vector<Vector<double>>& prev_xyz(agent_prev_xyz_[agent_id]);
  if (!prev_xyz.size()) {
    for (unsigned i(0); i < 2; ++i) {
      double x(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double z(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double y(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      prev_xyz.push_back(Vector<double>(x, y, z));
    }
  }
  update_new_agent_xyz(prev_xyz, agent_xyz); // now, we have at least 3 points
  Vector<double> p1(prev_xyz[0]); //earliest point
  Vector<double> p2(prev_xyz[1]);
  Vector<double> p3(prev_xyz[prev_xyz.size()-2]);
  Vector<double> p4(prev_xyz[prev_xyz.size()-1]); //latest point
  Vector<double> u(p3-p1);
  u = u.norm();
  Vector<double> v(p4-p2);
  v = v.norm();
  Vector<double> y(u.cross_product(v));
  y = y.norm();

  Vector<double> d(draw_direction_from_vmf2_intercept());
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(atan2(d.y, d.z));

  const Quaternion yQ(Quaternion(0, y).passive_rotate(roll, v));
  const Vector<double> _y(yQ.v_);
  const Quaternion wQ(Quaternion(0, v).passive_rotate(turn, _y));
  const double agent_speed = get_new_agent_speed(agent_id);
  //convert interval in s to min
  const double disp(agent_speed*get_interval()/60.0);
  const Vector<double> w(wQ.v_*disp);
  return w;
}

const Vector<double>
Walker::get_exp_window_vmf2_intercept_fit_displacement(
                                                 const Vector<double> agent_xyz,
                                                 const unsigned agent_id) {
  std::vector<Vector<double>>& prev_xyz(agent_prev_xyz_[agent_id]);
  if (!prev_xyz.size()) {
    for (unsigned i(0); i < 2; ++i) {
      double x(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double z(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double y(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      prev_xyz.push_back(Vector<double>(x, y, z));
    }
  }
  update_new_agent_xyz(prev_xyz, agent_xyz); // now, we have at least 3 points
  Vector<double> p1(prev_xyz[0]); //earliest point
  Vector<double> p2(prev_xyz[1]);
  Vector<double> p3(prev_xyz[prev_xyz.size()-2]);
  Vector<double> p4(prev_xyz[prev_xyz.size()-1]); //latest point
  Vector<double> u(p3-p1);
  u = u.norm();
  Vector<double> v(p4-p2);
  v = v.norm();
  Vector<double> y(u.cross_product(v));
  y = y.norm();

  Vector<double> d(draw_direction_from_vmf2_intercept());
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(atan2(d.y, d.z));

  const Quaternion yQ(Quaternion(0, y).passive_rotate(roll, v));
  const Vector<double> _y(yQ.v_);
  const Quaternion wQ(Quaternion(0, v).passive_rotate(turn, _y));
  const double agent_speed = get_new_agent_speed(agent_id);
  //convert interval in s to min
  const double disp(agent_speed*get_interval()/60.0);
  const Vector<double> w(wQ.v_*disp);
  return w;
}

const Vector<double>
Walker::get_exp_window_vmf2_von_mises_intercept_fit_displacement(
                                                 const Vector<double> agent_xyz,
                                                 const unsigned agent_id) {
  std::vector<Vector<double>>& prev_xyz(agent_prev_xyz_[agent_id]);
  if (!prev_xyz.size()) {
    for (unsigned i(0); i < 2; ++i) {
      double x(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double z(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double y(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      prev_xyz.push_back(Vector<double>(x, y, z));
    }
  }
  update_new_agent_xyz(prev_xyz, agent_xyz); // now, we have at least 3 points
  Vector<double> p1(prev_xyz[0]); //earliest point
  Vector<double> p2(prev_xyz[1]);
  Vector<double> p3(prev_xyz[prev_xyz.size()-2]);
  Vector<double> p4(prev_xyz[prev_xyz.size()-1]); //latest point
  Vector<double> u(p3-p1);
  u = u.norm();
  Vector<double> v(p4-p2);
  v = v.norm();
  Vector<double> y(u.cross_product(v));
  y = y.norm();

  Vector<double> d(draw_direction_from_vmf2_intercept());
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(draw_roll_from_vector_von_mises_intercept());

  const Quaternion yQ(Quaternion(0, y).passive_rotate(roll, v));
  const Vector<double> _y(yQ.v_);
  const Quaternion wQ(Quaternion(0, v).passive_rotate(turn, _y));
  const double agent_speed = get_new_agent_speed(agent_id);
  //convert interval in s to min
  const double disp(agent_speed*get_interval()/60.0);
  const Vector<double> w(wQ.v_*disp);
  return w;
}

const Vector<double> Walker::get_synthetic_window_walk_3D_displacement(
                                                 const Vector<double> agent_xyz,
                                                 const unsigned agent_id) {
  std::vector<Vector<double>>& prev_xyz(agent_prev_xyz_[agent_id]);
  if (!prev_xyz.size()) {
    for (unsigned i(0); i < 2; ++i) {
      double x(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double z(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      double y(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
      prev_xyz.push_back(Vector<double>(x, y, z));
    }
  }
  update_new_agent_xyz(prev_xyz, agent_xyz); //now, we have at least 3 positions
  Vector<double> p1(prev_xyz[0]); //earliest position
  Vector<double> p2(prev_xyz[1]);
  Vector<double> p3(prev_xyz[prev_xyz.size()-2]);
  Vector<double> p4(prev_xyz[prev_xyz.size()-1]); //latest position
  Vector<double> u(p3-p1);
  u = u.norm();
  Vector<double> v(p4-p2);
  v = v.norm();
  Vector<double> y(u.cross_product(v));
  y = y.norm();

  const double kappa(4);
  Vector<double> mu(1, 0, 0);
  Vector<double> d(sample_vmf(kappa, mu));
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(atan2(d.y, d.z));

  const Quaternion yQ(Quaternion(0, y).passive_rotate(roll, v));
  const Vector<double> _y(yQ.v_);
  const Quaternion wQ(Quaternion(0, v).passive_rotate(turn, _y));
  const double agent_speed = get_new_agent_speed(agent_id);
  //convert interval in s to min
  const double disp(agent_speed*get_interval()/60.0);
  const Vector<double> w(wQ.v_*disp);
  return w;
}

const Vector<double> Walker::get_window_vmf2_intercept_walk_3D_displacement(
                                                 const Vector<double> agent_xyz,
                                                 const unsigned agent_id) {
  std::vector<Vector<double>>& prev_xyz(agent_prev_xyz_[agent_id]);
  if (!prev_xyz.size()) {
    prev_xyz.push_back(Vector<double>(0, 0, 0));
    prev_xyz.push_back(Vector<double>(1, 0, 0));
  }
  update_new_agent_xyz(prev_xyz, agent_xyz); // now, we have at least 3 points
  Vector<double> p1(prev_xyz[0]); //earliest point
  Vector<double> p2(prev_xyz[1]);
  Vector<double> p3(prev_xyz[prev_xyz.size()-1]); //latest point
  Vector<double> u(p2-p1);
  u = u.norm();
  Vector<double> v(p3-p2);
  v = v.norm();
  Vector<double> y(u.cross_product(v));
  y = y.norm();

  Vector<double> d(draw_direction_from_vmf2_intercept());
  const double turn(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(atan2(d.y, d.z));

  const Quaternion yQ(Quaternion(0, y).passive_rotate(roll, v));
  const Vector<double> _y(yQ.v_);
  const Quaternion wQ(Quaternion(0, v).passive_rotate(turn, _y));
  const double agent_speed = get_new_agent_speed(agent_id);
  //convert interval in s to min
  const double disp(agent_speed*get_interval()/60.0);
  const Vector<double> w(wQ.v_*disp);
  return w;
}

const Vector<double> Walker::get_vmf2_intercept_walk_3D_displacement(
                                                     const unsigned agent_id) {
  double agent_speed = get_new_agent_speed(agent_id);
  Vector<double> d(draw_direction_from_vmf2_intercept());
  const double pitch(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  const double roll(atan2(d.y, d.z));
  const Vector<double> x_axis(1, 0, 0);
  Quaternion rotateQ(Quaternion::represent_rotation(roll, x_axis)); 
  const Vector<double> y_axis(0, 1, 0);
  Quaternion pitchQ(Quaternion::represent_rotation(pitch, y_axis));
  Quaternion& orientation(agent_orientations_[agent_id]); 
  orientation = orientation.multiply(rotateQ).normalise();
  orientation = orientation.multiply(pitchQ).normalise();
  Vector<double> facing(orientation.transform(x_axis));		
  double disp(agent_speed*get_interval()/60.0); //convert interval in s to min
  const Vector<double> vector_disp(double(facing.x*disp), double(facing.y*disp),
                                  double(facing.z*disp));			
  return vector_disp;
}

const Vector<double> Walker::get_walk_3D_displacement(const unsigned agent_id) {
  double agent_speed = get_new_agent_speed(agent_id);
  double x(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
  double z(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
  double y(1 + gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
  Vector<double> facing({x, y, z});
  facing = facing.norm();
  double disp(agent_speed*get_interval()/60.0); //convert interval in s to min
  const Vector<double> vector_disp(double(facing.x*disp), double(facing.y*disp),
                                  double(facing.z*disp));			
  return vector_disp;
}

const Vector<double> Walker::get_walk_2D_displacement(const unsigned agent_id) {
  double agent_speed = get_new_agent_speed(agent_id);
  double x(gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
  double y(1 + gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
  Vector<double> facing({x, y, 0});
  facing = facing.norm();
  double disp(agent_speed*get_interval()/60.0); //convert interval in s to min
  const Vector<double> vector_disp(double(facing.x*disp), double(facing.y*disp),
                                  double(facing.z*disp));			
  return vector_disp;
}

const Vector<double> Walker::get_brownian_walk_displacement(
                                                     const unsigned agent_id) {
  double agent_speed = gsl_ran_flat(get_model().get_gsl_rng(), 0, 30);
  //const double agent_speed(draw_next_speed_vs_none_gamma(agent_id));
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  Vector<double> facing(get_random_direction());
  double disp(agent_speed*get_interval()/60.0); //convert interval in s to min
  const Vector<double> vector_disp(facing.x*disp, facing.y*disp, facing.z*disp);
  return vector_disp;
}

const Vector<double> Walker::get_brownian_uniform_turn_displacement(
                                                const Vector<double> agent_xyz,
                                                const unsigned agent_id) {
  double agent_speed = gsl_ran_flat(get_model().get_gsl_rng(), 0, 30);
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  const double turn(gsl_ran_flat(get_model().get_gsl_rng(), 0, M_PI));
  const double roll(gsl_ran_flat(get_model().get_gsl_rng(), -M_PI, M_PI));
  Vector<double> v, y;
  get_vectors_v_and_y(agent_xyz, agent_id, v, y);
  const Quaternion y_nod(Quaternion(0, y).passive_rotate(roll, v));
  const Quaternion w(Quaternion(0, v).passive_rotate(turn, y_nod.v_));
  const double disp(agent_speed*get_interval()/60.0);
  return w.v_*disp;
}

const Vector<double> Walker::get_random_directed_walk_displacement(
                                                     const unsigned agent_id) {
  double agent_speed = get_new_agent_speed(agent_id);
  Vector<double> facing(facing_[agent_id]);
  double disp(agent_speed*get_interval()/60.0); //convert interval in s to min
  const Vector<double> vector_disp(facing.x*disp, facing.y*disp, facing.z*disp);
  return vector_disp;
}

const Vector<double> Walker::get_uniform_directed_walk_displacement(
                                                     const unsigned agent_id) {
  double agent_speed = get_new_agent_speed(agent_id);
  Vector<double> facing({0, 1, 0});
  double disp(agent_speed*get_interval()/60.0); //convert interval in s to min
  const Vector<double> vector_disp(double(facing.x*disp), double(facing.y*disp),
                                  double(facing.z*disp));			
  return vector_disp;
}

const double Walker::get_agent_mean_speed(std::vector<double>& prev_speeds) {
  double mean_speed(0);
  unsigned size(prev_speeds.size());
  for (unsigned i(0); i < size; ++i) {
    mean_speed += prev_speeds[i];
  }
  return mean_speed/size;
}

const double Walker::get_new_mean_speed(const unsigned agent_id) {
  if (agent_prev_speeds_[agent_id].size() < window_size_) {
    return agent_params_[mean_speed_][agent_id];
  }
  else {
    double mean_speed(0);
    for (unsigned i(0); i < window_size_; ++i) {
      mean_speed += agent_prev_speeds_[agent_id][i];
    }
    return mean_speed/window_size_;
  }
}

const double Walker::get_new_std_mean_speed(const unsigned agent_id,
                                           const double mean_speed) {
  double alpha(0);
  double loc(0);
  double scale(0);
  unsigned N(params_[alpha_][std_at_mean_speed_]);
  const double x(mean_speed);
  for (unsigned i(1); i < N+1; ++i) {
    //p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    alpha += params_[alpha_][std_at_mean_speed_+i]*pow(x, N-i);
  }
  N = params_[sigma_][std_at_mean_speed_];
  for (unsigned i(1); i < N+1; ++i) {
    //p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    scale += params_[sigma_][std_at_mean_speed_+i]*pow(x, N-i);
  }
  /*std::cout << "alpha:" << alpha << " scale:" << scale <<
    " mean_speed:" << x << " std:" << 
    gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale) << std::endl;
    */
  double std_mean_speed;
  do {
    std_mean_speed = loc + 
      gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale);
    /*
    std::cout << "std_mean_speed:" << std_mean_speed << " loc:" << loc <<
      " alpha:" << alpha << " scale:" << scale << " mean_speed:" << mean_speed 
      << std::endl;
      */
  }
  while (std_mean_speed < 0);
  return std_mean_speed;
}

void Walker::update_new_agent_speed(std::vector<double>& prev_speeds,
                                    const double agent_speed) {
  if (prev_speeds.size() < window_size_) {
    prev_speeds.push_back(agent_speed);
  }
  else {
    for (unsigned i(1); i < window_size_; ++i) { 
      prev_speeds[i-1] = prev_speeds[i];
    }
    prev_speeds[window_size_-1] = agent_speed;
  }
}

void Walker::update_new_agent_xyz(std::vector<Vector<double>>& prev_xyz,
                                  const Vector<double> agent_xyz) {
  if (prev_xyz.size() < window_size_+2) {
    prev_xyz.push_back(agent_xyz);
  }
  else {
    for (unsigned i(1); i < prev_xyz.size(); ++i) { 
      prev_xyz[i-1] = prev_xyz[i];
    }
    prev_xyz[prev_xyz.size()-1] = agent_xyz;
  }
}

const double Walker::draw_poly_gamma_distribution_at_col(const unsigned col,
                                                        const double x,
                                                        const double max_value) {
  double alpha(0);
  double loc(0);
  double scale(0);
  unsigned N(params_[alpha_][col]);
  for (unsigned i(1); i < N+1; ++i) {
    //p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    alpha += params_[alpha_][col+i]*pow(x, N-i);
  }
  N = params_[sigma_][col];
  for (unsigned i(1); i < N+1; ++i) {
    //p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    scale += params_[sigma_][col+i]*pow(x, N-i);
  }
  double value;
  do {
    value = loc + gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale);
  }
  while (value < 0 || value > max_value);
  return value;
}

const double Walker::draw_next_speed_xcorrelated_beta(const unsigned agent_id) {
  const double max_speed(30);
  const double prev_speed(agent_prev_speeds_[agent_id].back()/max_speed);
  double prev_turn(0);
  if (agent_prev_turns_[agent_id].size()) {
    prev_turn = agent_prev_turns_[agent_id].back();
  }
  const double mean_speed(get_agent_mean_speed(agent_prev_speeds_[agent_id])/
                         max_speed);
  const double mean_turn(get_agent_mean_speed(agent_prev_turns_[agent_id]));
  const unsigned col(next_speed_);

  const double s_intercept(params_[alpha_][col+1]);
  const double s_speed(params_[alpha_][col+2]);
  const double s_turn(params_[alpha_][col+3]);
  const double s_mean_speed(params_[alpha_][col+4]);
  const double s_mean_turn(params_[alpha_][col+5]);
  const double s_mean_speed_turn(params_[alpha_][col+6]);
  const double s_speed_squared(params_[alpha_][col+7]);
  const double shape(exp(s_intercept +
                        s_speed*prev_speed +
                        s_turn*prev_turn +
                        s_mean_speed*mean_speed +
                        s_mean_turn*mean_turn +
                        s_mean_speed_turn*(mean_speed/mean_turn) +
                        s_speed_squared*(prev_speed*prev_speed)));

  const double r_intercept(params_[beta_][col+1]);
  const double r_speed(params_[beta_][col+2]);
  const double r_turn(params_[beta_][col+3]);
  const double r_mean_speed(params_[beta_][col+4]);
  const double r_mean_turn(params_[beta_][col+5]);
  const double r_mean_speed_turn(params_[beta_][col+6]);
  const double r_speed_squared(params_[beta_][col+7]);
  const double rate(exp(r_intercept +
                       r_speed*prev_speed +
                       r_turn*prev_turn +
                       r_mean_speed*mean_speed +
                       r_mean_turn*mean_turn +
                       r_mean_speed_turn*(mean_speed/mean_turn) +
                       r_speed_squared*(prev_speed*prev_speed)));

  return gsl_ran_beta(get_model().get_gsl_rng(), rate, shape)*max_speed;
}

const double Walker::draw_next_speed_vs_s_t_r_s2_beta(const unsigned agent_id) {
  const double max_speed(30);
  const double prev_speed(agent_prev_speeds_[agent_id].back()/max_speed);
  double prev_turn(0);
  if (agent_prev_turns_[agent_id].size()) {
    prev_turn = agent_prev_turns_[agent_id].back();
  }
  double prev_roll(0);
  if (agent_prev_rolls_[agent_id].size()) {
    prev_roll = agent_prev_rolls_[agent_id].back();
  }
  const unsigned col(next_speed_);

  const double s_intercept(params_[alpha_][col+1]);
  const double s_speed(params_[alpha_][col+2]);
  const double s_turn(params_[alpha_][col+3]);
  const double s_roll(params_[alpha_][col+4]);
  const double s_speed_squared(params_[alpha_][col+5]);
  const double shape(exp(s_intercept +
                        s_speed*prev_speed +
                        s_turn*prev_turn +
                        s_roll*prev_roll +
                        s_speed_squared*(prev_speed*prev_speed)));

  const double r_intercept(params_[beta_][col+1]);
  const double r_speed(params_[beta_][col+2]);
  const double r_turn(params_[beta_][col+3]);
  const double r_roll(params_[beta_][col+4]);
  const double r_speed_squared(params_[beta_][col+5]);
  const double rate(exp(r_intercept +
                       r_speed*prev_speed +
                       r_turn*prev_turn +
                       r_roll*prev_roll +
                       r_speed_squared*(prev_speed*prev_speed)));

  return gsl_ran_beta(get_model().get_gsl_rng(), rate, shape)*max_speed;
}

const double Walker::draw_next_speed_vs_s_beta(const unsigned agent_id) {
  const double max_speed(30);
  const double prev_speed(agent_prev_speeds_[agent_id].back()/max_speed);
  const unsigned col(next_speed_);

  const double s_intercept(params_[alpha_][col+1]);
  const double s_speed(params_[alpha_][col+2]);
  const double s_speed_squared(params_[alpha_][col+3]);
  const double shape(exp(s_intercept +
                        s_speed*prev_speed +
                        s_speed_squared*(prev_speed*prev_speed)));

  const double r_intercept(params_[beta_][col+1]);
  const double r_speed(params_[beta_][col+2]);
  const double r_speed_squared(params_[beta_][col+3]);
  const double rate(exp(r_intercept +
                       r_speed*prev_speed +
                       r_speed_squared*(prev_speed*prev_speed)));

  return gsl_ran_beta(get_model().get_gsl_rng(), rate, shape)*max_speed;
}



const double Walker::draw_next_speed_vs_s_t_r_s2_gamma(
                                                   const unsigned agent_id) {
  const double max_speed(30);
  const double prev_speed(agent_prev_speeds_[agent_id].back()/max_speed);
  double prev_turn(0);
  if (agent_prev_turns_[agent_id].size()) {
    prev_turn = agent_prev_turns_[agent_id].back();
  }
  //double prev_roll(0);
  //if (agent_prev_rolls_[agent_id].size()) {
  //  prev_roll = fabs(agent_prev_rolls_[agent_id].back());
  //}
  const unsigned col(next_speed_);

  const double s_intercept(params_[alpha_][col+1]);
  const double s_speed(params_[alpha_][col+2]);
  const double s_turn(params_[alpha_][col+3]);
  //const double s_roll(params_[alpha_][col+4]);
  const double s_speed_squared(params_[alpha_][col+4]);
  const double shape(exp(s_intercept +
                        s_speed*prev_speed +
                        s_turn*prev_turn +
                        //s_roll*prev_roll +
                        s_speed_squared*(prev_speed*prev_speed)));

  const double r_intercept(params_[beta_][col+1]);
  const double r_speed(params_[beta_][col+2]);
  const double r_turn(params_[beta_][col+3]);
  //const double r_roll(params_[beta_][col+4]);
  const double r_speed_squared(params_[beta_][col+4]);
  const double rate(exp(r_intercept +
                       r_speed*prev_speed +
                       r_turn*prev_turn +
                       //r_roll*prev_roll +
                       r_speed_squared*(prev_speed*prev_speed)));
  double speed;
  do {
    speed = gsl_ran_gamma(get_model().get_gsl_rng(), shape, 1/rate)*max_speed;
    //std::cout << "speed:" << speed << std::endl;
  }
  while (speed < 0 || speed > max_speed);
  return speed;
}

const double Walker::draw_next_speed_vs_sp_gamma(const unsigned agent_id,
                                                 double prev_speed) {
  const double max_speed(30);
  prev_speed /= max_speed;
  const unsigned col(next_speed_);
  const double s_intercept(params_[alpha_][col+1]);
  const double s_speed(params_[alpha_][col+2]);
  const double s_speed_squared(params_[alpha_][col+3]);
  const double shape(exp(s_intercept +
                        s_speed*prev_speed +
                        s_speed_squared*(prev_speed*prev_speed)));

  const double r_intercept(params_[beta_][col+1]);
  const double r_speed(params_[beta_][col+2]);
  const double r_speed_squared(params_[beta_][col+3]);
  const double rate(exp(r_intercept +
                       r_speed*prev_speed +
                       r_speed_squared*(prev_speed*prev_speed)));
  double speed;
  do {
    speed = gsl_ran_gamma(get_model().get_gsl_rng(), shape, 1/rate)*max_speed;
    //std::cout << "speed:" << speed << std::endl;
  }
  while (speed <= 0 || speed >= max_speed);
  return speed;
}


const double Walker::draw_next_speed_vs_none_gamma(const unsigned agent_id) {
  const double max_speed(30);
  const unsigned col(next_speed_);
  const double s_intercept(params_[alpha_][col+1]);
  const double shape(s_intercept);
  const double r_intercept(params_[beta_][col+1]);
  const double rate(r_intercept);
  double speed;
  do {
    speed = gsl_ran_gamma(get_model().get_gsl_rng(), shape, 1/rate)*max_speed;
  }
  while (speed < 0 || speed > max_speed);
  return speed;
}

const double Walker::draw_next_speed_xcorrelated_gamma(const unsigned agent_id) {
  const double max_speed(30);
  const double prev_speed(agent_prev_speeds_[agent_id].back()/max_speed);
  double prev_turn(0);
  if (agent_prev_turns_[agent_id].size()) {
    prev_turn = agent_prev_turns_[agent_id].back();
  }
  const double mean_speed(get_agent_mean_speed(agent_prev_speeds_[agent_id])/
                         max_speed);
  const double mean_turn(get_agent_mean_speed(agent_prev_turns_[agent_id]));
  const unsigned col(next_speed_);

  const double s_intercept(params_[alpha_][col+1]);
  const double s_speed(params_[alpha_][col+2]);
  const double s_turn(params_[alpha_][col+3]);
  const double s_mean_speed(params_[alpha_][col+4]);
  const double s_mean_turn(params_[alpha_][col+5]);
  const double s_mean_speed_turn(params_[alpha_][col+6]);
  const double s_speed_squared(params_[alpha_][col+7]);
  const double shape(exp(s_intercept +
                        s_speed*prev_speed +
                        s_turn*prev_turn +
                        s_mean_speed*mean_speed +
                        s_mean_turn*mean_turn +
                        s_mean_speed_turn*(mean_speed/mean_turn) +
                        s_speed_squared*(prev_speed*prev_speed)));

  const double r_intercept(params_[beta_][col+1]);
  const double r_speed(params_[beta_][col+2]);
  const double r_turn(params_[beta_][col+3]);
  const double r_mean_speed(params_[beta_][col+4]);
  const double r_mean_turn(params_[beta_][col+5]);
  const double r_mean_speed_turn(params_[beta_][col+6]);
  const double r_speed_squared(params_[beta_][col+7]);
  const double rate(exp(r_intercept +
                       r_speed*prev_speed +
                       r_turn*prev_turn +
                       r_mean_speed*mean_speed +
                       r_mean_turn*mean_turn +
                       r_mean_speed_turn*(mean_speed/mean_turn) +
                       r_speed_squared*(prev_speed*prev_speed)));

  double speed;
  do {
    speed = gsl_ran_gamma(get_model().get_gsl_rng(), shape, 1/rate)*max_speed;
    //std::cout << "speed:" << speed << std::endl;
  }
  while (speed < 0 || speed > max_speed);
  return speed;
}

const Vector<double> Walker::sample_vmf(const double kappa,
                                       const Vector<double> mu) {
  const unsigned dim(3);
  //sample offset from center (on sphere) with spread kappa
  const double w(sample_weight(kappa, dim));
  //sample a point v on the unit sphere that's orthogonal to mu
  const Vector<double> v(sample_orthonormal_to(mu));
  //compute new point
  return v*sqrt(1.-w*w) + mu*w;
}

const double Walker::sample_von_mises(const double mu, const double kappa) {
  double s;
  double U, V, W, Y, Z;
  double result, mod;
  int neg;
  if (kappa < 1e-8) {
    return M_PI * (2 * gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) - 1);
  } else {
    /* with double precision rho is zero until 1.4e-8 */
    if (kappa < 1e-5) {
      /*
       * second order taylor expansion around kappa = 0
       * precise until relatively large kappas as second order is 0
       */
      s = (1. / kappa + kappa);
    } else {
      double r = 1 + sqrt(1 + 4 * kappa * kappa);
      double rho = (r - sqrt(2 * r)) / (2 * kappa);
      s = (1 + rho * rho) / (2 * rho);
    }

    while (1) {
      U = gsl_ran_flat(get_model().get_gsl_rng(), 0, 1);
      Z = cos(M_PI * U);
      W = (1 + s * Z) / (s + Z);
      Y = kappa * (s - W);
      V = gsl_ran_flat(get_model().get_gsl_rng(), 0, 1);
      /*
       * V==0.0 is ok here since Y >= 0 always leads
       * to accept, while Y < 0 always rejects
       */
      if ((Y * (2 - Y) - V >= 0) || (log(Y / V) + 1 - Y >= 0)) {
        break;
      }
    }

    U = gsl_ran_flat(get_model().get_gsl_rng(), 0, 1);

    result = acos(W);
    if (U < 0.5) {
      result = -result;
    }
    result += mu;
    neg = (result < 0);
    mod = fabs(result);
    mod = (fmod(mod + M_PI, 2 * M_PI) - M_PI);
    if (neg) {
      mod *= -1;
    }
    return mod;
  }
}



const double Walker::sample_weight(const double kappa, unsigned dim) {
  //Rejection sampling scheme for sampling distance from center on
  //surface of the sphere.
  dim -= 1;
  double b(dim / (sqrt(4. * kappa*kappa + dim*dim) + 2 * kappa));
  double x((1. - b) / (1. + b));
  double c(kappa * x + dim * log(1 - x*x)); 
  while (1) {
    double z(gsl_ran_beta(get_model().get_gsl_rng(), dim/2., dim / 2.));
    double w((1. - (1. + b) * z) / (1. - (1. - b) * z));
    double u(gsl_ran_flat(get_model().get_gsl_rng(), 0, 1));
    if (kappa*w + dim*log(1. - x*w) - c >= log(u)) {
      return w;
    }
  }
}

const Vector<double> Walker::sample_orthonormal_to(const Vector<double> mu) {
  //Sample point on sphere orthogonal to mu
  Vector<double> v(gsl_ran_gaussian(get_model().get_gsl_rng(), 1),
                  gsl_ran_gaussian(get_model().get_gsl_rng(), 1),
                  gsl_ran_gaussian(get_model().get_gsl_rng(), 1));
  Vector<double> proj_mu_v(mu*mu.dot_product(v) / mu.magnitude());
  Vector<double> orthto(v - proj_mu_v);
  return orthto.norm();
}

const Vector<double> Walker::draw_poly_vmf_distribution_at_col(
                                                      const unsigned col,
                                                      const double x,
                                                      const double prev_x) {
  const double kappa_uniform(0.0001);
  double kappa(0);
  double weight(0);
  Vector<double> mu(1, 0, 0);
  unsigned N(params_[alpha_][col]);
  for (unsigned i(1); i < N+1; ++i) {
    //p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    kappa += params_[alpha_][col+i]*pow(x, N-i);
  }
  N = params_[beta_][col];
  for (unsigned i(1); i < N+1; ++i) {
    //p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    weight += params_[beta_][col+i]*pow(x, N-i);
  }
  mu = Vector<double>(params_[loc_][col+1], params_[loc_][col+2],
                     params_[loc_][col+3]);

  //kappa = std::max(kappa, double(4));
  /*
  if (prev_x*180/M_PI < 50) {
    //kappa = std::max(kappa, double(4));
    weight = std::max(weight, double(1.0));
  }
  */
  /*
  if (x < 10) {
    weight = std::max(weight, double(0.3));
    kappa = std::max(kappa, double(6));
  }
  */
  //weight = std::max(weight, double(0.9));
  //if (prev_x*180/M_PI < 50 || 
  //    gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= weight) { 
  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= weight) { 
    return sample_vmf(kappa, mu);
  }
  return sample_vmf(kappa_uniform, mu);
}

const double Walker::draw_roll_from_vector_von_mises_intercept() {
  const unsigned col(next_roll_);
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa1(std::min(1e+7, exp(kappa_intercept1)));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa2(std::min(1e+7, exp(kappa_intercept2)));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta(std::max(0.0, std::min(1.0, (theta_intercept + 2.0)/4.0)));

  double mu1(params_[sigma_][col+1]);
  double mu2(params_[max_][col+1]);

  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= theta) { 
    return sample_von_mises(mu1, kappa1);
  }
  return sample_von_mises(mu2, kappa2);
}

const double Walker::draw_roll_from_vector_von_mises(const double agent_speed,
                                                    const double agent_turn,
                                                    const double agent_roll) {
  const double speed(agent_speed/17);
  double zz = fabs(sin(agent_turn)*cos(agent_roll));
  double yy = fabs(sin(agent_turn)*sin(agent_roll));
  double xx = cos(agent_turn);

  const unsigned col(next_roll_);
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa_speed1(params_[alpha_][col+2]);
  const double kappa_xx1(params_[alpha_][col+3]);
  const double kappa_yy1(params_[alpha_][col+4]);
  const double kappa_zz1(params_[alpha_][col+5]);
  const double kappa1(exp(kappa_intercept1 + kappa_speed1*speed + kappa_xx1*xx +
                         kappa_yy1*yy + kappa_zz1*zz));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa_speed2(params_[beta_][col+2]);
  const double kappa_xx2(params_[beta_][col+3]);
  const double kappa_yy2(params_[beta_][col+4]);
  const double kappa_zz2(params_[beta_][col+5]);
  const double kappa2(exp(kappa_intercept2 + kappa_speed2*speed + kappa_xx2*xx +
                         kappa_yy2*yy + kappa_zz2*zz));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta_speed(params_[loc_][col+2]);
  const double theta_xx(params_[loc_][col+3]);
  const double theta_yy(params_[loc_][col+4]);
  const double theta_zz(params_[loc_][col+5]);
  const double weight(1/(1+exp(-theta_intercept - theta_speed*speed - 
                              theta_xx*xx - theta_yy*yy - theta_zz*zz)));

  double mu1(params_[sigma_][col+1]);
  double mu2(params_[max_][col+1]);

  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= weight) { 
    double angle(sample_von_mises(mu1, kappa1));
    std::cout << "0:" << " " << angle*180/M_PI << " kappa:" << kappa1 << 
      " mu:" << mu1*180/M_PI << std::endl;
    return angle;
  }
  double angle(sample_von_mises(mu2, kappa2));
  std::cout << "180:" << " " << angle*180/M_PI << " kappa:" << kappa2 << 
    " mu:" << mu2*180/M_PI << std::endl;
  return angle;
}

const double Walker::draw_roll_from_von_mises(const double agent_speed,
                                             const double turn) {
  const unsigned col(next_roll_);
  const double speed(agent_speed/17); // 0 < prev_speed < 30
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa_speed1(params_[alpha_][col+2]);
  const double kappa_turn1(params_[alpha_][col+3]);
  const double kappa1(exp(kappa_intercept1 + kappa_speed1*speed +
                         kappa_turn1*turn));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa_speed2(params_[beta_][col+2]);
  const double kappa_turn2(params_[beta_][col+3]);
  const double kappa2(exp(kappa_intercept2 + kappa_speed2*speed +
                         kappa_turn2*turn));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta_speed(params_[loc_][col+2]);
  const double theta_turn(params_[loc_][col+3]);
  const double weight(1/(1+exp(-theta_intercept - theta_speed*speed - 
                              theta_turn*turn)));

  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= weight) { 
    double angle(sample_von_mises(0, kappa1));
    //std::cout << "0:" << " " << angle*180/M_PI << std::endl;
    return angle;
  }
  double angle(sample_von_mises(M_PI, kappa2));
  //std::cout << "180:" << " " << angle*180/M_PI << " " << std::endl;
  return angle;
}

const double Walker::draw_next_roll_vs_none_von_mises() {
  const unsigned col(next_roll_);
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa1(std::min(1e+7, exp(kappa_intercept1)));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa2(std::min(1e+7, exp(kappa_intercept2)));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta(1/(1+exp(-theta_intercept))); 

  double mu1(params_[sigma_][col+1]);
  double mu2(params_[max_][col+1]);

  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= theta) { 
    return sample_von_mises(mu1, kappa1);
  }
  return sample_von_mises(mu2, kappa2);
}

const double Walker::draw_next_roll_vs_t_von_mises(const double turn) {
  const unsigned col(next_roll_);
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa_turn1(params_[alpha_][col+2]);
  const double kappa1(std::min(1e+7, exp(kappa_intercept1 +
                                         kappa_turn1*turn)));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa_turn2(params_[beta_][col+2]);
  const double kappa2(std::min(1e+7, exp(kappa_intercept2 +
                                         kappa_turn2*turn)));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta_turn(params_[loc_][col+2]);
  const double theta(1/(1+exp(-theta_intercept -
                              theta_turn*turn))); 

  double mu1(params_[sigma_][col+1]);
  double mu2(params_[max_][col+1]);

  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= theta) { 
    return sample_von_mises(mu1, kappa1);
  }
  return sample_von_mises(mu2, kappa2);
}

const double Walker::draw_next_roll_vs_s_von_mises(const double agent_speed) {
  const double speed(agent_speed/17); // 0 < prev_speed < 30
  const unsigned col(next_roll_);
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa_speed1(params_[alpha_][col+2]);
  const double kappa1(std::min(1e+7, exp(kappa_intercept1 +
                                         kappa_speed1*speed)));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa_speed2(params_[beta_][col+2]);
  const double kappa2(std::min(1e+7, exp(kappa_intercept2 +
                                         kappa_speed2*speed)));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta_speed(params_[loc_][col+2]);
  const double theta(1/(1+exp(-theta_intercept -
                              theta_speed*speed))); 

  double mu1(params_[sigma_][col+1]);
  double mu2(params_[max_][col+1]);

  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= theta) { 
    return sample_von_mises(mu1, kappa1);
  }
  return sample_von_mises(mu2, kappa2);
}

const Vector<double> Walker::draw_next_direction_vs_s_tp_vmf(
                                            const double agent_speed,
                                            const double prev_turn) {
  const double speed(agent_speed/17); // 0 < prev_speed < 30
  const unsigned col(next_direction_);
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa_speed1(params_[alpha_][col+2]);
  const double kappa_turn1(params_[alpha_][col+3]);
  const double kappa1(std::min(1e+7, exp(kappa_intercept1 +
                                         kappa_speed1*speed +
                                         kappa_turn1*prev_turn)));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa_speed2(params_[beta_][col+2]);
  const double kappa_turn2(params_[beta_][col+3]);
  const double kappa2(std::min(1e+7, exp(kappa_intercept2 +
                                         kappa_speed2*speed +
                                         kappa_turn2*prev_turn)));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta_speed(params_[loc_][col+2]);
  const double theta_turn(params_[loc_][col+3]);
  const double theta(1/(1+exp(-theta_intercept -
                              theta_speed*speed -
                              theta_turn*prev_turn))); 

  const Vector<double> mu1(Vector<double>(params_[sigma_][col+1],
                                        params_[sigma_][col+2],
                                        params_[sigma_][col+3]));
  const Vector<double> mu2(Vector<double>(params_[max_][col+1],
                                        params_[max_][col+2],
                                        params_[max_][col+3]));

  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= theta) { 
    return sample_vmf(kappa1, mu1);
  }
  return sample_vmf(kappa2, mu2);
}

const Vector<double> Walker::draw_next_direction_vs_s_vmf(
                                            const double agent_speed) {
  const double speed(agent_speed/17); // 0 < prev_speed < 30
  const unsigned col(next_direction_);
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa_speed1(params_[alpha_][col+2]);
  const double kappa1(std::min(1e+7, exp(kappa_intercept1 +
                                         kappa_speed1*speed)));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa_speed2(params_[beta_][col+2]);
  const double kappa2(std::min(1e+7, exp(kappa_intercept2 +
                                         kappa_speed2*speed)));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta_speed(params_[loc_][col+2]);
  const double theta(1/(1+exp(-theta_intercept -
                              theta_speed*speed))); 

  const Vector<double> mu1(Vector<double>(params_[sigma_][col+1],
                                        params_[sigma_][col+2],
                                        params_[sigma_][col+3]));
  const Vector<double> mu2(Vector<double>(params_[max_][col+1],
                                        params_[max_][col+2],
                                        params_[max_][col+3]));

  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= theta) { 
    return sample_vmf(kappa1, mu1);
  }
  return sample_vmf(kappa2, mu2);
}

const Vector<double> Walker::draw_next_direction_vs_none_vmf() {
  const unsigned col(next_direction_);
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa1(std::min(1e+7, exp(kappa_intercept1)));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa2(std::min(1e+7, exp(kappa_intercept2)));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta(1/(1+exp(-theta_intercept))); 

  const Vector<double> mu1(Vector<double>(params_[sigma_][col+1],
                                        params_[sigma_][col+2],
                                        params_[sigma_][col+3]));
  const Vector<double> mu2(Vector<double>(params_[max_][col+1],
                                        params_[max_][col+2],
                                        params_[max_][col+3]));

  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= theta) { 
    return sample_vmf(kappa1, mu1);
  }
  return sample_vmf(kappa2, mu2);
}

const Vector<double> Walker::draw_direction_from_vmf2_intercept() {
  const unsigned col(next_direction_);
  const double kappa_intercept1(params_[alpha_][col+1]);
  const double kappa1(std::min(1e+7, exp(kappa_intercept1)));

  const double kappa_intercept2(params_[beta_][col+1]);
  const double kappa2(std::min(1e+7, exp(kappa_intercept2)));

  const double theta_intercept(params_[loc_][col+1]);
  const double theta(std::max(0.0, std::min(1.0, (theta_intercept + 2.0)/4.0)));

  const Vector<double> mu1(Vector<double>(params_[sigma_][col+1],
                                        params_[sigma_][col+2],
                                        params_[sigma_][col+3]));

  const Vector<double> mu2(Vector<double>(params_[max_][col+1],
                                        params_[max_][col+2],
                                        params_[max_][col+3]));
  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= theta) { 
    return sample_vmf(kappa1, mu1);
  }
  return sample_vmf(kappa2, mu2);
}

const Vector<double> Walker::draw_direction_from_vector_vmf(
                                                     const double agent_speed,
                                                     const double agent_turn,
                                                     const double agent_roll) {

  const double speed(agent_speed/17);
  double zz = sin(agent_turn)*cos(agent_roll)+1;
  double yy = sin(agent_turn)*sin(agent_roll)+1;
  double xx = cos(agent_turn)+1;

  const unsigned col(next_direction_);
  const double kappa_intercept(params_[alpha_][col+1]);
  const double kappa_speed(params_[alpha_][col+2]);
  const double kappa_xx(params_[alpha_][col+3]);
  const double kappa_yy(params_[alpha_][col+4]);
  const double kappa_zz(params_[alpha_][col+5]);
  const double kappa(exp(kappa_intercept + kappa_speed*speed + kappa_xx*xx +
                        kappa_yy*yy + kappa_zz*zz));

  const double theta_intercept(params_[beta_][col+1]);
  const double theta_speed(params_[beta_][col+2]);
  const double theta_xx(params_[beta_][col+3]);
  const double theta_yy(params_[beta_][col+4]);
  const double theta_zz(params_[beta_][col+5]);
  const double weight(1/(1+exp(-theta_intercept - theta_speed*speed - 
                              theta_xx*xx - theta_yy*yy - theta_zz*zz)));

  const Vector<double> mu1(Vector<double>(params_[loc_][col+1],
                                        params_[loc_][col+2],
                                        params_[loc_][col+3]));

  const Vector<double> mu2(Vector<double>(params_[sigma_][col+1],
                                        params_[sigma_][col+2],
                                        params_[sigma_][col+3]));

  /*
  std::cout << "weight:" << weight << " kappa:" << kappa << " speed:" << 
    speed*17 << " turn:" << turn*100 << " roll:" << roll*180 << std::endl;
    */
  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= weight) { 
    return sample_vmf(kappa, mu1);
  }
  const double kappa_uniform(0.0001);
  return sample_vmf(kappa_uniform, mu2);
}

const Vector<double> Walker::draw_exp_poly_vmf_distribution_at_col(
                                                      const unsigned col,
                                                      const double prev_speed,
                                                      const double prev_turn,
                                                      const double prev_roll) {
  const double speed(prev_speed/17); // 0 < prev_speed < 30
  const double turn(prev_turn*180/M_PI/100); // 0 < prev_turn < PI
  const double roll(fabs(prev_roll)/M_PI); //-PI < prev_roll < PI
  const double kappa_intercept(params_[alpha_][col+1]);
  const double kappa_speed(params_[alpha_][col+2]);
  const double kappa_turn(params_[alpha_][col+3]);
  const double kappa_roll(params_[alpha_][col+4]);
  const double kappa(exp(kappa_intercept + kappa_speed*speed + kappa_turn*turn +
                        kappa_roll*roll));

  const double theta_intercept(params_[beta_][col+1]);
  const double theta_speed(params_[beta_][col+2]);
  const double theta_turn(params_[beta_][col+3]);
  const double theta_roll(params_[beta_][col+4]);
  const double weight(1/(1+exp(-theta_intercept - theta_speed*speed - 
                              theta_turn*turn - theta_roll*roll)));

  const Vector<double> mu(Vector<double>(params_[loc_][col+1],
                                       params_[loc_][col+2],
                                       params_[loc_][col+3]));

  /*
  std::cout << "weight:" << weight << " kappa:" << kappa << " speed:" << 
    speed*17 << " turn:" << turn*100 << " roll:" << roll*180 << std::endl;
    */
  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= weight) { 
    return sample_vmf(kappa, mu);
  }
  const double kappa_uniform(0.0001);
  return sample_vmf(kappa_uniform, mu);
}



/*
// uncorrelated agent speed
const double Walker::get_new_agent_speed(const unsigned agent_id) {
  const double alpha(2.16);
  const double loc(0.00);
  const double scale(2.78);
  speed_max_ = (28);
  double agent_speed(loc + 
                    gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale));
  while (agent_speed > speed_max_ || agent_speed < 0) {
    agent_speed = loc + gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale);
  }
  return agent_speed;
}
*/

/*
// autocorrelated agent speed
const double Walker::get_new_agent_speed(const unsigned agent_id) {
  double agent_speed(0);
  if (!agent_prev_speeds_[agent_id].size()) {
    agent_speed = agent_params_[instant_speed_][agent_id];
    agent_prev_speeds_[agent_id].push_back(agent_speed);
  }
  else {
    agent_speed = draw_poly_gamma_distribution_at_col(next_speed_,
                                              agent_prev_speeds_[agent_id][0],
                                              speed_max_);
    agent_prev_speeds_[agent_id][0] = agent_speed;
  }
  return agent_speed;
}
*/

//cross-correlated agent speed
const double Walker::get_new_agent_speed(const unsigned agent_id) {
  double agent_speed(0);
  if (!agent_prev_speeds_[agent_id].size()) {
    agent_speed = agent_params_[instant_speed_][agent_id];
  }
  else {
    agent_speed = draw_next_speed_vs_s_t_r_s2_gamma(agent_id);
    //agent_speed = agent_params_[instant_speed_][agent_id];
    //agent_speed = draw_next_speed_xcorrelated_gamma(agent_id);
  }
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  return agent_speed;
}


/*
// uncorrelated agent turn with speed 
void Walker::update_correlated_orientation(const unsigned agent_id,
                                           const double agent_speed) {

  const double alpha(1.94);
  const double loc(0.00);
  const double scale(30.39);
  double turn_max(180);
  double turn(loc + gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale));
  while (turn > turn_max || turn < 0) {
    turn = loc + gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale);
  }

  const double roll(gsl_ran_flat(get_model().get_gsl_rng(), 0, 2.0*M_PI));
  const double pitch(turn*M_PI/180);

  const Vector<double> x_axis(1, 0, 0);
  Quaternion rotateQ(Quaternion::represent_rotation(roll, x_axis.x, x_axis.y,
                                                   x_axis.z)); 
  const Vector<double> y_axis(0, 1, 0);
  Quaternion pitchQ(Quaternion::represent_rotation(pitch, y_axis.x, y_axis.y,
                                                  y_axis.z));
  Quaternion& orientation(agent_orientations_[agent_id]); 
  orientation = orientation.multiply(rotateQ).normalise();
  orientation = orientation.multiply(pitchQ).normalise();
}
*/


// correlated agent turn with previous speed, turn and roll
void Walker::update_correlated_orientation(const unsigned agent_id,
                                           const double agent_speed) {
  double prev_speed(agent_speed);
  double prev_turn(1/agent_speed*M_PI);
  double prev_roll(gsl_ran_flat(get_model().get_gsl_rng(), -M_PI, M_PI));
  if (agent_prev_speeds_[agent_id].size() > 1) { 
    prev_speed = agent_prev_speeds_[agent_id][
      agent_prev_speeds_[agent_id].size()-2];

  }
  if (agent_prev_turns_[agent_id].size()) {
    prev_turn = agent_prev_turns_[agent_id].back();
  }
  if (agent_prev_rolls_[agent_id].size()) {
    prev_roll = agent_prev_rolls_[agent_id].back();
  }

  Vector<double> d(draw_direction_from_vector_vmf(prev_speed, prev_turn,
                                                 prev_roll));
  //Vector<double> d(draw_exp_poly_vmf_distribution_at_col(next_direction_,
  //                                        agent_prev_speeds_[agent_id].back(),
  //                                        prev_turn, prev_roll));
  //const double roll(M_PI);
  const double pitch(atan2(sqrt(d.z*d.z + d.y*d.y), d.x));
  //const double pitch(0);
  //const double roll(atan2(d.y, d.z));
  //const double roll(draw_roll_from_von_mises(agent_speed, pitch));
  const double roll(draw_roll_from_vector_von_mises(prev_speed, prev_turn,
                                                    prev_roll));
  update_new_agent_speed(agent_prev_turns_[agent_id], pitch);
  update_new_agent_speed(agent_prev_rolls_[agent_id], roll);

  const Vector<double> x_axis(1, 0, 0);
  Quaternion rotateQ(Quaternion::represent_rotation(roll, x_axis)); 
  const Vector<double> y_axis(0, 1, 0);
  Quaternion pitchQ(Quaternion::represent_rotation(pitch, y_axis));
  Quaternion& orientation(agent_orientations_[agent_id]); 
  orientation = orientation.multiply(rotateQ).normalise();
  orientation = orientation.multiply(pitchQ).normalise();
}

const Quaternion Walker::get_updated_agent_orientation(const unsigned agent_id,
                                                     const double agent_speed) {
  Quaternion& orientation(agent_orientations_[agent_id]);
  if (agent_redirection_magnitudes_[agent_id]) {
    orientation = Quaternion::face_vector(Vector<double>(
                                         agent_redirections_[agent_id].x,
                                         -agent_redirections_[agent_id].y,
                                         -agent_redirections_[agent_id].z));
  }
  else {
    switch (walk_type_) {
    case correlated_walk_:
      update_correlated_orientation(agent_id, agent_speed); 
      break;
    case uncorrelated_walk_:
      update_uncorrelated_orientation(agent_id);
      break;
    case brownian_walk_:
      update_brownian_orientation(agent_id);
      break;
    case corkscrew_:
      update_corkscrew_orientation(agent_id);
      break;
    }
  }
  return orientation;
}

void Walker::update_corkscrew_orientation(const unsigned agent_id) {
  const double pitch(5*M_PI/180);
  const double roll(-5*M_PI/180);
  const Vector<double> x_axis(1, 0, 0);
  Quaternion rotateQ(Quaternion::represent_rotation(roll, x_axis)); 
  const Vector<double> y_axis(0, 1, 0);
  Quaternion pitchQ(Quaternion::represent_rotation(pitch, y_axis));
  Quaternion& orientation(agent_orientations_[agent_id]); 
  orientation = orientation.multiply(rotateQ).normalise();
  orientation = orientation.multiply(pitchQ).normalise();
}

void Walker::update_brownian_orientation(const unsigned agent_id) {
  Quaternion& orientation(agent_orientations_[agent_id]); 
  Vector<double> agent_direction(get_random_direction());
  orientation = Quaternion::face_vector(agent_direction);
}

void Walker::update_uncorrelated_orientation(const unsigned agent_id) {
  Quaternion& orientation(agent_orientations_[agent_id]); 
  // first roll the agent:
  const double mean_roll_rate(agent_params_[mean_roll_rate_][agent_id]*
                             M_PI/180);
  const double std_roll_rate(agent_params_[std_mean_roll_rate_][agent_id]*
                            M_PI/180);
  double roll;
  if (mean_roll_rate < 0) {
    // if mean roll rate is negative, assume this indicates a uniform
    // distribution should be used. 
    roll = gsl_ran_flat(get_model().get_gsl_rng(), 0, 2.0*M_PI);
  } else {
    roll = gsl_ran_gaussian(get_model().get_gsl_rng(), std_roll_rate) +
      mean_roll_rate;
    // randomly invert roll direction. Avoids corkscrewing. 
    if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= 0.5) {
      roll *= -1.0;		// cells can roll in either direction.
    }
    roll *= get_interval()/60.0;
  }
  // roll as a quaternion.
  const Vector<double> x_axis(1, 0, 0);
  Quaternion rotateQ(Quaternion::represent_rotation(roll, x_axis)); 
  // multiply orientation by rotateQ, because rotateQ is calculated relative
  // to cell, not in absolute space. 
  // alter the cell's orientation. 
  orientation = orientation.multiply(rotateQ).normalise();

  // then pitch:
  const double mean_pitch_rate(agent_params_[mean_pitch_rate_][agent_id]*
                              M_PI/180);
  const double std_pitch_rate(agent_params_[std_mean_pitch_rate_][agent_id]*
                             M_PI/180);
  double pitch = gsl_ran_gaussian(get_model().get_gsl_rng(), std_pitch_rate) +
    mean_pitch_rate;
  // change cell pitch (roll along the y axis). Pitch can be changed in both
  // positive and negative directions.
  pitch *= get_interval()/60.0;		// account for timestep.
  const Vector<double> y_axis(0, 1, 0);
  Quaternion pitchQ(Quaternion::represent_rotation(pitch, y_axis));
  // multiply orientation by rotateQ, because pitchQ is calculated relative
  // to cell, not in absolute space.
  orientation = orientation.multiply(pitchQ).normalise();
}

Vector<double> Walker::get_random_direction() {
  std::uniform_real_distribution<> uni_z(-1, 1);
  std::uniform_real_distribution<> uni_t(0, 2*M_PI);
  double z(uni_z(get_model().get_rng()));
  double t(uni_t(get_model().get_rng()));
  Vector<double> position(sqrt(1-pow(z,2))*cos(t),
                          sqrt(1-pow(z,2))*sin(t),
                          z);
  return position.norm();
}

std::vector<Vector<double>>& Walker::get_agent_redirections() {
  return agent_redirections_;
}

std::vector<unsigned>& Walker::get_agent_redirection_magnitudes() {
  return agent_redirection_magnitudes_;
}

double Walker::step() {
  motile_species_.walk();
  return Process::step();
}
/*
double Walker::draw_gamma_distribution_at_col(const int col,
                                             const unsigned agent_id) {
  double alpha(0);
  double loc(0);
  double scale(0);
  if (col == std_at_mean_speed_) {
    const double x(agent_params_[mean_speed_][agent_id]);
    double x_mean(params_[alpha_][std_at_mean_speed_]);
    double x_std(params_[alpha_][std_at_mean_speed_+1]);
    double y_mean(params_[alpha_][std_at_mean_speed_+2]);
    double y_std(params_[alpha_][std_at_mean_speed_+3]);
    double a(params_[alpha_][std_at_mean_speed_+4]);
    double t(params_[alpha_][std_at_mean_speed_+5]);
    double c(params_[alpha_][std_at_mean_speed_+6]);
    double x_norm((x - x_mean)/x_std);
    const double alpha_norm(a*exp(-t*x_norm) + c);
    alpha = alpha_norm*y_std + y_mean;
    std::cout << "x:" << x << std::endl;
    std::cout << x_mean << " " << x_std << " " << y_mean << " " <<
      y_std <<
      " " << a << " " << t << " " << c << " " << x_norm << " " << alpha_norm <<
      " " <<  " " << alpha << std::endl;
    x_mean = (params_[loc_][std_at_mean_speed_]);
    x_std = (params_[loc_][std_at_mean_speed_+1]);
    y_mean = (params_[loc_][std_at_mean_speed_+2]);
    y_std = (params_[loc_][std_at_mean_speed_+3]);
    a = (params_[loc_][std_at_mean_speed_+4]);
    t = (params_[loc_][std_at_mean_speed_+5]);
    c = (params_[loc_][std_at_mean_speed_+6]);
    x_norm = ((x - x_mean)/x_std);
    const double loc_norm(a*exp(-t*x_norm) + c);
    loc = loc_norm*y_std + y_mean;
    std::cout << x_mean << " " << x_std << " " << y_mean << " " << y_std <<
      " " << a << " " << t << " " << c << " " << x_norm << " " << loc_norm <<
      " " <<  " " << loc << std::endl;

    unsigned N(params_[sigma_][std_at_mean_speed_]);
    for (unsigned i(1); i < N+1; ++i) {
      //p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
      scale += params_[sigma_][std_at_mean_speed_+i]*pow(x, N-i);
    }
    std::cout << "alpha:" << alpha << " loc:" << loc << " scale:" << scale <<
      " mean_speed:" << x << " std:" << 
      loc + gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale) << std::endl;

  }
  else {
    alpha = (params_[alpha_][col]);
    loc = (params_[loc_][col]);
    scale = (params_[sigma_][col]);
  }
  return loc + gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale);
}
*/


/*
// zero correlation with gamma instant fit
const double Walker::get_new_agent_speed(const unsigned agent_id) {
  double agent_speed(0);
  if (!agent_prev_speeds_[agent_id].size()) {
    agent_speed = agent_params_[instant_speed_][agent_id];
    agent_prev_speeds_[agent_id].push_back(agent_speed);
  }
  else {
    agent_speed = draw_gamma_distribution_at_col(instant_speed_, agent_id);
    agent_prev_speeds_[agent_id][0] = agent_speed;
  }
  return agent_speed;
}
*/


/*
const double Walker::get_new_agent_speed(const unsigned agent_id) {
  const double mean_speed(get_new_mean_speed(agent_id));
  const double std_mean_speed(get_new_std_mean_speed(agent_id, mean_speed));
  double agent_speed(gsl_ran_gaussian(get_model().get_gsl_rng(),
                                           std_mean_speed) + mean_speed);
  std::cout << "agent id:" << agent_id << " mean speed:" << mean_speed
    << " std mean speed:" << std_mean_speed << " speed:" << agent_speed <<
    std::endl;
  // units in um/min
  while (agent_speed > speed_max_ || agent_speed < 0) {
    agent_speed = (gsl_ran_gaussian(get_model().get_gsl_rng(),
                                           std_mean_speed) + mean_speed);
  }
  update_new_agent_speed(agent_prev_speeds_[agent_id], agent_speed);
  return agent_speed;
}
*/

/*
// uncorrelated agent speed
const double Walker::get_new_agent_speed(const unsigned agent_id) {
  const double alpha(2.16);
  const double loc(0.00);
  const double scale(2.78);
  speed_max_ = (25);
  double agent_speed(loc + 
                    gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale));
  while (agent_speed > speed_max_ || agent_speed < 0) {
    agent_speed = loc + gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale);
  }
  return agent_speed;
}
*/

/*
//correct formula but not working autocorrelation time
const double Walker::get_new_agent_speed(const unsigned agent_id) {
  double agent_speed(0);
  if (!agent_prev_speeds_[agent_id].size()) {
    const double alpha(2.16);
    const double loc(0.0);
    const double scale(2.78);
    speed_max_ = (25);
    const double std(0.5);
    double mean_speed(loc + 
                     gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale));
    while (mean_speed > speed_max_ || mean_speed < 0) {
      mean_speed = loc + gsl_ran_gamma(get_model().get_gsl_rng(), alpha, scale);
    }
    double agent_speed(gsl_ran_gaussian(get_model().get_gsl_rng(),
                                       std*mean_speed) + mean_speed);
    while (agent_speed > speed_max_ || agent_speed < 0) {
     agent_speed = gsl_ran_gaussian(get_model().get_gsl_rng(),
                                    std*mean_speed) + mean_speed;
    }
    agent_prev_speeds_[agent_id].push_back(agent_speed);
    agent_mean_speed_[agent_id] = mean_speed;
  }
  else {
    //double track_mean_speed(agent_sum_speed_[agent_id]/
    //                       agent_speed_cnt_[agent_id]);
    const double prev_speed(agent_prev_speeds_[agent_id][0]);
    const double agent_mean_speed(agent_mean_speed_[agent_id]);
    const double mean_speed(6.02);
    //const double variance(pow(prev_speed-mean_speed, 2));
    const double variance(16.56);
    double mean_instant_speed(variance*((1-0.18)*exp(-30/22.5)+0.18)/fabs(prev_speed-mean_speed) + mean_speed);
    double agent_variance(variance*((1-0.18)*exp(-30/22.5)+0.18));
    agent_speed = gsl_ran_gaussian(get_model().get_gsl_rng(),
                 agent_variance) + mean_instant_speed;
    std::cout << "agent speed:" << agent_speed << " mean_instant_speed:" << mean_instant_speed << " prev_speed:" << prev_speed << std::endl;
    while (agent_speed > speed_max_ || agent_speed < 0) {
      agent_speed = gsl_ran_gaussian(get_model().get_gsl_rng(),
                 agent_variance) + mean_instant_speed;
    }
    agent_prev_speeds_[agent_id][0] = agent_speed;
  }
  //agent_sum_speed_[agent_id] += agent_speed;
  //agent_speed_cnt_[agent_id] += 1;
  return agent_speed;
}
*/

/*
void Walker::update_correlated_orientation(const unsigned agent_id,
                                           const double agent_speed) {
  double roll(-20*M_PI/180);
  const Vector<double> x_axis(1, 0, 0);
  Quaternion rotateQ(Quaternion::represent_rotation(roll, x_axis.x, x_axis.y, x_axis.z)); 
  double pitch(-90*M_PI/180);
  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= 0.5) {
    pitch = -pitch;
  }

  const Vector<double> y_axis(0, 1, 0);
  Quaternion pitchQ(Quaternion::represent_rotation(pitch, y_axis.x, y_axis.y, y_axis.z));
  pitchQ = rotateQ.multiply(pitchQ);
  Quaternion& orientation(agent_orientations_[agent_id]); 
  Quaternion inv(orientation.inverse());
  orientation = orientation.multiply(pitchQ).normalise();
  q_ = q_.multiply(pitchQ).normalise();
  Quaternion P(inv.multiply(orientation).normalise());
  double r(atan2(P.z, P.y));
  double final_r(atan2(sin(2*r), cos(2*r)));
  final_r *= 180/M_PI;
  double p(atan2(P.y, P.w));
  double final_p(atan2(sin(2*p), cos(2*p)));
  final_p *= 180/M_PI;
  std::cout << "roll:" << final_r << " real_roll:" <<  roll*180/M_PI << std::endl;
  std::cout << "pitch:" << final_p << " real_pitch:" <<  pitch*180/M_PI << std::endl;
}
*/

/*
void Walker::update_correlated_orientation(const unsigned agent_id,
                                           const double agent_speed) {
  double roll(draw_poly_gamma_distribution_at_col(next_roll_, agent_speed,
                                                  360));
  if (gsl_ran_flat(get_model().get_gsl_rng(), 0, 1) <= 0.5) {
    roll = 360-roll;
  }
  roll *= get_interval()/60.0*M_PI/180;
  // roll as a quaternion.
  const Vector<double> x_axis(1, 0, 0);
  Quaternion rotateQ(Quaternion::represent_rotation(roll,
                                             x_axis.x, x_axis.y, x_axis.z)); 

  Quaternion& orientation(agent_orientations_[agent_id]); 
  // multiply orientation by rotateQ, because rotateQ is calculated relative
  // to cell, not in absolute space. 
  // alter the cell's orientation. 
  orientation = orientation.multiply(rotateQ).normalise();


  double pitch(draw_poly_gamma_distribution_at_col(next_pitch_, agent_speed,
                                                  360)*M_PI/180);
  pitch *= get_interval()/60.0;		// account for timestep.
  const Vector<double> y_axis(0, 1, 0);
  Quaternion pitchQ(Quaternion::represent_rotation(pitch,
                                              y_axis.x, y_axis.y, y_axis.z));
  // multiply orientation by rotateQ, because pitchQ is calculated relative
  // to cell, not in absolute space.
  orientation = orientation.multiply(pitchQ).normalise();
}
*/

/*
const Vector<double> Walker::get_displacement(const unsigned agent_id) {
  const Vector<double> p(
        dynamic_cast<MicroSpecies*>(&motile_species_)->get_relative_positions()[
    agent_id]);
  const double agent_speed(get_new_agent_speed(agent_id));
  Quaternion orientation(get_updated_agent_orientation(agent_id, agent_speed));
  const Vector<double> x_axis(1, 0, 0);
  Vector<double> facing(orientation.transform(x_axis));		
  double disp(agent_speed*get_interval()/60.0); //convert interval in s to min
  disp = fabs(disp);
  Vector<double> nex(double(facing.x*disp), double(facing.y*disp),
                     double(facing.z*disp));			
  return nex;
}
*/

/*
const Vector<double> Walker::get_displacement(const unsigned agent_id) {
  const double agent_speed(get_new_agent_speed(agent_id));
  //const double theta(gsl_ran_flat(get_model().get_gsl_rng(), -M_PI, M_PI));
  //const double phi(gsl_ran_flat(get_model().get_gsl_rng(), 0, M_PI));
  const double theta(45*M_PI/180);
  const double phi(0*M_PI/180);
  //const double r(agent_speed*get_interval()/60.0); //interval in s to min
  const double r(5);
  const Vector<double> disp(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi),
                           r*cos(theta));
  std::cout << "disp:" << disp.x << " " << disp.y << " " << disp.z << std::endl;
  return disp;
}
*/

/*
//corkscrew
void Walker::update_correlated_orientation(const unsigned agent_id,
                                           const double agent_speed) {
  //double roll(gsl_ran_flat(get_model().get_gsl_rng(), -180, 180)*M_PI/180);
  double roll(175*M_PI/180);
  const Vector<double> x_axis(1, 0, 0);
  Quaternion rotateQ(Quaternion::represent_rotation(roll, x_axis.x, x_axis.y,
                                                   x_axis.z)); 
  double pitch(5*M_PI/180);
  //double rand(gsl_ran_flat(get_model().get_gsl_rng(), 0, 1));
  //if (rand < 0.5) {
  //  pitch = -pitch;
  //}
  const Vector<double> y_axis(0, 1, 0);
  Quaternion pitchQ(Quaternion::represent_rotation(pitch, y_axis.x, y_axis.y,
                                                  y_axis.z));
  //pitchQ = rotateQ.multiply(pitchQ);
  Quaternion& orientation(agent_orientations_[agent_id]); 
  //Quaternion inv(orientation.inverse());
  orientation = orientation.multiply(rotateQ).normalise();
  orientation = orientation.multiply(pitchQ).normalise();
  //q_ = q_.multiply(pitchQ).normalise();
  //Quaternion P(inv.multiply(orientation).normalise());
  //double r(atan2(P.z, P.y));
  //double final_r(atan2(sin(2*r), cos(2*r)));
  //final_r *= 180/M_PI;
  //double p(atan2(P.y, P.w));
  //double final_p(atan2(sin(2*p), cos(2*p)));
  //final_p *= 180/M_PI;
  //std::cout << "roll:" << final_r << " real_roll:" <<  roll*180/M_PI << 
  //  std::endl;
  //std::cout << "pitch:" << final_p << " real_pitch:" <<  pitch*180/M_PI << 
  //  std::endl;
}
*/

