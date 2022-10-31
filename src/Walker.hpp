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


#ifndef __Walker_hpp
#define __Walker_hpp

#include <fstream>
#include <random>
#include <Common.hpp>
#include <Process.hpp>
#include <Quaternion.hpp>

typedef const Vector<double> (Walker::*DisplaceMethod)(const unsigned);

struct Param {
};

class Walker: public Process { 
public: 
  Walker(const double, Model& model, MotileSpecies&, const double interval);
  ~Walker() {}
  const Vector<double> get_displacement(const Vector<double> agent_xyz,
                                        const unsigned agent_idx);
  void set_homogeneous_crw() {};
  double step();
  std::vector<Vector<double>>& get_agent_redirections();
  std::vector<unsigned>& get_agent_redirection_magnitudes();
  virtual void set_interval(double);
  virtual void initialize();
  void initialize_agent_parameters();
  void remove_agent(const unsigned agent_idx);
  void add_agent();
  void set_walk_type(const unsigned walk_type);
protected:
  const Vector<double> get_exp_model_a_displacement(
                                        const Vector<double> agent_xyz,
                                        const unsigned agent_id);
  const Vector<double> get_exp_model_b_vs_s_displacement(
                                        const Vector<double> agent_xyz,
                                        const unsigned agent_id);
  const Vector<double> get_exp_model_b_vs_t_displacement(
                                        const Vector<double> agent_xyz,
                                        const unsigned agent_id);
  const Vector<double> get_exp_model_c_vs_s_displacement(
                                        const Vector<double> agent_xyz,
                                        const unsigned agent_id);
  const Vector<double> get_exp_model_c_vs_t_displacement(
                                        const Vector<double> agent_xyz,
                                        const unsigned agent_id);
  const Vector<double> get_exp_model_d_vs_s_displacement(
                                        const Vector<double> agent_xyz,
                                        const unsigned agent_id);
  const Vector<double> get_exp_model_d_vs_t_displacement(
                                        const Vector<double> agent_xyz,
                                        const unsigned agent_id);
  void get_vectors_v_and_y(const Vector<double> agent_xyz,
                           const unsigned agent_id, Vector<double>& v,
                           Vector<double>& y);
  const Vector<double> get_synthetic_window_walk_3D_fit_displacement(
                                              const Vector<double> agent_xyz,
                                              const unsigned agent_id);
  const Vector<double> get_synthetic_window_walk_3D_displacement(
                                              const Vector<double> agent_xyz,
                                              const unsigned agent_id);
  void update_new_agent_xyz(std::vector<Vector<double>>& prev_xyz,
                            const Vector<double> agent_xyz);
  const Vector<double> get_brownian_walk_displacement(const unsigned agent_id);
  const Vector<double> get_brownian_uniform_turn_displacement(
                                                const Vector<double> agent_xyz,
                                                const unsigned agent_id);
  const Vector<double> get_walk_2D_displacement(const unsigned agent_id);
  const Vector<double> get_walk_3D_displacement(const unsigned agent_id);
  const Vector<double> get_vmf2_intercept_walk_3D_displacement(
                                                     const unsigned agent_id);
  const Vector<double> get_exp_window_vmf2_intercept_fit_displacement(
                                                 const Vector<double> agent_xyz,
                                                 const unsigned agent_id);
  const Vector<double> get_exp_window_vmf2_von_mises_intercept_fit_displacement(
                                                 const Vector<double> agent_xyz,
                                                 const unsigned agent_id);
  const Vector<double> get_window_vmf2_intercept_walk_3D_displacement(
                                                 const Vector<double> agent_xyz,
                                                 const unsigned agent_id);
  const Vector<double> get_random_directed_walk_displacement(
                                                     const unsigned agent_id);
  const Vector<double> get_uniform_directed_walk_displacement(
                                                     const unsigned agent_id);
  const Quaternion get_homogeneous_crw_orientation(const unsigned agent_id);
  Vector<double> get_random_direction();
  const Quaternion get_updated_agent_orientation(const unsigned agent_id,
                                                 const double agent_speed);
  void update_brownian_orientation(const unsigned agent_id);
  void update_corkscrew_orientation(const unsigned agent_id);
  void update_correlated_orientation(const unsigned agent_id,
                                     const double agent_speed);
  void update_uncorrelated_orientation(const unsigned agent_id);
  const double get_new_agent_speed(const unsigned agent_id);
  void parse_file();
  double draw_stable_distribution_at_col(const unsigned col);
  double draw_gamma_distribution_at_col(const int col, const unsigned agent_id);
  void update_new_agent_speed(std::vector<double>& prev_speeds,
                              const double agent_speed);
  const double get_new_std_mean_speed(const unsigned agent_id,
                                     const double mean_speed);
  const double get_new_mean_speed(const unsigned agent_id);
  const double draw_poly_gamma_distribution_at_col(const unsigned col,
                                                  const double x,
                                                  const double max_value);
  const Vector<double> sample_vmf(const double kappa, const Vector<double> mu);
  const double sample_von_mises(const double mu, const double kappa);
  const double sample_weight(const double kappa, unsigned dim);
  const Vector<double> sample_orthonormal_to(const Vector<double> mu);
  const Vector<double> draw_poly_vmf_distribution_at_col(const unsigned col,
                                                        const double x,
                                                        const double prev_x);
  const double draw_roll_from_von_mises(const double agent_speed,
                                       const double turn);
  const double draw_roll_from_vector_von_mises(const double agent_speed,
                                              const double agent_turn,
                                              const double agent_roll);
  const double draw_next_speed_xcorrelated_gamma(const unsigned agent_id);
  const double draw_next_speed_xcorrelated_beta(const unsigned agent_id);
  const double draw_next_speed_vs_s_t_r_s2_gamma(const unsigned agent_id);
  const double draw_next_speed_vs_s_t_r_s2_beta(const unsigned agent_id);
  const double draw_next_speed_vs_s_beta(const unsigned agent_id);
  const double draw_next_speed_vs_none_gamma(const unsigned agent_id);
  const double draw_next_speed_vs_sp_gamma(const unsigned agent_id,
                                           double prev_speed);
  const double draw_roll_from_vector_von_mises_intercept();
  const double draw_next_roll_vs_none_von_mises();
  const double draw_next_roll_vs_s_von_mises(const double agent_speed);
  const double draw_next_roll_vs_t_von_mises(const double agent_turn);
  const Vector<double> draw_direction_from_vmf2_intercept();
  const Vector<double> draw_next_direction_vs_none_vmf();
  const Vector<double> draw_next_direction_vs_s_vmf(const double agent_speed);
  const Vector<double> draw_next_direction_vs_s_tp_vmf(const double agent_speed,
                                                       const double prev_turn);
  const Vector<double> draw_direction_from_vector_vmf(const double agent_speed,
                                                     const double agent_turn,
                                                     const double agent_roll);
  const Vector<double> draw_exp_poly_vmf_distribution_at_col(
                                                      const unsigned col,
                                                      const double prev_speed,
                                                      const double prev_turn,
                                                      const double prev_roll);
  const double get_agent_mean_speed(std::vector<double>& prev_speeds);

private:
  Quaternion q_ = Quaternion::face_vector(Vector<double>(3, 5, 5));
  std::vector<Vector<double>> facing_;
  double pitch_ = 0;
  double roll_ = -180;
  unsigned walk_type_ = brownian_walk_;
  std::string filename_;
  std::ifstream param_file_;
  int next_speed_ = -1;
  int next_direction_ = -1;
  int instant_speed_ = -1;
  int mean_speed_ = -1;
  int std_mean_speed_ = -1;
  int mean_pitch_rate_ = -1;
  int std_mean_pitch_rate_ = -1;
  int mean_roll_rate_ = -1;
  int std_mean_roll_rate_ = -1;
  int next_pitch_ = -1;
  int next_roll_ = -1;
  int std_at_mean_speed_ = -1;
  unsigned window_size_ = 1;
  unsigned std_at_mean_speed_size_ = 0;
  unsigned param_col_size_ = 0;
  const unsigned alpha_ = 0;
  const unsigned beta_ = 1;
  const unsigned loc_ = 2;
  const unsigned sigma_ = 3;
  const unsigned max_ = 4;
  double speed_max_ = 16.0;
  double std_speed_max_ = 8.0;
  std::vector<std::vector<double>> params_;
  std::vector<double> pitch_upper_speeds_;
  std::vector<double> roll_upper_speeds_;
  std::vector<std::vector<double>> agent_prev_speeds_;
  std::vector<std::vector<double>> agent_prev_turns_;
  std::vector<std::vector<double>> agent_prev_rolls_;
  std::vector<std::vector<Vector<double>>> agent_prev_xyz_;
  std::vector<double> agent_mean_speed_;
  std::vector<unsigned> agent_speed_cnt_;
  /*
  const double speed_mean_ = 5.15;
  const double speed_std_ = 8.0;
  const double max_speed_ = 25;
  const double roll_rate_mean_ = -1.0; 
  const double roll_rate_std_ = 0.0; 
  const double pitch_rate_mean_ = 0.5;
  const double pitch_rate_std_ = 0.35;
  */
  //Jorge's data
  const double speed_mean_mean_ = 8.0; //unit um/min
  //const double speed_mean_std_ = 0.1;
  const double speed_mean_std_ = 2.0; //unit um/min
  //const double speed_std_mean_ = 0.076;
  const double speed_std_mean_ = 0.00001;
  //const double speed_std_std_ = 0.02;
  const double speed_std_std_ = 0.00001;
  double D_;
  MotileSpecies& motile_species_;
  std::vector<Vector<double>> agent_redirections_;
  std::vector<std::vector<double>> agent_params_; 
  std::vector<unsigned> agent_redirection_magnitudes_;
  std::vector<Quaternion> agent_orientations_; 
};

#endif /* __Walker_hpp */

