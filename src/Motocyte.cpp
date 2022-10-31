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

#include <iostream> 
#include <boost/date_time/posix_time/posix_time.hpp>
#include <Motocyte.hpp>
#include <Model.hpp>
#include <VisualWriter.hpp>
#include <NumberWriter.hpp>
#include <TrackReader.hpp>
#include <XMLModelReader.hpp>
#include <XMLModelWriter.hpp>
#include <TrackSpecies.hpp>
#include <TrackIDSpecies.hpp>
#include <PointSpecies.hpp>
#include <Reaction.hpp>
#include <Populator.hpp>





/*
average cells per frame:340.964
dims:1594.17 454.108 60.9877
min:225.971 1.00017 1.78676
max:1820.14 455.108 62.7745
*/


/*
//exp
//correct dimension, time and concentration with exp
int main() {
  const double interval(30.0);
  const Vector<double> dims(5000, 5000, 5000);
  const double agent_radius(5);
  Model model(dims, agent_radius);
  //VolumeSpecies A("A", model, 1600*3);
  VolumeSpecies A("A", model, 470); //2015 tracks
  //A.get_walker().set_walk_type(correlated_walk_);
  //A.get_walker().set_walk_type(exp_window_vmf2_intercept_fit_);
  //A.get_walker().set_walk_type(brownian_walk_);
  //A.get_walker().set_walk_type(exp_model_d_vs_t_);
  //A.get_walker().set_walk_type(exp_model_d_vs_s_);
  //A.get_walker().set_walk_type(exp_model_c_vs_t_);
  //A.get_walker().set_walk_type(exp_model_c_vs_s_);
  //A.get_walker().set_walk_type(exp_model_b_vs_t_);
  //A.get_walker().set_walk_type(exp_model_b_vs_s_);
  //A.get_walker().set_walk_type(exp_model_a_);
  //A.get_walker().set_walk_type(exp_window_vmf2_von_mises_intercept_fit_);
  //A.get_walker().set_walk_type(window_vmf2_intercept_walk_3D_);
  //A.get_walker().set_walk_type(brownian_walk_);
  A.get_walker().set_interval(interval);
  model.get_visual_writer().set_interval(interval);
  model.get_visual_writer().add(A);
  model.get_track_writer().set_skip_frames(25);
  model.get_track_writer().set_interval(interval);
  model.get_track_writer().set_padding(10);
  model.get_track_writer().add(A);
  model.initialize();
  model.step(249);
  return 0;
}
*/


//synthetic
int main() {
  const double interval(30.0);
  const Vector<double> dims(5000, 5000, 5000);
  const double agent_radius(5);
  Model model(dims, agent_radius);
  //VolumeSpecies A("A", model, 1600*3);
  VolumeSpecies A("A", model, 300); //2015 tracks
  //A.get_walker().set_walk_type(corkscrew_);
  A.set_check_collision(false); //no volume exclusion
  //A.get_walker().set_walk_type(synthetic_disp_auto_analytical_);
  //A.get_walker().set_walk_type(synthetic_window_walk_3D_fit_);
  //A.get_walker().set_walk_type(synthetic_window_walk_3D_);
  //A.get_walker().set_walk_type(window_vmf2_intercept_walk_3D_);
  //A.get_walker().set_walk_type(vmf2_intercept_walk_3D_);
  //A.get_walker().set_walk_type(walk_3D_);
  //A.get_walker().set_walk_type(walk_2D_);
  //A.get_walker().set_walk_type(brownian_walk_);
  A.get_walker().set_walk_type(brownian_uniform_turn_);
  //A.get_walker().set_walk_type(correlated_walk_);
  //A.get_walker().set_walk_type(uniform_directed_walk_);
  //A.get_walker().set_walk_type(random_directed_walk_);
  //A.get_walker().set_walk_type(brownian_walk_);
  A.get_walker().set_interval(interval);
  model.get_visual_writer().set_interval(interval);
  model.get_visual_writer().add(A);
  model.get_track_writer().set_interval(interval);
  model.get_track_writer().set_padding(10);
  model.get_track_writer().add(A);
  model.initialize();
  model.step(250);
  return 0;
}

/*
int main() {
  const double interval(30.0);
  const Vector<double> dims(1000, 1000, 1000);
  const double agent_radius(5);
  Model model(dims, agent_radius);
  //VolumeSpecies A("A", model, 1600*3);
  VolumeSpecies A("A", model, 10000); //2015 tracks
  A.get_walker().set_walk_type(brownian_walk_);
  A.get_walker().set_interval(interval);
  model.get_visual_writer().set_interval(interval);
  model.get_visual_writer().add(A);
  model.get_track_writer().set_interval(interval);
  model.get_track_writer().set_padding(100);
  model.get_track_writer().add(A);
  model.initialize();
  model.step(100);
  return 0;
}
*/




/*
//Brownian walk
int main() {
  std::string filename("model.xml");
  XMLModelWriter xml_model(filename);

  const Vector<double> dims(394.172, 394.752, 69.947);
  const double agent_radius(5);
  Model model(dims, agent_radius);
  VolumeSpecies A("A", model, 44, 0.115);
  A.get_walker().set_interval(20.63);
  PointSpecies As("As", A, 600);
  model.get_visual_writer().set_interval(20.63);
  model.get_visual_writer().add(A);
  model.get_visual_writer().add(As);
  model.get_track_writer().set_interval(20.63);
  model.get_track_writer().add(A);
  model.initialize();
  model.step(180);
  return 0;
}
*/

/*
int main() {
  std::string filename("20.628_Position.csv");
  TrackReader tracks(filename);
  const Vector<double> dims(tracks.get_dimensions());
  const double agent_radius(5);
  Model model(dims, agent_radius);
  TrackSpecies A("A", model, tracks);
  PointSpecies Ap("Ap", A, 50);
  
  MesoSpecies Am("Am", model, 0);
  Reaction exocytosis(Ap, 20, Am); //Ap -> Am
  exocytosis.set_interval(4);

  model.get_visual_writer().set_interval(tracks.get_interval());
  model.get_visual_writer().add(A);
  model.get_visual_writer().add(Ap);
  model.get_visual_writer().add(Am);

  model.initialize();
  //model.step(50);
  model.run(tracks.get_end_time());
  return 0;
}
*/

/*
int main() {
  std::string filename("20.628_Position.csv");
  TrackReader tracks(filename);
  const Vector<double> dims(tracks.get_dimensions());
  const double agent_radius(5);
  Model model(dims, agent_radius);
  TrackSpecies A("A", model, tracks);
  PointSpecies Ap("Ap", A, 200);

  MesoSpecies Ac("Ac", model, 1, 0);
  Reaction scavenging(Ap, Ac, 50, Ap); //Ap + Ac -> Ap
  scavenging.set_interval(5);

  model.get_visual_writer().set_interval(tracks.get_interval());
  model.get_visual_writer().add(A);
  model.get_visual_writer().add(Ap);
  model.get_visual_writer().add(Ac);

  model.initialize();
  //model.step(50);
  boost::posix_time::ptime start(
                     boost::posix_time::microsec_clock::universal_time()); 
  model.run(tracks.get_end_time());
  boost::posix_time::ptime end(
      boost::posix_time::microsec_clock::universal_time());
  boost::posix_time::time_duration duration(end-start);
  std::cout << "duration:" << duration.total_milliseconds()/1000.0 << std::endl; 
  return 0;
}
*/


/*
int main() {
  std::string filename("uniform_space_time_unfiltered_20.628_Position.csv");
  TrackReader tracks(filename);
  const Vector<double> dims(tracks.get_dimensions());
  const double agent_radius(5);
  Model model(dims, agent_radius);

  const unsigned track_size(tracks.get_size());
  model.get_visual_writer().set_interval(tracks.get_interval());
  model.get_number_writer().set_interval(0.1);

  std::vector<TrackIDSpecies*> A;
  std::vector<PointSpecies*> Ap;
  std::vector<MesoSpecies*> Am;
  std::vector<Reaction*> exocytosis;
  A.resize(track_size);
  Ap.resize(track_size);
  Am.resize(track_size);
  exocytosis.resize(track_size);
  for (unsigned i(0); i < track_size; ++i) {
    std::stringstream name_A;
    name_A << "A" << i;
    A[i] = new TrackIDSpecies(name_A.str(), model, tracks, i);
    name_A << "p";
    Ap[i] = new PointSpecies(name_A.str(), *A[i], 50);
    std::stringstream name_Am;
    name_Am << "Am";
    Am[i] = new MesoSpecies(name_Am.str(), model, 0);
    exocytosis[i] = new Reaction(*Ap[i], 20, *Am[i], 1); //Ap -> Am
    exocytosis[i]->set_interval(4);
    model.get_visual_writer().add(*A[i]);
    model.get_visual_writer().add(*Ap[i]);
    //model.get_visual_writer().add(*Am[i]); //for trail of agents
    //model.get_number_writer().add(*A[i]); //for concentration of agents
    //model.get_number_writer().add(*Am[i]); //for trail of agents
  }

  MesoSpecies Amo("Amo", model, 0);
  Reaction overlap(Ap, Am, 1, Amo); //Ap + Am -> Ao
  overlap.set_interval(0.2);
  model.get_visual_writer().add(Amo);
  model.get_number_writer().add(Amo); //for trail of agents

  model.initialize();

  boost::posix_time::ptime start(
                     boost::posix_time::microsec_clock::universal_time()); 
  model.run(tracks.get_end_time());
  boost::posix_time::ptime end(
      boost::posix_time::microsec_clock::universal_time());
  boost::posix_time::time_duration duration(end-start);
  std::cout << "duration:" << duration.total_milliseconds()/1000.0 << std::endl; 
  return 0;
}
*/

/*
//average cells per frame:340.964
//dims:1594.17 454.108 60.9877
//min:225.971 1.00017 1.78676
//max:1820.14 455.108 62.7745

//correct dimension, time and concentration with exp
int main() {
  const double interval(30.0);
  const Vector<double> dims(1620, 480, 85);
  const double agent_radius(5);
  Model model(dims, agent_radius);
  //VolumeSpecies A("A", model, 1600*3);
  VolumeSpecies A("A", model, 470); //2015 tracks
  A.get_walker().set_walk_type(correlated_walk_);
  //A.get_walker().set_walk_type(brownian_walk_);
  A.get_walker().set_interval(interval);
  model.get_visual_writer().set_interval(interval);
  model.get_visual_writer().add(A);
  model.get_track_writer().set_interval(interval);
  model.get_track_writer().set_padding(10);
  model.get_track_writer().add(A);
  model.initialize();
  model.step(249);
  return 0;
}
*/


/*
//synthetic experiments with large volume
int main() {
  const double interval(30.0);
  const Vector<double> dims(5000, 5000, 5000);
  const double agent_radius(5);
  Model model(dims, agent_radius);
  //VolumeSpecies A("A", model, 1600*3);
  VolumeSpecies A("A", model, 300); //2015 tracks
  //A.get_walker().set_walk_type(corkscrew_);
  A.set_check_collision(false); //no volume exclusion
  A.get_walker().set_walk_type(window_vmf2_intercept_walk_3D_);
  //A.get_walker().set_walk_type(vmf2_intercept_walk_3D_);
  //A.get_walker().set_walk_type(walk_3D_);
  //A.get_walker().set_walk_type(walk_2D_);
  //A.get_walker().set_walk_type(brownian_walk_);
  //A.get_walker().set_walk_type(correlated_walk_);
  //A.get_walker().set_walk_type(uniform_directed_walk_);
  //A.get_walker().set_walk_type(random_directed_walk_);
  //A.get_walker().set_walk_type(brownian_walk_);
  A.get_walker().set_interval(interval);
  model.get_visual_writer().set_interval(interval);
  model.get_visual_writer().add(A);
  model.get_track_writer().set_interval(interval);
  model.get_track_writer().set_padding(10);
  model.get_track_writer().add(A);
  model.initialize();
  model.step(100);
  return 0;
}
*/
