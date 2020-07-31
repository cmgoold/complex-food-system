/*
 * Implementation of the RK4 class
 */

#include "rk4.hpp"
#include "utilities.hpp"
#include <vector>

RK4::RK4(){ }

RK4::RK4(double& dt, CFSModel& cfs_model, std::vector<double>& initial_states){
  dt_ = dt;
  cfs_model_ = std::make_unique<CFSModel>(cfs_model);
  states_ = initial_states;
}

std::vector<double> RK4::integrate(){

  // vectors (1 per state variable) for the intermediate steps
  std::vector<double> K1, K2, K3, K4;
  // number of states
  int n_states = states_.size();
  // new state variables
  std::vector<double> new_states;
  double half = 0.5;

  for(int n = 0; n < n_states; ++n){
    // same as forward Euler integration
    K1.push_back(dt_ * cfs_model_->derivatives(states_)[n]);
  }

  for(int n = 0; n < n_states; ++n){
    // second step -- now with states + 0.5*K1
    // note: vector addition and vector-scalar multiplication
    std::vector<double> K1_states = vector_scalar(K1, half);
    K2.push_back(dt_ * cfs_model_->derivatives( vector_add(states_, K1_states) )[n]);
  }

  for(int n = 0; n < n_states; ++n){
    // third step -- now with states + 0.5*K2
    std::vector<double> K2_states = vector_scalar(K2, half);
    K3.push_back(dt_ * cfs_model_->derivatives( vector_add(states_, K2_states) )[n]);
  }

  for(int n = 0; n < n_states; ++n){
    // final step -- with states + K3
    K4.push_back(dt_ * cfs_model_->derivatives( vector_add(states_, K3))[n]);
  }

  // update new states
  double constant = 1.0/6.0;
  for(int n = 0; n < n_states; ++n){
    new_states.push_back(states_[n] + constant * (K1[n] + (K2[n] + K3[n]) * 2.0 + K4[n]));
  }

  // update the private data variable
  states_ = new_states;

  return states_;

}
