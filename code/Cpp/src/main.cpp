/*
 * Estimate the parameters of the complex food systems model
 * using approximate Bayesian computation (ABC) Sequential Monte Carlo (SMC)
 */

#include "rk4.hpp"
#include <vector>
#include <iostream>

int main(){

  // create a new model
  int n_ts = 3;
  double dt = 0.01;

  std::vector<double> parameters;
  std::vector<double> initial_states;
  int n_states = 4;

  parameters.push_back(1.0/52.0); //a
  parameters.push_back(80.0); //b
  parameters.push_back(1/(2.5*52)); //e
  parameters.push_back(1/52.0); //f
  parameters.push_back(110.0*24*0.75); //g
  parameters.push_back(0.3); //w
  parameters.push_back(1.0); //s
  parameters.push_back(0.8); //k
  parameters.push_back(30e6); //h
  parameters.push_back(0.1); //m
  parameters.push_back(150.0); //q
  parameters.push_back(0.05); //r

  initial_states.push_back(100e3);
  initial_states.push_back(30e6);
  initial_states.push_back(30e6);
  initial_states.push_back(140.0);

  // make model
  CFSModel cfs_model(parameters);

  // solver
  RK4 solve(dt, cfs_model, initial_states);

  std::cout << "--------------------------------------" << std::endl;

  std::vector<double> out;
  int sim_t = n_ts * 1/dt;

  // run the model
  std::cout << "Solving..." << std::endl;
  std::cout << "0: ";
  for(int n = 0; n < n_states; ++n){
    std::cout << initial_states[n] << " ";
  }
  std::cout << std::endl;
  for(int t = 0; t < sim_t; ++t){
    out = solve.integrate();
    std::cout << t * dt + dt << ": ";
    for(int n = 0; n < n_states; ++n){
      std::cout << out[n] << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
