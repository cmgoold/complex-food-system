/*
 * Implementation of the CFSModel class
 */

#include "cfs_model.hpp"

#include <vector>

CFSModel::CFSModel() { }

CFSModel::CFSModel(std::vector<double> parameters){
  a_ = parameters[0];
  b_ = parameters[1];
  e_ = parameters[2];
  f_ = parameters[3];
  g_ = parameters[4];
  w_ = parameters[5];
  s_ = parameters[6];
  k_ = parameters[7];
  h_ = parameters[8];
  m_ = parameters[9];
  q_ = parameters[10];
  r_ = parameters[11];
}

std::vector<double> CFSModel::derivatives(std::vector<double> states){

  double C = states[0];
  double I = states[1];
  double D = states[2];
  double P = states[3];

  // return vector
  std::vector<double> derivs;

  // C
  derivs.push_back(a_ * C * ( P/b_ - 1 ) - e_ * C);
  // I
  derivs.push_back(f_ * g_ * C - w_ * I - I*D/(D*s_ + I)  + k_ * (h_ - f_ * g_ * C));
  // D
  derivs.push_back(m_ * (h_ * q_/P - D));
  // P
  derivs.push_back(r_ * P * ( (s_*D)/I - 1));

  return derivs;
}
