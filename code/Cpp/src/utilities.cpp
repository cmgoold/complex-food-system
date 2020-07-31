/*
 * Implementation of utility functions
 */

#include "utilities.hpp"
#include <vector>

std::vector<double> vector_add(std::vector<double>& v1, std::vector<double>& v2){
    std::vector<double> out;
    int v_size = v1.size();

    for(int i = 0; i < v_size; ++i){
      out.push_back( v1[i] + v2[i]);
    }

    return out;
}

std::vector<double> vector_scalar(std::vector<double>& v, double& x){
  std::vector<double> out;
  int v_size = v.size();

  for(int i = 0; i < v_size; ++i){
    out.push_back( x * v[i] );
  }

  return out;
}
