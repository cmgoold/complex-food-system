/*
 * Definition of the Runge-Kutta 4th Order integration scheme
 */

#ifndef CFS_INCLUDE_RK4_
#define CFS_INCLUDE_RK4_

#include "cfs_model.hpp"
#include <vector>
#include <memory>

class RK4{

  public:
    // constructors
    RK4();
    RK4(double& dt, CFSModel& cfs_model, std::vector<double>& initial_states);

  private:
    // unique pointer (sole owner) of the ODE model
    std::unique_ptr<CFSModel> cfs_model_;
    // vector of state variables
    std::vector<double> states_;
    // integration interval
    double dt_;

  // member functions
  public:
    std::vector<double> integrate();


};

#endif
