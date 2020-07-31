/*
 * Definition of the complex food system (cfs) model class
 */

#ifndef CFS_INCLUDE_MODEL_
#define CFS_INCLUDE_MODEL_

#include <vector>

class CFSModel{

    public:
      CFSModel();
      CFSModel(std::vector<double> parameters);

    private:
      double a_, b_, e_, f_, g_, w_, s_, k_, h_, m_, q_, r_;

    public:
      std::vector<double> derivatives(std::vector<double> states);
};

#endif
