#pragma once 

#include <cmath>
#include <iostream>
#include <fstream>
#include <random>

namespace GlobalFlags
{
  extern bool verbose;
};

namespace Constants
{
//Maths:
  constexpr double TWO_PI         = 2.0 * M_PI;
  constexpr double PI_INV         = M_1_PI;
  constexpr double TWO_PI_INV     = 1.0 / (2.0 * M_PI);

//Physics:
  constexpr double ELECTRON_MASS      = 9.1093837015e-31;
  constexpr double ELEMENTARY_CHARGE    = 1.602176634e-19;

  constexpr double H_BAR              = 1.054571817e-34;

  constexpr double M_FACTOR_GaAs      = 0.067;
};

namespace MathsFunctionApprox
{
  double stirling(int n);
};

namespace GlobalCounters
{
  //extern double memsize;
};

namespace Globalrand
{
    extern std::mt19937_64 e1;
    extern std::mt19937_64 e2;
}


extern std::ofstream gl_tr_file;
extern std::ofstream gl_obs_file;    


