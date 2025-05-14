#include "globals.h"
#include <iostream>
#include <fstream>
#include <random>

bool GlobalFlags::verbose = false;

//double GlobalCounters::memsize = 0.0;


namespace MathsFunctionApprox
{
  double stirling(int n)
  {
    return n*log(n) - n + 0.5 * log ((6*n+1) * (M_PI/3.0));
  }
}


std::ofstream gl_tr_file;
std::ofstream gl_obs_file;

namespace Globalrand
{
    //globally accessible random number generator:
    
    //random device for seeds
    std::random_device rd;
    
    //seeding
    std::mt19937_64 e1(rd());
    std::mt19937_64 e2(rd());
    
    //manual seeding  REALLY DO NOT FORGET TO DIABLE THIS ON PRODUCTION !!!
    //std::mt19937_64 e1(0);
    //std::mt19937_64 e2(0);
}
