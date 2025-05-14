#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdint.h>

namespace Shapes
{
class Segment
{

  public:

    int n_tot_;

    double *x_;
    double *y_;

    double *length_;
    double *length_sq_;
    double *length_inv_;
    double max_length_;

  public:
    //constructors and destructors

    Segment()
    {
      n_tot_ = 0;
      x_ = NULL;
      y_ = NULL;
      length_ = NULL;
      length_sq_ = NULL;
      length_inv_ = NULL;
      max_length_ = 0.0;
    }

    Segment(int n_tot)
    {
      n_tot_ = n_tot;

      x_ = new double[n_tot_ + 1];
      y_ = new double[n_tot_ + 1];
      length_ = new double[n_tot_ + 1];
      length_sq_ = new double[n_tot_ + 1];
      length_inv_ = new double[n_tot_ + 1];
      max_length_ = 0.0;
    }

    ~Segment()
    {
      delete[] x_;
      delete[] y_;
      delete[] length_;
      delete[] length_inv_;
      delete[] length_sq_;
    }

//functions:

    void setVertex(double x, double y, int_fast32_t vertex_nr)
    {
      x_[vertex_nr] = x;
      y_[vertex_nr] = y;
    }

    void initPath()
    {
      double x, y;
      for (int_fast32_t i = 0; i < n_tot_; ++i)
      {
        x = x_[i + 1] - x_[i];
        y = y_[i + 1] - y_[i];
        length_[i] = std::sqrt(x * x + y * y);
        length_inv_[i] = 1.0 / length_[i];
        length_sq_[i] = length_[i] * length_[i];
      }
    }
};

}
