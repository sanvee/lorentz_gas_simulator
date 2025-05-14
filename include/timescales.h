#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdint.h>

#include"globals.h"


class Timestamp
{
public:
    double absolute_time_;    // absolute time
    int step_offset_;         // intermediate starting times
    int relative_step_;       // steps from offset

    Timestamp()
    {
        absolute_time_ = 0.0;
        step_offset_ = 0;
        relative_step_ = 0;
    }
};

bool timestampSorter ( Timestamp const &lhs, Timestamp const &rhs )
{
    if ( lhs.absolute_time_ == rhs.absolute_time_ )
    {
        return lhs.relative_step_ < rhs.relative_step_;
    }
    else
    {
        return lhs.absolute_time_ < rhs.absolute_time_;
    }
}

namespace TimeScale
{
class LogScale
{
public:
    double n_averages_;
    int64_t n_stamps_;
    double t_end_;
    double t_begin_;

    std::vector<Timestamp> timestamps_;
    std::vector<double> time_;

    LogScale ( double t_begin, double t_end, int n_bin_log, int n_averages )
    {
        n_averages_ = n_averages;

        if ( t_begin <= 0.0 )
        {
            t_begin = 0.01;
        }
        t_begin_ = t_begin;
        t_end_ = t_end;

        double t_basis = std::pow ( t_end / t_begin, 1.0 / ( double ) n_bin_log );
        double t_avg_shift = ( t_end - t_begin ) /n_averages;

        int relative_step;
        double time = 0.0;
        double time_relative;
        double time_offset;

        time_.reserve ( n_bin_log );
        Timestamp stamp;

        for ( auto i = 0; i < n_averages; i++ )
        {
            time_offset = i * t_avg_shift;
            time_relative = 0.0;
            time = time_offset + time_relative;

            //first timestamp relative to 0
            relative_step = 0;

            stamp.absolute_time_ = time;
            stamp.relative_step_ = relative_step;
            stamp.step_offset_ = i;
            timestamps_.push_back ( stamp );

            if ( i == 0 )
            {
                time_.push_back ( time );
            }

            while ( true )
            {
                relative_step ++;
                time_relative = std::pow ( t_basis, relative_step ) * t_begin;
                time = time_offset + time_relative;

                if ( time_relative > t_end_ )
                {
                    break;
                }

                stamp.absolute_time_ = time;
                stamp.relative_step_ = relative_step;
                stamp.step_offset_ = i;
                timestamps_.push_back ( stamp );

                if ( i == 0 )
                {
                    //as we take the initial position at t_begin is 0
                    time_.push_back ( time );
                }
            }
        }
        n_stamps_ = timestamps_.size();
        std::sort ( timestamps_.begin(), timestamps_.end(), &timestampSorter );

        /*
        for ( int i = 0; i < timestamps_.size(); i++ )
        {
            std::cout << i << ' ' << timestamps_[i].step_offset_ << ' ' << timestamps_[i].relative_step_ << ' ' << timestamps_[i].absolute_time_ << std::endl;
        }
        */

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Timescale---" << std::endl;
            std::cout << "t_begin " << t_begin_ << std::endl;
            std::cout << "t_end " << t_end_ << std::endl;
            std::cout << "n_averages " << n_averages_ << std::endl;
            std::cout << "LogScale has " << time_.size() << " datapoints" << std::endl;
            std::cout << "LogScale has " << n_stamps_ << " timestamps" << std::endl;
            std::cout << "and " << sizeof ( Timestamp ) * n_stamps_ / 1000000.0 << " Mb heapsize" << std::endl << std::endl;
        }
    }
};

/*
  class LogLinScale
  {
    public:
      double n_averages_;
      int64_t n_stamps_;

      double t_begin_;
      double t_end_;

      std::vector<Timestamp> timestamps_;
      std::vector<double> time_;

      LogLinScale ( double t_begin, double t_log_end, int n_bins_log, double dt_max, double t_end, int n_averages )
      {
        n_averages_ = n_averages;

        if ( t_begin == 0 ) t_begin = 0.001;

        t_begin_ = t_begin;
        t_end_ = t_end;

        double t_basis_small = std::pow ( t_log_end / t_begin, 1.0 / ( double ) n_bins_log );
        double t_avg_shift = ( t_end - t_begin ) /n_averages;

        int partial_step;
        double time = 0.0;
        double time_relative;
        double time_offset;

        time_.reserve ( n_bins_log + std::ceil ( ( t_end - t_log_end ) / dt_max ) );
        Timestamp stamp;

        for ( auto i = 0; i < n_averages; i++ )
        {
          time_relative = 0.0;
          time_offset = i * t_avg_shift;
          time = time_offset + time_relative;

          partial_step = 0;

          stamp.absolute_time_ = time;
          stamp.relative_step_ = partial_step;
          stamp.step_offset_ = i;
          timestamps_.push_back ( stamp );

          if ( i == 0 ) time_.push_back ( time );

          while ( time < t_end )
          {
            if ( time_relative < t_log_end )
            {
              time_relative = std::pow ( t_basis_small, partial_step ) * t_begin;
              time = time_offset + time_relative;
              partial_step ++;
              stamp.absolute_time_ = time;
              stamp.relative_step_ = partial_step;
              stamp.step_offset_ = i;
              timestamps_.push_back ( stamp );
            }
            else
            {
              time_relative += dt_max;
              time = time_offset + time_relative;
              partial_step++;
              stamp.absolute_time_ = time;
              stamp.relative_step_ = partial_step;
              stamp.step_offset_ = i;
              timestamps_.push_back ( stamp );
            }
            if ( i == 0 ) time_.push_back ( time );
          }
        }
        n_stamps_ = timestamps_.size();
        std::sort ( timestamps_.begin(), timestamps_.end(), &timestampSorter );

        if ( GlobalFlags::verbose )
        {
          std::cout << std::endl << "---Timescale---" << std::endl;
          std::cout << "LogLinScale has " << time_.size() << " datapoints" << std::endl;
          std::cout << "LogLinScale has " << n_stamps_ << " timestamps" << std::endl;
          std::cout << "and " << sizeof ( Timestamp ) * n_stamps_ / 1000000.0 << " Mb heapsize" << std::endl << std::endl;
        }
      }
  };
*/

class LinScale
{
public:
    double t_begin_;
    double t_end_;
    double dt_;
    int n_averages_;
    int n_timesteps_;
    int n_partial_;

    LinScale ( double t_begin, double t_end, double dt, int n_averages )
    {
        t_begin_ = t_begin;
        t_end_ = t_end;

        n_timesteps_ = static_cast<int> ( ( t_end - t_begin ) / dt );
        n_averages_ = n_averages;
        n_partial_ = std::ceil ( static_cast<double> ( n_timesteps_ ) / n_averages_ );
        dt_ = dt;
    }
};
};
