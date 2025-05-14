#pragma once

#include <iostream>

#include "globals.h"

namespace Obstacles
{
class Circles
{
public:

    int n_obstacles_;
    int n_ghosts_;
    int n_tot_;

    //double total_surface_;
    std::vector<double> x_;
    std::vector<double> y_;

    std::vector<double> radius_;
    //double *radius_sq_;
    double max_radius_;

    static constexpr int sides_ = 1;
    static constexpr int type_ = 0;

public:

    //constructors destructors:

    Circles()
    {
        n_obstacles_ = 0;
        n_ghosts_ = 0;
        n_tot_ = 0;
        max_radius_ = 0.0;
    }

    Circles ( int n_obstacles, int n_ghosts )
    {
        n_obstacles_ = n_obstacles;
        if ( n_obstacles > 1500000000 )
        {
            std::cout << "too many obsatcles ... ... ... aborting" << std::endl;
            std::terminate();
        }

        n_ghosts_ = n_ghosts;
        n_tot_ = n_ghosts + n_obstacles;
        //total_surface_ = 0.0;

        x_.resize ( n_tot_ );
        y_.resize ( n_tot_ );
        radius_.resize ( n_tot_ );
        max_radius_ = 0.0;

        double mem_size = sizeof ( double ) * n_tot_ * 4.0;
        //GlobalCounters::memsize += mem_size;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Obstacles (Circles)---" << std::endl;
            //std::cout << "radius[0] " << radius_[0] << std::endl;
            std::cout << "n Obstacles " << n_obstacles_ << std::endl;
            std::cout << "n Ghosts " << n_ghosts_ << std::endl;

            printf ( "total heapsize of the Obstacles is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    void setOne ( double x, double y, double r, int obs_nr )
    {
        x_[obs_nr] = x;
        y_[obs_nr] = y;
        radius_[obs_nr] = r;
        //total_surface_ += M_PI * r * r;
        max_radius_ = r > max_radius_ ? r : max_radius_;
    }

    void set ( int n, double x, double y )
    {
        x_[n] = x;
        y_[n] = y;
    }

    void setSizesBinary ( double r1, double r2, int n1, int n2 )
    {
        for ( auto i = 0; i < n1; ++i )
        {
            radius_[i] = r1;
        }

        for ( auto i = n1; i < n2; ++i )
        {
            radius_[i] = r2;
        }
        max_radius_ = r1 > r2 ? r1 : r2;

        if ( n1==0 ) max_radius_ = r2;
        if ( n2==0 ) max_radius_ = r1;
    }
    
    void setSizes ( double r1)
    {
        for ( auto i = 0; i < n_tot_; ++i )
        {
            radius_[i] = r1;
        }
        max_radius_ = r1;
    }
    

    inline bool isInside ( int *obstacles_list, int last, double x, double y )
    {
        int n = 0;
        for ( auto i = 0; i < last; i++ )
        {
            n = obstacles_list[i];
            if ( ( x - x_[n] ) * ( x - x_[n] ) + ( y - y_[n] ) * ( y - y_[n] ) < radius_[n] * radius_[n] )
            {
                return true;
            }
        }
        return false;
    }

    inline bool isInside ( int n, double x, double y )
    {
        if ( ( x - x_[n] ) * ( x - x_[n] ) + ( y - y_[n] ) * ( y - y_[n] ) < radius_[n] * radius_[n] )
        {
            return true;
        }
        return false;
    }

    inline int overlapp ( int nr_obs1, double shift_x, double shift_y )
    {
        for ( int j = 0 ; j < nr_obs1 ; j++ )
        {
            double dx = x_[j] + shift_x - x_[nr_obs1];
            double dy = y_[j] + shift_y - y_[nr_obs1];
            double dr = radius_[j] + radius_[nr_obs1];

            if ( dx*dx + dy*dy < dr*dr )
            {
                return -1;
            }
        }
        return 1;
    }

    void writeObstacles ( std::string filename, bool append )
    {
        std::ofstream obstaclesFile;

        if ( append == true )
        {
            obstaclesFile.open ( filename, std::ofstream::app );
        }
        else
        {
            obstaclesFile.open ( filename );
        }

        obstaclesFile << type_<< std::endl;
        obstaclesFile << n_tot_ << std::endl;
        obstaclesFile << 1 << std::endl;
        obstaclesFile.precision ( 8 );

        for ( int i = 0; i < n_tot_; ++i )
        {
            obstaclesFile << radius_[i] << ' ' << x_[i] << ' ' << y_[i] << std::endl;
        }
        obstaclesFile.close();
    }



};
};
