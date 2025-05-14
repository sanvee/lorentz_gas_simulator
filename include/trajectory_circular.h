#pragma once

#include "globals.h"

#include<iostream>
#include<cmath>

#include "mytime.h"

//SoA layout of the particle trajectories.
//Two cases are considered: with magnetic field on and off.
//Circular(int nr_trajectories, double lifetime, double radius)
namespace Trajectory
{
class Circular // v = 1
{
public:

    //number of maximum circle segment
    double rx_;
    double ry_;

    double radius_;
    double radius_inv_;

    double rx_old_;
    double ry_old_;

    double lifetime_;
    int n_tot_coll_;

    TIME_type current_time_;

    //-------------------------------------------
    TIME_type last_time_;
    TIME_type time_correction_;

    double rx_unfolded_;
    double ry_unfolded_;

    double rx_unfolded_old_;
    double ry_unfolded_old_;

    double cx_begin_;
    double cy_begin_;

    //-------------------------------------------

    double theta_begin_;
    double theta_begin_old_;

    //coordinates of the center of the trajectory circles at t=0;
    double *rx_initial_;
    double *ry_initial_;
    double *cx_initial_;
    double *cy_initial_;

    //number of different trajectory starting points
    int n_trajectories_;
    bool written_ = false;

public:
    Circular()
    {
        rx_initial_ = NULL;
        ry_initial_ = NULL;
        cx_initial_ = NULL;
        cy_initial_ = NULL;
    }

    ///Circular ( int nr_trajectories, double lifetime, double radius )
    Circular ( int nr_trajectories, double lifetime, double radius )
    {
        n_trajectories_ = nr_trajectories;
        lifetime_ = lifetime;
        n_tot_coll_ = 0;

        radius_ = radius;
        radius_inv_ = 1.0 / radius;

        current_time_ = 0.0;
        last_time_ = 0.0;
        time_correction_ = 0.0;

        rx_ = 0.0;
        ry_ = 0.0;

        rx_old_ = 0.0;
        ry_old_ = 0.0;

        rx_unfolded_ = 0.0;
        ry_unfolded_ = 0.0;

        rx_unfolded_old_ = 0.0;
        ry_unfolded_old_ = 0.0;

        cx_begin_ = 0.0;
        cy_begin_ = 0.0;

        theta_begin_ = 0.0;
        theta_begin_old_ = 0.0;

        rx_initial_ = new double[nr_trajectories]();
        ry_initial_ = new double[nr_trajectories]();
        cx_initial_ = new double[nr_trajectories]();
        cy_initial_ = new double[nr_trajectories]();

        double mem_size = sizeof ( double ) * nr_trajectories * 4.;
        //GlobalCounters::memsize += mem_size;
        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Circular Trajectory---" << std::endl;
            std::cout << "Bfield " << radius_inv_ << std::endl;
            std::cout << "radius " << radius << std::endl;
            std::cout << "lifetime " << lifetime_ << std::endl;
            printf ( "total heapsize of the trajectories is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    ~Circular()
    {
        delete[] rx_initial_;
        delete[] ry_initial_;
        delete[] cx_initial_;
        delete[] cy_initial_;
    }

    void engage ( int n )
    {
        n_tot_coll_ = 0;
        rx_ = rx_initial_[n];
        ry_ = ry_initial_[n];

        rx_old_ = rx_;
        ry_old_ = ry_;

        rx_unfolded_ = rx_;
        ry_unfolded_ = ry_;

        rx_unfolded_old_ = rx_;
        ry_unfolded_old_ = ry_;

        cx_begin_ = cx_initial_[n];
        cy_begin_ = cy_initial_[n];

        theta_begin_ = atan2 ( cy_begin_ - ry_, cx_begin_ - rx_ );
        theta_begin_old_ = theta_begin_;

        current_time_ = 0.0;
        last_time_ = 0.0;
        time_correction_ = 0.0;
    }

    inline void get_current_position_and_velocity ( double time, double &cx_t, double &cy_t, double &vx_t, double &vy_t )
    {
        //actual angle at the time
        double theta_t = ( time - last_time_ ) * radius_inv_ + theta_begin_old_;

        cx_t = rx_unfolded_old_ + radius_ * cos ( theta_t );
        cy_t = ry_unfolded_old_ + radius_ * sin ( theta_t );

        vx_t = ( ry_unfolded_old_ - cy_t ) * radius_inv_;
        vy_t = ( cx_t - rx_unfolded_old_ ) * radius_inv_;
    }

    inline void get_current_position ( double time, double &cx_t, double &cy_t )
    {
        //actual angle at the time
        double theta_t = ( time - last_time_ ) * radius_inv_ + theta_begin_old_;

        cx_t = rx_unfolded_old_ + radius_ * cos ( theta_t );
        cy_t = ry_unfolded_old_ + radius_ * sin ( theta_t );
    }

    inline void get_current_position_and_velocity_inbox ( double time, double &cx_t, double &cy_t, double &vx_t, double &vy_t )
    {
        //actual angle at the time
        double theta_t = ( time - last_time_ ) * radius_inv_ + theta_begin_old_;

        cx_t = rx_old_ + radius_ * cos ( theta_t );
        cy_t = ry_old_ + radius_ * sin ( theta_t );

        vx_t = ( ry_old_ - cy_t ) * radius_inv_;
        vy_t = ( cx_t - rx_old_ ) * radius_inv_;
    }

    inline void normalizeAngle ( double &theta )
    {
        theta = theta - std::floor ( theta * Constants::TWO_PI_INV ) * Constants::TWO_PI;
    }

    inline void kahan_time_addition ( double addend )
    {
        last_time_ = current_time_;
        double y = addend - time_correction_;
        current_time_ = last_time_ + y;
        time_correction_ = ( current_time_ - last_time_ ) - y;
    }
    void writeStartpositions(std::string filename)
    {
                std::ofstream obstaclesFile;
        obstaclesFile.open ( filename );
                obstaclesFile.precision ( 8 );
        for(int i = 0; i < n_trajectories_; ++i)
        {
            obstaclesFile << cx_initial_[i] << ' ' << cy_initial_[i] << std::endl;
        }
        
    }

};
};
