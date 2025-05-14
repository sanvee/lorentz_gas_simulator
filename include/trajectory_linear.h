#pragma once

#include <iostream>

#include "container.h"


class LinearPropagator;
class LinearPropagator_0;
class WindtreePropagator;

namespace Trajectory
{


class Linear //v=1
{
    friend class ::LinearPropagator;
    friend class ::LinearPropagator_0;
    friend class ::WindtreePropagator;

public:

    double cx_;
    double cy_;

    double vx_;
    double vy_;

    double cx_unfolded_;
    double cy_unfolded_;

    double cx_unfolded_old_;
    double cy_unfolded_old_;

    //--Cacheline--------------------

    //last segment end point
    double cx_old_;
    double cy_old_;

    //last segment tangent vector (speed = 1)
    double vx_old_;
    double vy_old_;

    double lifetime_;
    
private:
    
    double current_time_;
    double time_correction_;
    double last_time_;
public:

    //-- Cacheline----------------------

    int n_trajectories_;
    double windtree_angle_;

    bool written_ = false;

    //starting positions
    double * cx_initial_;
    double * cy_initial_;
    double * vx_initial_;
    double * vy_initial_;

    Linear()
    {
        n_trajectories_ = 0.0;
        lifetime_ = 0.0;
        current_time_ = 0.0;
        last_time_ = 0.0;

        cx_ = cy_ = cx_old_ = cy_old_ = 0.0;
        cx_unfolded_ = cy_unfolded_ = cx_unfolded_old_ = cy_unfolded_old_ = 0.0;
        vx_ = vy_ = vx_old_ = vy_old_ = 0.0;
        //angle_ = angle_correction_ = 0.0;

        windtree_angle_ = 0.0;

        cx_initial_ =  NULL;
        cy_initial_ =  NULL;
        vx_initial_ =  NULL;
        vy_initial_ =  NULL;
        //angle_inital_= NULL;
    }

    Linear ( int nr_trajectories, double lifetime )
    {
        //kahanFile.open ( "kahan.txt");
        //kahanFile.precision ( 18 );

        n_trajectories_ = nr_trajectories;
        lifetime_ = lifetime;
        current_time_ = 0.0;
        last_time_ = 0.0;
        time_correction_ = 0.0;
        //angle_ = angle_correction_ = 0.0;


        cx_ = cy_ = cx_old_ = cy_old_ = 0.0;
        cx_unfolded_ = cy_unfolded_ = cx_unfolded_old_ = cy_unfolded_old_ = 0.0;
        vx_ = vy_ = vx_old_ = vy_old_ = 0.0;

        windtree_angle_ = 0.0;

        cx_initial_ = new double[n_trajectories_]();
        cy_initial_ = new double[n_trajectories_]();
        vx_initial_ = new double[n_trajectories_]();
        vy_initial_ = new double[n_trajectories_]();

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Linear Trajectory---" << std::endl;
            std::cout << "angle " << windtree_angle_ << std::endl;
            std::cout << "lifetime " << lifetime_ << std::endl;
        }

    }

    Linear ( int nr_trajectories, double lifetime, double windtree_angle )
    {
        //kahanFile.open ( "kahan.txt");
        //kahanFile.precision ( 18 );

        n_trajectories_ = nr_trajectories;
        lifetime_ = lifetime;
        current_time_ = 0.0;
        last_time_ = 0.0;
        time_correction_ = 0.0;

        cx_ = cy_ = cx_old_ = cy_old_ = 0.0;
        cx_unfolded_ = cy_unfolded_ = cx_unfolded_old_ = cy_unfolded_old_ = 0.0;
        vx_ = vy_ = vx_old_ = vy_old_ = 0.0;
        //angle_ = angle_correction_ = 0.0;

        windtree_angle_ = windtree_angle;

        cx_initial_ = new double[n_trajectories_]();
        cy_initial_ = new double[n_trajectories_]();
        vx_initial_ = new double[n_trajectories_]();
        vy_initial_ = new double[n_trajectories_]();
        //angle_inital_ = new double[n_trajectories_]();

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Linear Trajectory---" << std::endl;
            std::cout << "angle " << windtree_angle_ << std::endl;
            std::cout << "lifetime " << lifetime_ << std::endl;
        }
    }

    ~Linear()
    {
        delete[] cx_initial_;
        delete[] cy_initial_;
        delete[] vx_initial_;
        delete[] vy_initial_;
        //delete[] angle_inital_;
    }



    void reverse ( double time, Container &container )
    {
        get_current_position_inbox ( time, cx_, cy_ );
        container.applyBoundary ( cx_, cy_ );

        cx_old_ = cx_unfolded_ = cx_unfolded_old_ = cx_;
        cy_old_ = cy_unfolded_ = cy_unfolded_old_ = cy_;

        vx_ =  -vx_old_;
        vy_ =  -vy_old_;

        vx_old_ = vx_;
        vy_old_ = vy_;

        current_time_ = last_time_ = time_correction_ = 0.0;
    }


    inline void get_current_position_and_velocity ( double time, double &cx_t, double &cy_t, double &vx_t, double &vy_t )
    {
        double k = ( time - last_time_ );

        cx_t = cx_unfolded_old_ + k * vx_old_;
        cy_t = cy_unfolded_old_ + k * vy_old_;

        vx_t = vx_old_;
        vy_t = vy_old_;
    }

    inline void get_current_position ( double time, double &cx_t, double &cy_t )
    {
        double k = ( time - last_time_ );

        cx_t = cx_unfolded_old_ + k * vx_old_;
        cy_t = cy_unfolded_old_ + k * vy_old_;
    }

    inline void get_current_position_and_velocity_inbox ( double time, double &cx_t, double &cy_t, double &vx_t, double &vy_t )
    {
        double k = ( time - last_time_ );

        cx_t = cx_old_ + k * vx_old_;
        cy_t = cy_old_ + k * vy_old_;

        vx_t = vx_old_;
        vy_t = vy_old_;
    }

    inline void get_current_position_inbox ( double time, double &cx_t, double &cy_t )
    {
        double k = ( time - last_time_ );

        cx_t = cx_old_ + k * vx_old_;
        cy_t = cy_old_ + k * vy_old_;
    }

    inline void kahan_time_addition ( double k )
    {
        last_time_ = current_time_;
        double y = k - time_correction_;
        current_time_ = last_time_ + y;
        time_correction_ = ( current_time_ - last_time_ ) - y;
        // Here time_correction hast the negative value of the error:
        // if (current_time_ - last_time_) > y ----> the result is bigger than the real result
        //
    }

private:

    void engage ( int n )
    {
        cx_ = cx_old_ = cx_unfolded_ = cx_unfolded_old_ = cx_initial_[n];
        cy_ = cy_old_ = cy_unfolded_ = cy_unfolded_old_ = cy_initial_[n];

        vx_ =  vx_initial_[n];
        vy_ =  vy_initial_[n];

        vx_old_ = vx_initial_[n];
        vy_old_ = vy_initial_[n];

        current_time_ = last_time_ = time_correction_ = 0.0;
    }

};
};

