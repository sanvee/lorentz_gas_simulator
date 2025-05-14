#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "simulationbox.h"
#include "trajectory.h"
#include "neighbourlist.h"
#include "globals.h"
#include "io.h"

class LinearPropagator
{
    friend class Testing;

private:

    double cell_l_;
    int ncells_;
    int nobtacles_;

    double *collision_buffer_;
    int *collided_obstacles_;

    Trajectory::Linear *trajectory_;
    Neighbourlist::Linear *neighlist_;

    int cell_nx_r_;
    int cell_ny_r_;

    int last_collided_item_;
    int nr_collisions_;

public:
    Container *container_;
    int nr_collision_points_;

    LinearPropagator()
    {
        collision_buffer_ = NULL;
        collided_obstacles_ = NULL;
    }

    ~LinearPropagator()
    {
        delete[] collision_buffer_;
        delete[] collided_obstacles_;
        //delete[] collision_points_;
    }

    //***********//
    //**Squares**//
    //***********//
    LinearPropagator ( SimulationBox::PeriodicBoxSquareObstacles &simbox, Trajectory::Linear &trajectory, Neighbourlist::Linear &neighlist )
    {
        container_ = simbox.container_;
        trajectory_ = &trajectory;
        neighlist_ = &neighlist;

        ncells_ = 0;
        nobtacles_=simbox.obstacles_->n_tot_;

        nr_collision_points_ = 4 * neighlist.nlist_length_;

        collision_buffer_ = new double[nr_collision_points_];
        collided_obstacles_ = new int[nr_collision_points_];

        last_collided_item_ = -10;
        double mem_size = sizeof ( double ) * nr_collision_points_ * 2.0 ;
        //GlobalCounters::memsize += mem_size;

        cell_l_ = neighlist.celllist_->l_;

        cell_nx_r_ = neighlist.celllist_->nx_ - 1;
        cell_ny_r_ = neighlist.celllist_->ny_ - 1;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Propagator(linear)---" << std::endl;
            printf ( "total heapsize of the Propagator is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    //***********//
    //**Crosses**//
    //***********//
    LinearPropagator ( SimulationBox::PeriodicBoxCrossObstacles &simbox, Trajectory::Linear &trajectory, Neighbourlist::Linear &neighlist )
    {
        container_ = simbox.container_;
        trajectory_ = &trajectory;
        neighlist_ = &neighlist;

        ncells_ = 0;

        nr_collision_points_ = 8 * neighlist.nlist_length_;

        collision_buffer_ = new double[nr_collision_points_];
        collided_obstacles_ = new int[nr_collision_points_];

        last_collided_item_ = -10;
        double mem_size = sizeof ( double ) * nr_collision_points_ * 2.0 ;
        //GlobalCounters::memsize += mem_size;

        cell_l_ = neighlist.celllist_->l_;



        cell_nx_r_ = neighlist.celllist_->nx_ - 1;
        cell_ny_r_ = neighlist.celllist_->ny_ - 1;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Propagator(linear)---" << std::endl;
            printf ( "total heapsize of the Propagator is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    //***********//
    //**Circles**//
    //***********//

    LinearPropagator ( SimulationBox::PeriodicBoxCircularObstacles &simbox, Trajectory::Linear &trajectory, Neighbourlist::Linear &neighlist )
    {
        container_ = simbox.container_;
        trajectory_ = &trajectory;
        neighlist_ = &neighlist;

        ncells_ = 0;

        nr_collision_points_ = 4 * neighlist.nlist_length_;

        collision_buffer_ = new double[nr_collision_points_];
        collided_obstacles_ = new int[nr_collision_points_];

        last_collided_item_ = -10;
        double mem_size = sizeof ( double ) * nr_collision_points_ * 2.0 ;
        //GlobalCounters::memsize += mem_size;

        cell_l_ = neighlist.celllist_->l_;

        cell_nx_r_ = neighlist.celllist_->nx_ - 1;
        cell_ny_r_ = neighlist.celllist_->ny_ - 1;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Propagator(linear)---" << std::endl;
            printf ( "total heapsize of the Propagator is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    void engage ( int n )
    {
        trajectory_->engage ( n );
        last_collided_item_ = -10;
    }

private:

    inline void getTwoSmallest ( double *__restrict array, int lenght, double &low, double &sec_low, int &index_low, int &index_sec_low )
    {
        if ( lenght < 2 )
        {
            low = array[0];
            index_low = 0;

            sec_low = array[0];
            index_sec_low = 0;
            return;
        }
        //get the two smallest:
        low = array[0];
        index_low = 0;

        sec_low = array[1];
        index_sec_low = 1;

        if ( sec_low < low )
        {
            low = array[1];
            index_low = 1;

            sec_low = array[0];
            index_sec_low = 0;
        }
        for ( auto i = 2; i < lenght;  i++ )
        {
            if ( array[i] < low )
            {
                sec_low = low;
                index_sec_low = index_low;

                low = array[i];
                index_low = i;
            }
            else
            {
                if ( array[i] < sec_low )
                {
                    sec_low = array[i];
                    index_sec_low = i;
                }
            }
        }
    }

public:
    /*
        inline void setStartPositions()
        {
          //initializing random number generator
          gsl_rng *my_rng;
          gsl_rng_env_setup();
          std::ifstream urandom_FILE;
          urandom_FILE.open("/dev/urandom", std::ios::in | std::ios::binary);
          if (urandom_FILE.is_open())
          {
            urandom_FILE.read(reinterpret_cast < char *>(&gsl_rng_default_seed), sizeof(gsl_rng_default_seed));
          }
          urandom_FILE.close();
          my_rng = gsl_rng_alloc(gsl_rng_taus);

          double x, y;
          int index_end;

          //initializing random positions in the box and setting the radius
          for (auto i = 0; i < trajectory_->n_trajectories_; ++i)
          {
            do
            {
              x = gsl_rng_uniform(my_rng) * container_->lx_;
              y = gsl_rng_uniform(my_rng) * container_->ly_;

              index_end = neighlist_->getPossibleCoveringNeighbours(x, y);
            }
            while (obstacles_->isInside(neighlist_->neighbour_obstacles_, index_end, x, y));

            trajectory_->cx_initial_[i] = x;
            trajectory_->cy_initial_[i] = y;
          }

          for (auto i = 0; i < trajectory_->n_trajectories_; ++i)
          {
            double theta = gsl_rng_uniform(my_rng) * 2 * M_PI;
            trajectory_->vx_initial_[i] = cos(theta);
            trajectory_->vy_initial_[i] = sin(theta);
          }
          gsl_rng_free(my_rng);
        }

        */

    template<class T_Obstacles>
    inline void setStartPositions ( T_Obstacles *obstacles )
    {
        std::random_device rd;
        std::mt19937_64 engine ( rd() );
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );
        std::uniform_int_distribution<int> rand_int ( 0, 3 );

        double x, y;
        int index_end;
        double theta;

        //initializing random positions in the box and setting the radius
        for ( auto i = 0; i < trajectory_->n_trajectories_; ++i )
        {
            do
            {
                x = uni_dist ( engine ) * container_->lx_;
                y = uni_dist ( engine ) * container_->ly_;

                index_end = neighlist_->getPossibleCoveringNeighbours ( x, y );
            }
            while ( obstacles->isInside ( neighlist_->neighbour_obstacles_, index_end, x, y ) );

            trajectory_->cx_initial_[i] = x;
            trajectory_->cy_initial_[i] = y;
        }

        for ( auto i = 0; i < trajectory_->n_trajectories_; ++i )
        {
            int direction = rand_int ( engine );
            switch ( direction )
            {
            case 0:
            {
                trajectory_->vx_initial_[i] = 1.0;
                trajectory_->vy_initial_[i] = 0.0;
                break;
            }
            case 1:
            {
                trajectory_->vx_initial_[i] = 0.0;
                trajectory_->vy_initial_[i] = 1.0;
                break;
            }
            case 2:
            {
                trajectory_->vx_initial_[i] = -1.0;
                trajectory_->vy_initial_[i] =  0.0;
                break;
            }
            case 3:
            {
                trajectory_->vx_initial_[i] =  0.0;
                trajectory_->vy_initial_[i] = -1.0;
                break;
            }
            }

            //afterward rotation
            if ( trajectory_->windtree_angle_ != 0.0 )
            {
                if ( trajectory_->windtree_angle_ < 0.0 ) //random angle
                {
                    theta = uni_dist ( engine ) * M_PI_2;
                }
                else
                {
                    theta = trajectory_->windtree_angle_;
                }

                double c = std::cos ( theta );
                double s = std::sin ( theta );

                double vx_new = trajectory_->vx_initial_[i] * c - trajectory_->vy_initial_[i] * s;
                double vy_new = trajectory_->vx_initial_[i] * s + trajectory_->vy_initial_[i] * c;

                //trajectory_->angle_inital_[i] = theta + M_PI_2 * direction;

                trajectory_->vx_initial_[i] = vx_new;
                trajectory_->vy_initial_[i] = vy_new;
            }
        }
    }

    //Trajectory = Line : also returns the cell_s where the intersection auccurs.
    inline double reEnter ( int &cell_x, int &cell_y )
    {
        int n = 0;
        for ( int nr_seg = 0; nr_seg < 4; nr_seg++ )
        {
            double p1x = container_->segments->x_[nr_seg];
            double p1y = container_->segments->y_[nr_seg];

            double Ax = container_->segments->x_[nr_seg + 1] - p1x;
            double Ay = container_->segments->y_[nr_seg + 1] - p1y;

            double Bx = - trajectory_->vx_;
            double By = - trajectory_->vy_;

            double Cx = p1x - trajectory_->cx_;
            double Cy = p1y - trajectory_->cy_;

            double den = Ay * Bx - Ax * By;
            double num1 = By * Cx - Bx * Cy;
            double num2 = Ax * Cy - Ay * Cx;

            if ( den > 0.0 )
            {
                if ( num1 <= 0.0 or num1 > den )
                {
                    continue;
                }
                else
                {
                    if ( num2 <= 0.0 )
                    {
                        continue;
                    }
                }
            }
            else
            {
                if ( num1 > 0.0 or num1 < den )
                {
                    continue;
                }
                else
                {
                    if ( num2 > 0 )
                    {
                        continue;
                    }
                }
            }
            collision_buffer_[n] = num2 / den;
            collided_obstacles_[n] = nr_seg;
            n++;
        }

        int side;
        double k;
        //std::cout << n << std::endl;
        //std::cout << trajectory_->cx_ << ' ' << trajectory_->cy_ << std::endl;
        //if ( n == 0 )
        //{
        //  std::cout << "-----------------------------------Should not happen !!!!! " << std::endl;
        //}
        if ( n == 1 )   //sortie
        {
            k = collision_buffer_[0];
            side = collided_obstacles_[0];
        }
        else     //entree (prender le deuxieme plus grand)
        {
            if ( collision_buffer_[0] > collision_buffer_[1] )
            {
                side = collided_obstacles_[0];
                k = collision_buffer_[0];
            }
            else
            {
                side = collided_obstacles_[1];
                k = collision_buffer_[1];
            }
        }
        //std::cout << "n " << n << std::endl;
        //std::cout << "k " << k << std::endl;

        switch ( side )
        {
        case 0 :
        {
            neighlist_->celllist_->getCellx ( trajectory_->cx_ + trajectory_->vx_ * k, cell_x );
            cell_y = cell_ny_r_;
            trajectory_->cy_ += container_->ly_;
            //std::cout << "added ly_" << std::endl;
            //std::cout << "cell_x " << cell_x << std::endl;
            //std::cout << "cell_y " << cell_y << std::endl;
            return k;
        }
        case 1 :
        {
            cell_x = 0;
            neighlist_->celllist_->getCelly ( trajectory_->cy_ + trajectory_->vy_ * k, cell_y );
            trajectory_->cx_ -= container_->lx_;
            //std::cout << "sub lx_" << std::endl;
            //std::cout << "cell_x " << cell_x << std::endl;
            //std::cout << "cell_y " << cell_y << std::endl;
            return k;
        }
        case 2 :
        {
            neighlist_->celllist_->getCellx ( trajectory_->cx_ + trajectory_->vx_ * k, cell_x );
            cell_y = 0;
            trajectory_->cy_ -= container_->ly_;
            //std::cout << "sub ly_" << std::endl;
            //std::cout << "cell_x " << cell_x << std::endl;
            //std::cout << "cell_y " << cell_y << std::endl;
            return k;
        }
        case 3 :
        {
            cell_x = cell_nx_r_;
            neighlist_->celllist_->getCelly ( trajectory_->cy_ + trajectory_->vy_ * k, cell_y );
            trajectory_->cx_ += container_->lx_;
            //std::cout << "added lx_" << std::endl;
            //std::cout << "cell_x " << cell_x << std::endl;
            //std::cout << "cell_y " << cell_y << std::endl;
            return k;
        }
        }

        std::cout << "something Wrong in Container intersections !!!" << std::endl;
        return 0.0;
    }

    void inline placeReenter ( double k_reenter )
    {
        double dx = trajectory_->vx_ * k_reenter;
        double dy = trajectory_->vy_ * k_reenter;

        trajectory_->cx_old_ = trajectory_->cx_;
        trajectory_->cy_old_ = trajectory_->cy_;

        trajectory_->cx_ += dx;
        trajectory_->cy_ += dy;

        container_->remapInside ( trajectory_->cx_, trajectory_->cy_ );

        trajectory_->cx_unfolded_old_ = trajectory_->cx_unfolded_;
        trajectory_->cy_unfolded_old_ = trajectory_->cy_unfolded_;

        trajectory_->cx_unfolded_ += dx;
        trajectory_->cy_unfolded_ += dy;

        trajectory_->vx_old_ = trajectory_->vx_;
        trajectory_->vy_old_ = trajectory_->vy_;

        trajectory_->last_time_ = trajectory_->current_time_;
        trajectory_->current_time_ += k_reenter; // again we suppose that v = 1;

    }

    template<class T_Obstacles>
    inline void propagate ( double end_time, T_Obstacles *obstacles )
    {
        if ( trajectory_->current_time_ >= end_time )
        {
            return;
        }

        int reenter = 0;
        double k_reenter = 0.0;
        /*will run until next collision. If not, then care schould ne taken if
        the propagator aborts and the particle is still outside of the container_->
        */
        int cell_x1 = 0;
        int cell_y1 = 0;

        while ( trajectory_->current_time_ <= end_time )
        {
            double klo1;
            double klo2;

            klo1 = klo2 = std::numeric_limits<double>::infinity();

            int iklo1 = -2;
            int iklo2 = -3;

            if ( reenter == 0 )
            {
                neighlist_->celllist_->getCell ( trajectory_->cx_, trajectory_->cy_, cell_x1, cell_y1 );
            }

            ncells_ = neighlist_->getNeighbourCells ( cell_x1, cell_y1, trajectory_->vx_, trajectory_->vy_ );
            nr_collisions_ = getIntersections ( obstacles );

            int task = 0;

            if ( nr_collisions_ > 0 )
            {
                getTwoSmallest ( collision_buffer_, nr_collisions_, klo1, klo2, iklo1, iklo2 );

                if ( last_collided_item_ == collided_obstacles_[iklo1]/obstacles->sides_ )
                {
                    if ( nr_collisions_ == 1 )
                    {
                        task = 0;
                    }
                    else
                    {
                        task = 2;
                    }
                }
                else
                {
                    if ( klo1 > 10e-4 )
                    {
                        task = 1;
                    }
                    else
                    {
                        if ( isNotValidIntersection ( *obstacles, collided_obstacles_[iklo1], trajectory_->cx_ + trajectory_->vx_ * klo1, trajectory_->cy_ + trajectory_->vy_ * klo1 ) )
                        {
                            if ( nr_collisions_ == 1 )
                            {
                                task = 0;
                            }
                            else
                            {
                                task = 2;
                            }
                        }
                        else
                        {
                            task = 1;
                        }
                    }
                }
            }
            else
            {
                task = 0;
            }
            switch ( task )
            {
            case 0:   // reenter
            {
                k_reenter = reEnter ( cell_x1, cell_y1 );
                reenter ++;
                last_collided_item_ = -10;
                if ( reenter > 2 )
                {
                    //std::cout << "Place Reenter" << std::endl;
                    placeReenter ( k_reenter );
                    //need to take possible outside of the container ...
                    k_reenter = 0.0;
                    reenter = 0;
                }

                break;
            }
            case 1:   // mirror on first intersection
            {
                last_collided_item_ = mirror ( *obstacles, collided_obstacles_[iklo1], klo1 );
                reenter = 0;
                k_reenter = 0.0;
                break;
            }
            case 2:   // mirror on second intersection
            {
                last_collided_item_ = mirror ( *obstacles, collided_obstacles_[iklo2], klo2 );
                reenter = 0;
                k_reenter = 0.0;
                break;
            }
            }
        }
        return;
    };

//---------------Collisions----------------
//-----------------------------------------


//-----Traj: line; Obstacle: Crosses
//----------------------------------

    inline int getIntersections ( Obstacles::Crosses *obstacles )
    {
        int n = 0;
        double Bx = -trajectory_->vx_;
        double By = -trajectory_->vy_;
        double cx = trajectory_->cx_;
        double cy = trajectory_->cy_;
        bool found = false;

        for ( int curr_cell = 0; curr_cell < ncells_; curr_cell++ )
        {
            int index_last = neighlist_->getNeighboursInStencilCell ( curr_cell );

            for ( auto i = 0; i < index_last; i++ )

            {
                int nr_obst = neighlist_->neighbour_obstacles_[i];
                int n_inter = 0;

                double p1x, p1y, p2x, p2y;
                double Ax, Ay, Cx, Cy;
                double den, num1, num2, k;
                int index = 13 * nr_obst;

                for ( int nr_seg = 0; nr_seg < 12; nr_seg++ )
                {
                    p1x = obstacles->px_[index];
                    p1y = obstacles->py_[index];
                    index++;
                    p2x = obstacles->px_[index];
                    p2y = obstacles->py_[index];

                    Ax = p2x - p1x;
                    Ay = p2y - p1y;

                    Cx = p1x - cx;
                    Cy = p1y - cy;

                    den = Ay * Bx - Ax * By;
                    num1 = By * Cx - Bx * Cy;

                    if ( den > 0.0 )
                    {
                        if ( num1 <= 0.0 or num1 > den )
                        {
                            continue;
                        }
                        else
                        {
                            num2 = Ax * Cy - Ay * Cx;
                            if ( num2 <= 0.0 )
                            {
                                continue;
                            }
                        }
                    }
                    else
                    {
                        if ( num1 > 0.0 or num1 < den )
                        {
                            continue;
                        }
                        else
                        {
                            num2 = Ax * Cy - Ay * Cx;
                            if ( num2 > 0 )
                            {
                                continue;
                            }
                        }
                    }
                    k = num2 / den;

                    if ( container_->isInside ( cx - k * Bx, cy - k * By ) )
                    {
                        //std::cout << n << std::endl;
                        collision_buffer_[n] = k; // just saving the stretch of the unit vector ve
                        collided_obstacles_[n] = 12 * nr_obst + nr_seg;
                        n++;

                        if ( found == false and n == 3 )
                        {
                            found = true;
                            if ( ncells_ > curr_cell + 15 )
                            {
                                ncells_ = curr_cell + 15;
                            }
                        }
                        if ( ++n_inter > 4 )
                        {
                            break;
                        }
                        //only max 4 collisions are possible with a cross and a line
                    }
                }
            }
        }
        return n;
    }


//-----Traj: line; Obstacle: squares
//----------------------------------

    inline int getIntersections ( Obstacles::Squares *obstacles )
    {
        int n = 0;
        double Bx = -trajectory_->vx_;
        double By = -trajectory_->vy_;
        double cx = trajectory_->cx_;
        double cy = trajectory_->cy_;
        bool found = false;

        for ( int curr_cell = 0; curr_cell < ncells_; curr_cell++ )
        {
            int index_last = neighlist_->getNeighboursInStencilCell ( curr_cell );

            for ( auto i = 0; i < index_last; i++ )
            {
                int nr_obst = neighlist_->neighbour_obstacles_[i];
                int n_inter = 0;

                double p1x, p1y, p2x, p2y;
                double Ax, Ay, Cx, Cy;
                double den, num1, num2, k;
                int index = 5 * nr_obst;

                for ( int nr_seg = 0; nr_seg < 4; nr_seg++ )
                {
                    p1x = obstacles->px_[index];
                    p1y = obstacles->py_[index];
                    index++;
                    p2x = obstacles->px_[index];
                    p2y = obstacles->py_[index];

                    Ax = p2x - p1x;
                    Ay = p2y - p1y;

                    Cx = p1x - cx;
                    Cy = p1y - cy;

                    den = Ay * Bx - Ax * By;
                    num1 = By * Cx - Bx * Cy;

                    if ( den > 0.0 )
                    {
                        if ( num1 <= 0.0 or num1 > den )
                        {
                            continue;
                        }
                        else
                        {
                            num2 = Ax * Cy - Ay * Cx;
                            if ( num2 <= 0.0 )
                            {
                                continue;
                            }
                        }
                    }
                    else
                    {
                        if ( num1 > 0.0 or num1 < den )
                        {
                            continue;
                        }
                        else
                        {
                            num2 = Ax * Cy - Ay * Cx;
                            if ( num2 > 0 )
                            {
                                continue;
                            }
                        }
                    }
                    k = num2 / den;

                    if ( container_->isInside ( cx - k * Bx, cy - k * By ) )
                    {
                        //std::cout << n << std::endl;
                        collision_buffer_[n] = k; // just saving the stretch of the unit vector ve
                        collided_obstacles_[n] = 4 * nr_obst + nr_seg;
                        n++;

                        if ( found == false and n == 3 )
                        {
                            found = true;
                            if ( ncells_ > curr_cell + 15 )
                            {
                                ncells_ = curr_cell + 15;
                            }
                        }
                        if ( ++n_inter > 1 )
                        {
                            break;
                        }
                        //only max 2 collisions are possible as we have a convex shape
                    }
                }
            }
        }
        return n;
    }

//----Traj: Line; Obstacle: Circle
//--------------------------------

    inline int getIntersections ( Obstacles::Circles *obstacles )
    {
        auto evx = trajectory_->vx_;
        auto evy = trajectory_->vy_;

        auto cx = trajectory_->cx_;
        auto cy = trajectory_->cy_;

        auto etx =  evy;
        auto ety = -evx;

        int n = 0;
        bool found = false;

        for ( int curr_cell = 0; curr_cell < ncells_; curr_cell++ )
        {
            int n_neight = neighlist_->getNeighboursInStencilCell ( curr_cell );

            for ( auto i = 0; i < n_neight; i++ )
            {
                int n_obs = neighlist_->neighbour_obstacles_[i];

                double r_obs = obstacles->radius_[n_obs];

                double obs_x = obstacles->x_[n_obs];
                double obs_y = obstacles->y_[n_obs];

                double h = ( trajectory_->cx_ - obs_x ) * etx + ( trajectory_->cy_ - obs_y ) * ety;

                if ( h >= r_obs or - h >= r_obs )
                {
                    continue;
                }

                double b = std::sqrt ( ( r_obs + h ) * ( r_obs - h ) );
                //instead of:
                //double b = std::sqrt(r_obs * r_obs - h * h);

                double dx = obs_x + etx * h + evx * b;
                double dy = obs_y + ety * h + evy * b;
                double skp = ( dx - cx ) * evx + ( dy - cy ) * evy;

                if ( container_->isInside ( dx, dy ) and skp > 0 )
                {
                    collision_buffer_[n] = skp; // just saving the stretch of the unit vector ve
                    collided_obstacles_[n] = n_obs;
                    n++;

                    if ( found == false and n == 2 )
                    {
                        found = true;
                        if ( ncells_ > curr_cell + 15 )
                        {
                            ncells_ = curr_cell + 15;
                        }
                    }
                }

                dx -= 2 * evx * b;
                dy -= 2 * evy * b;
                skp = ( dx - cx ) * evx + ( dy - cy ) * evy;

                if ( container_->isInside ( dx, dy ) and skp > 0 )
                {
                    collision_buffer_[n] = skp; // just saving the stretch of the unit vector ve
                    collided_obstacles_[n] = n_obs;
                    n++;

                    if ( found == false and n == 2 )
                    {
                        found = true;
                        if ( ncells_ > curr_cell + 15 )
                        {
                            ncells_ = curr_cell + 15;
                        }
                    }
                }
            }
        }
        return n;
    }


    //----Mirroring
    //-------------

    //Traj: Line; Obstacle: Circle
    inline int mirror ( Obstacles::Circles &obstacles, int n_obs, double k )
    {
        double dx = trajectory_->vx_ * k;
        double dy = trajectory_->vy_ * k;

        trajectory_->cx_old_ = trajectory_->cx_;
        trajectory_->cy_old_ = trajectory_->cy_;

        trajectory_->cx_ += dx;
        trajectory_->cy_ += dy;

        trajectory_->cx_unfolded_old_ = trajectory_->cx_unfolded_;
        trajectory_->cy_unfolded_old_ = trajectory_->cy_unfolded_;

        trajectory_->cx_unfolded_ += dx;
        trajectory_->cy_unfolded_ += dy;

        trajectory_->vx_old_ = trajectory_->vx_;
        trajectory_->vy_old_ = trajectory_->vy_;

        //Vector perpendicular to the reflection axis
        double r_obs_inv = 1.0 / obstacles.radius_[n_obs];
        
        double nox = ( trajectory_->cy_ - obstacles.y_[n_obs] ) * r_obs_inv;
        double noy = ( obstacles.x_[n_obs] - trajectory_->cx_ ) * r_obs_inv;
        
        double H11 = ( 2.0 * nox * nox - 1.0 );
        double H12 = 2 * nox * noy;

        trajectory_->vx_ = ( H11 * trajectory_->vx_old_ + H12 * trajectory_->vy_old_ );
        trajectory_->vy_ = ( H12 * trajectory_->vx_old_ - H11 * trajectory_->vy_old_ );

        //normalisation----------                
        H12 = 1.0/ (std::sqrt( trajectory_->vx_ * trajectory_->vx_ + trajectory_->vy_ * trajectory_->vy_) );
        
        trajectory_->vx_ *= H12;
        trajectory_->vy_ *= H12;
        
        //-----------------------
        trajectory_->kahan_time_addition ( k );
        return n_obs;
    };

    //Traj: Line; Obstacle: square
    inline int mirror ( Obstacles::Squares &obstacles, int n_encoded, double k )
    {
        int nr_obst = n_encoded /4;

        double dx = trajectory_->vx_ * k;
        double dy = trajectory_->vy_ * k;

        trajectory_->cx_old_ = trajectory_->cx_;
        trajectory_->cy_old_ = trajectory_->cy_;

        trajectory_->cx_ += dx;
        trajectory_->cy_ += dy;

        trajectory_->cx_unfolded_old_ = trajectory_->cx_unfolded_;
        trajectory_->cy_unfolded_old_ = trajectory_->cy_unfolded_;

        trajectory_->cx_unfolded_ += dx;
        trajectory_->cy_unfolded_ += dy;

        trajectory_->vx_old_ = trajectory_->vx_;
        trajectory_->vy_old_ = trajectory_->vy_;

        //Vector perpendicular to the reflection axis

        double nox, noy, norm, Hdiag, H12;

        nox = obstacles.ve_x_[n_encoded];
        noy = obstacles.ve_y_[n_encoded];

        Hdiag = ( nox - noy ) * ( nox + noy );

        H12 = 2.0 * nox * noy;

        trajectory_->vx_ = Hdiag * trajectory_->vx_old_ + H12 * trajectory_->vy_old_;
        trajectory_->vy_ = H12 * trajectory_->vx_old_ - Hdiag * trajectory_->vy_old_;

        //normalisation----------
        H12 = 1.0/ (std::sqrt( trajectory_->vx_ * trajectory_->vx_ + trajectory_->vy_ * trajectory_->vy_));
        trajectory_->vx_ *= H12;
        trajectory_->vy_ *= H12;
        //-----------------------

        trajectory_->kahan_time_addition ( k );

        return nr_obst;
    };

    //Traj: Line; Obstacle: crosses
    inline int mirror ( Obstacles::Crosses &obstacles, int n_encoded, double k )
    {
        int nr_obst = n_encoded / 12;

        double dx = trajectory_->vx_ * k;
        double dy = trajectory_->vy_ * k;

        trajectory_->cx_old_ = trajectory_->cx_;
        trajectory_->cy_old_ = trajectory_->cy_;

        trajectory_->cx_ += dx;
        trajectory_->cy_ += dy;

        trajectory_->cx_unfolded_old_ = trajectory_->cx_unfolded_;
        trajectory_->cy_unfolded_old_ = trajectory_->cy_unfolded_;

        trajectory_->cx_unfolded_ += dx;
        trajectory_->cy_unfolded_ += dy;

        trajectory_->vx_old_ = trajectory_->vx_;
        trajectory_->vy_old_ = trajectory_->vy_;

        //Vector perpendicular to the reflection axis

        double nox, noy, Hdiag, H12;

        nox = obstacles.ve_x_[n_encoded];
        noy = obstacles.ve_y_[n_encoded];

        Hdiag = ( nox - noy ) * ( nox + noy );

        H12 = 2.0 * nox * noy;

        trajectory_->vx_ = Hdiag * trajectory_->vx_old_ + H12 * trajectory_->vy_old_;
        trajectory_->vy_ = H12 * trajectory_->vx_old_ - Hdiag * trajectory_->vy_old_;

        //normalisation----------
        H12 = 1.0/ (std::sqrt( trajectory_->vx_ * trajectory_->vx_ + trajectory_->vy_ * trajectory_->vy_) );
        trajectory_->vx_ *= H12;
        trajectory_->vy_ *= H12;
        //-----------------------

        trajectory_->kahan_time_addition ( k );

        return nr_obst;
    };


    //Traj: Line; Obstacle: square
    inline int mirror_0 ( Obstacles::Squares &obstacles, int n_encoded, double k )
    {
        int nr_obst = n_encoded / 4;
        int nr_seg = n_encoded - 4 * nr_obst;

        double dx = trajectory_->vx_ * k;
        double dy = trajectory_->vy_ * k;

        trajectory_->cx_old_ = trajectory_->cx_;
        trajectory_->cy_old_ = trajectory_->cy_;

        trajectory_->cx_ += dx;
        trajectory_->cy_ += dy;

        trajectory_->cx_unfolded_old_ = trajectory_->cx_unfolded_;
        trajectory_->cy_unfolded_old_ = trajectory_->cy_unfolded_;

        trajectory_->cx_unfolded_ += dx;
        trajectory_->cy_unfolded_ += dy;

        trajectory_->vx_old_ = trajectory_->vx_;
        trajectory_->vy_old_ = trajectory_->vy_;

        //depending on the side we just reverse one component
        if ( nr_seg % 2==0 )
        {
            trajectory_->vy_ = -trajectory_->vy_;
        }
        else
        {
            trajectory_->vx_ = -trajectory_->vx_;
        }
        trajectory_->kahan_time_addition ( k );

        return nr_obst;
    };

    //Traj: Line; Obstacle: square
    inline int mirror_0 ( Obstacles::Crosses &obstacles, int n_encoded, double k )
    {
        int nr_obst = n_encoded / 12;
        int nr_seg = n_encoded - 12 * nr_obst;

        double dx = trajectory_->vx_ * k;
        double dy = trajectory_->vy_ * k;

        trajectory_->cx_old_ = trajectory_->cx_;
        trajectory_->cy_old_ = trajectory_->cy_;

        trajectory_->cx_ += dx;
        trajectory_->cy_ += dy;

        trajectory_->cx_unfolded_old_ = trajectory_->cx_unfolded_;
        trajectory_->cy_unfolded_old_ = trajectory_->cy_unfolded_;

        trajectory_->cx_unfolded_ += dx;
        trajectory_->cy_unfolded_ += dy;

        trajectory_->vx_old_ = trajectory_->vx_;
        trajectory_->vy_old_ = trajectory_->vy_;

        //depending on the side we just reverse one component
        if ( nr_seg % 2==0 )
        {
            trajectory_->vx_ = -trajectory_->vx_;
        }
        else
        {
            trajectory_->vy_ = -trajectory_->vy_;
        }
        trajectory_->kahan_time_addition ( k );

        return nr_obst;
    };



    //----validity check
    //------------------

    //---Traj: line; Obstacle: Circle
    //-------------------------------
    inline bool isNotValidIntersection ( Obstacles::Circles &obstacles, int n_obs, double intersection_x, double intersection_y )
    {
        double ox = obstacles.x_[n_obs] - intersection_x;
        double oy = obstacles.y_[n_obs] - intersection_y;

        if ( trajectory_->vx_ * ox + trajectory_->vy_ * oy > 0.0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    //----Traj: line; Obstacle: Squares
    //--------------------------------
    inline bool isNotValidIntersection ( Obstacles::Squares &obstacles, int n_encoded, double intersection_x, double intersection_y )
    {
        //auto nr_obst = n_encoded / 4;
        //int edge = n_encoded % 4;

        //segment normal that is directed outside of the shape
        //double nx =  obstacles.ve_y_[4 * nr_obst + edge];
        //double ny = -obstacles.ve_x_[4 * nr_obst + edge];

        double nx =  obstacles.ve_y_[n_encoded];
        double ny = -obstacles.ve_x_[n_encoded];

        if ( trajectory_->vx_ * nx + trajectory_->vy_ * ny < 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    //----Traj: line; Obstacle: Crosses
    //--------------------------------
    inline bool isNotValidIntersection ( Obstacles::Crosses &obstacles, int n_encoded, double intersection_x, double intersection_y )
    {
        //auto nr_obst = n_encoded / 12;
        //int edge = n_encoded % 12;

        //segment normal that is directed outside of the shape
        double nx =  obstacles.ve_y_[n_encoded];
        double ny = -obstacles.ve_x_[n_encoded];

        if ( trajectory_->vx_ * nx + trajectory_->vy_ * ny < 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }
};



