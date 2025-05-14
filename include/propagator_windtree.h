#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

#include "simulationbox.h"
#include "trajectory.h"
#include "neighbourlist.h"
#include "globals.h"

class WindtreePropagator
{
    friend class Testing;

private:

    double cell_l_;
    int ncells_;

    double *collision_buffer_;
    int *collided_obstacles_;

    Container * container_;
    Obstacles::Squares * obstacles_;
    Trajectory::Linear *trajectory_;
    Neighbourlist::Linear *neighlist_;


    int cell_nx_r;
    int cell_ny_r;

    int last_collided_item_;
    int nr_collisions_;

public:

    int nr_collision_points_;

public:

    WindtreePropagator()
    {
        collision_buffer_ = NULL;
        collided_obstacles_ = NULL;
    }

    ~WindtreePropagator()
    {
        delete[] collision_buffer_;
        delete[] collided_obstacles_;
        //delete[] collision_points_;
    }

    WindtreePropagator ( SimulationBox::PeriodicBoxSquareObstacles &simbox, Trajectory::Linear &trajectory, Neighbourlist::Linear &neighlist )
    {
        ncells_ = 0;
        container_= simbox.container_;
        obstacles_ = simbox.obstacles_;
        trajectory_ = &trajectory;
        neighlist_ = &neighlist;

        //caution as this propagator assumes angle = 0.0
        if ( trajectory_->windtree_angle_ != 0 )
        {
            std::cout << "Error Trajectory angle =! 0 in windtree Propagator" << std::endl;
            std::exit ( -1 );
        }
        trajectory_->windtree_angle_ = 0.0;

        auto len = neighlist.nlist_length_;

        nr_collision_points_ = 4 * len;

        collision_buffer_ = new double[nr_collision_points_];
        collided_obstacles_ = new int[nr_collision_points_];
        //collision_points_ = new double[len];

        last_collided_item_ = -10;
        double mem_size = sizeof ( double ) * nr_collision_points_ * 2.0 ;

        cell_l_ = neighlist.celllist_->l_;

        cell_nx_r = neighlist.celllist_->nx_ - 1;
        cell_ny_r = neighlist.celllist_->ny_ - 1;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Propagator(windtree)---" << std::endl;
            printf ( "total heapsize of the Propagator is %.2f kB \n", mem_size / 1000.0 );
        }
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

    void setStartPositions()
    {
        std::random_device rd;
        std::mt19937_64 engine ( rd() );
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );
        std::uniform_int_distribution<int> rand_int ( 0, 3 );

        double x, y;
        int index_end;

        //initializing random positions in the box and setting the radius
        for ( auto i = 0; i < trajectory_->n_trajectories_; ++i )
        {
            do
            {
                x = uni_dist ( engine ) * container_->lx_;
                y = uni_dist ( engine ) * container_->ly_;

                index_end = neighlist_->getPossibleCoveringNeighbours ( x, y );
            }
            while ( obstacles_->isInside ( neighlist_->neighbour_obstacles_, index_end, x, y ) );

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
        }
    }

    //Trajectory = Line : also returns the cell_s where the intersection auccurs.
    inline double reEnter ( int &cell_x, int &cell_y )
    {
        int side = -1;
        double k = 0.0;

        if ( trajectory_->vx_ == 1.0 )
        {
            side = 1;
            k = container_->lx_ - trajectory_->cx_;
        }

        if ( trajectory_->vx_ == -1.0 )
        {
            side = 3;
            k = trajectory_->cx_;
        }

        if ( trajectory_->vy_ == 1.0 )
        {
            side = 2;
            k = container_->ly_ - trajectory_->cy_;
        }

        if ( trajectory_->vy_ == -1.0 )
        {
            side = 0;
            k = trajectory_->cy_;
        }

        switch ( side )
        {
        case 0 :
        {
            cell_y = cell_ny_r;
            trajectory_->cy_ += container_->ly_;

            return k;
        }
        case 1 :
        {
            cell_x = 0;
            trajectory_->cx_ -= container_->lx_;
            return k;
        }
        case 2 :
        {
            cell_y = 0;
            trajectory_->cy_ -= container_->ly_;
            return k;
        }
        case 3 :
        {
            cell_x = cell_nx_r;
            trajectory_->cx_ += container_->lx_;
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
        trajectory_->kahan_time_addition ( k_reenter );
    }

    void engage ( int n )
    {
        trajectory_->engage ( n );
        last_collided_item_ = -10;
    }

    void inline propagate ( double end_time, Obstacles::Squares * dummy )
    {
        if ( trajectory_->current_time_ >= end_time )
        {
            return ;
        }
        
        int reenter = 0;
        double k_reenter = 0.0;

        int cell_x1 = 0;
        int cell_y1 = 0;

        /*will run until next collision. If not, then care schould ne taken if
         *      the propagator aborts and the particle is still outside of the container.
         */
        while ( trajectory_->current_time_ <= end_time )
        {
            double klo1 ( std::numeric_limits<double>::infinity() );
            double klo2 ( std::numeric_limits<double>::infinity() );

            int iklo1 = -1000;
            int iklo2 = -1001;

            if ( reenter == 0 )
            {
                neighlist_->celllist_->getCell ( trajectory_->cx_, trajectory_->cy_, cell_x1, cell_y1 );
            }

            ncells_ = neighlist_->getNeighbourCellsWindtree ( cell_x1, cell_y1, trajectory_->vx_, trajectory_->vy_ );

            nr_collisions_ = getIntersections_();

            int task = 0;

            if ( nr_collisions_ > 0 )
            {
                getTwoSmallest ( collision_buffer_, nr_collisions_, klo1, klo2, iklo1, iklo2 );

                if ( last_collided_item_ == collided_obstacles_[iklo1]/4)
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
                        if ( isNotValidIntersection_ ( collided_obstacles_[iklo1] ) )
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
                    placeReenter ( k_reenter );
                    k_reenter = 0.0;
                    reenter = 0;
                }
                break;
            }
            case 1:   // mirror on first intersection
            {
                last_collided_item_ = mirrorWindTree_ ( collided_obstacles_[iklo1], klo1 );
                reenter = 0;
                k_reenter = 0.0;
                break;
            }
            case 2:   // mirror on scond intersection
            {
                last_collided_item_ = mirrorWindTree_ ( collided_obstacles_[iklo2],klo2 );
                reenter = 0;
                k_reenter = 0.0;
                break;
            }
            }
        }
        return;
    };

    //Traj: line; Obstacle: squares
    inline int getIntersections_()
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

            for ( auto i = 0; i < index_last ; i++ )
            {
                int nr_obst = neighlist_->neighbour_obstacles_[i];
                int n_inter = 0;

                double p1x, p1y, p2x, p2y;
                double Ax, Ay, Cx, Cy;
                double den, num1, num2, k;
                int index = 5 * nr_obst;

                for ( int nr_seg = 0; nr_seg < 4; nr_seg++ )
                {
                    p1x = obstacles_->px_[index];
                    p1y = obstacles_->py_[index];
                    index++;
                    p2x = obstacles_->px_[index];
                    p2y = obstacles_->py_[index];

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
                        collision_buffer_[n] = k; // just saving the stretch of the unit vector ve
                        collided_obstacles_[n] = 4 * nr_obst + nr_seg;
                        n++;

                        if ( found == false and n == 2 )
                        {
                            found = true;
                            if ( ncells_ > curr_cell + 9 )
                            {
                                ncells_ = curr_cell + 9;
                            }
                        }

                        if ( ++n_inter > 1 )
                        {
                            break;
                        }
                    }
                }
            }//loop over obstacles
        }//loop over cells
        return n;
    }

//Traj: square; Obstacle: Line
    inline bool isNotValidIntersection_ ( int n_encoded )
    {
        //segment normal that is directed outside of the shape
        double nx =  obstacles_->ve_y_[n_encoded];
        double ny = -obstacles_->ve_x_[n_encoded];

        if ( trajectory_->vx_ * nx + trajectory_->vy_ * ny < 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

// Obstacles rotated by Pi/4 ccw

    inline int mirrorWindTree_ ( int n_encoded, double k )
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

        double temp;

        switch ( nr_seg )
        {
        case 0:
        {
            temp = trajectory_->vx_;
            trajectory_->vx_ = trajectory_->vy_;
            trajectory_->vy_ = temp;
            break;
        }
        case 1:
        {
            temp = trajectory_->vx_;
            trajectory_->vx_ = -trajectory_->vy_;
            trajectory_->vy_ = -temp;
            break;
        }
        case 2:
        {
            temp = trajectory_->vx_;
            trajectory_->vx_ = trajectory_->vy_;
            trajectory_->vy_ = temp;
            break;
        }
        case 3:
        {
            temp = trajectory_->vx_;
            trajectory_->vx_ = -trajectory_->vy_;
            trajectory_->vy_ = -temp;
            break;
        }
        }
        trajectory_->kahan_time_addition ( k );

        return nr_obst;
    };
};
