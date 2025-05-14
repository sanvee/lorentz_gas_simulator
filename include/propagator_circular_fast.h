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

class CircularPropagator
{
    friend class Testing;
private:

    //std::ofstream debug_file_;

    double *collision_buffer_;
    int *collided_obstacles_;
    int ncells_;
    double *collision_angles_;

    double *cont_inters_;
    double *cont_inters_angles_;
    int *cont_inters_seg_;

    Trajectory::Circular *trajectory_;
    Neighbourlist::Circular *neighlist_;

    bool filter_circling_;


public:
    Container *container_;
    //int nr_collision_points_;

    CircularPropagator()
    {
        collision_buffer_ = NULL;
        collided_obstacles_ = NULL;
        collision_angles_ = NULL;

        cont_inters_ = NULL;
        cont_inters_angles_ = NULL;
        cont_inters_seg_ = NULL;
    }

    ~CircularPropagator()
    {
        delete[] collision_buffer_;
        delete[] collided_obstacles_;
        delete[] collision_angles_;

        delete[] cont_inters_;
        delete[] cont_inters_angles_;
        delete[] cont_inters_seg_;
    }

    template<class Simbox>
    CircularPropagator ( Simbox &simbox, Trajectory::Circular &trajectory, Neighbourlist::Circular &neighlist, bool filter )
    {
        filter_circling_= filter;
        //debug_file_.open ( "debug_inters.txt" );
        container_ = simbox.container_;
        trajectory_ = &trajectory;
        neighlist_ = &neighlist;

        /// really a maximum could be decreased ...
        int nr_collision_points =  2 * simbox.obstacles_->sides_ * neighlist.nlist_length_;

        /// for consecutive x and y coordinate
        collision_buffer_ = new double[2 * nr_collision_points];
        collided_obstacles_ = new int[2 * nr_collision_points];
        collision_angles_ = new double[nr_collision_points];

        cont_inters_ = new double[18];
        cont_inters_angles_ = new double[8];
        cont_inters_seg_ = new int [8];

        double mem_size = sizeof ( double ) * nr_collision_points * 5.0 ;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Propagator(circular directonal version)---" << std::endl;
            printf ( "total heapsize of the Propagator is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    template <class T_Obstacles>
    inline void setStartPositions ( T_Obstacles *obstacles )
    {
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        double x, y, theta;
        int index_end;

        //initializing random positions in the box
        for ( auto i = 0; i < trajectory_->n_trajectories_; ++i )
        {
            do
            {
                x = uni_dist ( Globalrand::e2 ) * container_->lx_;
                y = uni_dist ( Globalrand::e2 ) * container_->ly_;

                index_end = neighlist_->getPossibleCoveringNeighbours ( x, y );
            }
            while ( obstacles->isInside ( neighlist_->neighbour_obstacles_, index_end, x, y ) );

            trajectory_->cx_initial_[i] = x;
            trajectory_->cy_initial_[i] = y;

            // sampling the circle origin (velocity direction)
            theta = uni_dist ( Globalrand::e2 ) * 2 * M_PI;
            trajectory_->rx_initial_[i] = cos ( theta ) * trajectory_->radius_ + x;
            trajectory_->ry_initial_[i] = sin ( theta ) * trajectory_->radius_ + y;

            if ( filter_circling_ )
            {
                // now check if this trajectory is colliding with something
                trajectory_->engage ( i );

                if ( !checkIntersections ( *obstacles ) )
                {
                    //do not increment as the trajectory is not valid
                    i--;
                }
            }
        }
    }

    template <class T_Obstacles>
    inline void setStartPositions_no_check ( T_Obstacles *obstacles )
    {

        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        double x, y, theta;
        int index_end;

        //initializing random positions in the box
        for ( auto i = 0; i < trajectory_->n_trajectories_; ++i )
        {
            do
            {
                x = uni_dist ( Globalrand::e2 ) * container_->lx_;
                y = uni_dist ( Globalrand::e2 ) * container_->ly_;

                index_end = neighlist_->getPossibleCoveringNeighbours ( x, y );
            }
            while ( obstacles->isInside ( neighlist_->neighbour_obstacles_, index_end, x, y ) );

            trajectory_->cx_initial_[i] = x;
            trajectory_->cy_initial_[i] = y;

            // sampling the circle origin (velocity direction)
            theta = uni_dist ( Globalrand::e2 ) * 2 * M_PI;
            trajectory_->rx_initial_[i] = cos ( theta ) * trajectory_->radius_ + x;
            trajectory_->ry_initial_[i] = sin ( theta ) * trajectory_->radius_ + y;
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
    template<class T_Obstacles>
    void propagate ( double end_time, T_Obstacles &obstacles )
    {
        //Assumption: Every trajectory has at least one collision.
        // do nothing and go out if
        /**********************/

        if ( trajectory_->current_time_ >= end_time )
        {
            return;
        }

        double theta_low;
        double theta_low2;
        double theta_low_c;
        double theta_low2_c;

        double theta_ratchet;
        double theta_compare;

        double this_cx;
        double this_cy;

        int nr_collisions;
        int i_low = 0;
        int i_low2 = 0;
        int i_low_c = 0;
        int i_low2_c = 0;
        int nr_trepassing;
        int nr_reenter;

        nr_reenter = 0;
        theta_ratchet = 0.0;

        while ( trajectory_->current_time_ < end_time )
        {
            nr_trepassing = 0;

            //first we need to know if the trajectory leaves the container
            if ( canLeaveContainer ( *container_, *trajectory_ ) )
            {
                nr_trepassing = getContainerIntersections ();
                //std::cout << "reenter: " << nr_reenter << std::endl;

                for ( auto i = 0; i < nr_trepassing ;  ++i )
                {
                    cont_inters_angles_[i] = atan2 ( cont_inters_[2 * i + 1] - trajectory_->ry_, cont_inters_[2 * i] - trajectory_->rx_ ) - trajectory_->theta_begin_;
                    trajectory_->normalizeAngle ( cont_inters_angles_[i] );

                    if ( cont_inters_angles_[i] < theta_ratchet )
                    {
                        cont_inters_angles_[i] += Constants::TWO_PI;
                    }
                }

                getTwoSmallest ( cont_inters_angles_, nr_trepassing, theta_low_c, theta_low2_c, i_low_c, i_low2_c );

                //now compute theta ratchet the comparison angle
                if ( nr_reenter == 0 )
                {
                    //std::cout << "case 1" << std::endl;
                    theta_ratchet = theta_low_c * 0.5;
                    theta_compare = theta_low_c;
                    ncells_ = neighlist_->getNeighbourCells ( trajectory_->cx_begin_, trajectory_->cy_begin_, cont_inters_[2*i_low_c], cont_inters_[2*i_low_c + 1], trajectory_->rx_, trajectory_->ry_ );
                }
                else
                {
                    //std::cout << "case 2" << std::endl;
                    theta_ratchet = ( theta_low_c + theta_low2_c ) * 0.5;
                    theta_compare = theta_low2_c;
                    ncells_ = neighlist_->getNeighbourCells ( cont_inters_[2*i_low_c], cont_inters_[2*i_low_c+1],cont_inters_[2*i_low2_c], cont_inters_[2*i_low2_c+1],  trajectory_->rx_, trajectory_->ry_ );
                }
            }
            else
            {
                //std::cout << "case 3" << std::endl;
                theta_ratchet = Constants::TWO_PI;
                theta_compare = Constants::TWO_PI;
                ncells_ = neighlist_->getNeighbourCells ( trajectory_->cx_begin_, trajectory_->cy_begin_, trajectory_->rx_, trajectory_->ry_ );
            }

            nr_collisions = getIntersections ( obstacles );

            if ( nr_collisions == 0 )
            {
                if ( nr_reenter > 0 ) // select second edge
                {
                    reEnter ( cont_inters_seg_[i_low2_c] );
                }
                else
                {
                    reEnter ( cont_inters_seg_[i_low_c] );
                }
                ++nr_reenter;
            }
            else
            {
                for ( auto i = 0; i < nr_collisions; ++i )
                {
                    collision_angles_[i] = atan2 ( collision_buffer_[2 * i + 1] - trajectory_->ry_, collision_buffer_[2 * i] - trajectory_->rx_ ) - trajectory_->theta_begin_;
                    trajectory_->normalizeAngle ( collision_angles_[i] );
                }

                if ( nr_collisions == 1 ) // Dummy
                {
                    collision_angles_[1] = Constants::TWO_PI;
                    collision_buffer_[2] = 0.0;
                    collision_buffer_[3] = 0.0;
                    collided_obstacles_[1] = 0;
                    nr_collisions = 2;
                }

                getTwoSmallest ( collision_angles_, nr_collisions, theta_low, theta_low2, i_low, i_low2 );

                //check the validity of the collision
                this_cx = collision_buffer_[i_low * 2];
                this_cy = collision_buffer_[i_low * 2 + 1];

                double dr = ( trajectory_->cx_begin_ - this_cx ) * ( trajectory_->cx_begin_ - this_cx ) + ( trajectory_->cy_begin_ - this_cy ) * ( trajectory_->cy_begin_ - this_cy );

                if ( dr < 1.0e-4 )
                {
                    if ( isNotValidIntersection ( obstacles, collided_obstacles_[i_low], this_cx, this_cy ) )
                    {
                        this_cx = collision_buffer_[i_low2 * 2];
                        this_cy = collision_buffer_[i_low2 * 2 + 1];
                        i_low = i_low2;
                        theta_low = theta_low2;
                    }
                }

                if ( theta_low < theta_compare ) //mirror
                {
                    mirror ( obstacles, collided_obstacles_[i_low], this_cx, this_cy, theta_low );
                    trajectory_->n_tot_coll_ ++;
                    theta_ratchet = 0.0;
                    nr_reenter = 0;
                }
                else //shift
                {
                    if ( nr_reenter > 0 ) // select second edge
                    {
                        reEnter ( cont_inters_seg_[i_low2_c] );
                    }
                    else
                    {
                        reEnter ( cont_inters_seg_[i_low_c] );
                    }
                    ++nr_reenter;
                }
            }
        }
        return;
    };

    template<class T_Obstacles>
    bool test_coll ( T_Obstacles &obstacles )
    {
        //Assumption: Every trajectory has at least one collision.
        // do nothing and go out if
        /**********************/
        double theta_low_c;
        double theta_low2_c;
        double theta_ratchet;
        double theta_compare;
        int nr_collisions;
        int i_low_c = 0;
        int i_low2_c = 0;
        int nr_trepassing;
        int nr_reenter;

        nr_reenter = 0;
        theta_ratchet = 0.0;

        while ( true )
        {
            nr_trepassing = 0;

            //first we need to know if the trajectory leaves the container
            if ( canLeaveContainer ( *container_, *trajectory_ ) )
            {
                nr_trepassing = getContainerIntersections ();
                //std::cout << "reenter: " << nr_reenter << std::endl;

                for ( auto i = 0; i < nr_trepassing ;  ++i )
                {
                    cont_inters_angles_[i] = atan2 ( cont_inters_[2 * i + 1] - trajectory_->ry_, cont_inters_[2 * i] - trajectory_->rx_ ) - trajectory_->theta_begin_;
                    trajectory_->normalizeAngle ( cont_inters_angles_[i] );

                    if ( cont_inters_angles_[i] < theta_ratchet )
                    {
                        cont_inters_angles_[i] += Constants::TWO_PI;
                    }
                }

                getTwoSmallest ( cont_inters_angles_, nr_trepassing, theta_low_c, theta_low2_c, i_low_c, i_low2_c );

                //now compute theta ratchet the comparison angle
                if ( nr_reenter == 0 )
                {
                    //std::cout << "case 1" << std::endl;
                    theta_ratchet = theta_low_c * 0.5;
                    theta_compare = theta_low_c;
                    ncells_ = neighlist_->getNeighbourCells ( trajectory_->cx_begin_, trajectory_->cy_begin_, cont_inters_[2*i_low_c], cont_inters_[2*i_low_c + 1], trajectory_->rx_, trajectory_->ry_ );
                }
                else
                {
                    //std::cout << "case 2" << std::endl;
                    theta_ratchet = ( theta_low_c + theta_low2_c ) * 0.5;
                    theta_compare = theta_low2_c;
                    ncells_ = neighlist_->getNeighbourCells ( cont_inters_[2*i_low_c], cont_inters_[2*i_low_c+1],cont_inters_[2*i_low2_c], cont_inters_[2*i_low2_c+1],  trajectory_->rx_, trajectory_->ry_ );
                }
            }
            else
            {
                //std::cout << "case 3" << std::endl;
                theta_ratchet = Constants::TWO_PI;
                theta_compare = Constants::TWO_PI;
                ncells_ = neighlist_->getNeighbourCells ( trajectory_->cx_begin_, trajectory_->cy_begin_, trajectory_->rx_, trajectory_->ry_ );
            }

            nr_collisions = getIntersections ( obstacles );

            if ( nr_collisions == 0 )
            {
                if ( nr_trepassing == 0 )
                {
                    return false;
                }

                if ( theta_ratchet >= 2*M_PI )
                {
                    return false;
                }

                if ( nr_reenter > 0 ) // select second edge
                {
                    reEnter ( cont_inters_seg_[i_low2_c] );
                }
                else
                {
                    reEnter ( cont_inters_seg_[i_low_c] );
                }
                ++nr_reenter;
            }
            else
            {
                return true;
            }
        }
        return false;
    };


    //helper functions:

    //---------------------------Container Intersections-----------------------------//


    //checks if the trajectory can leave the container
    inline bool canLeaveContainer ( Container &container, Trajectory::Circular &trajectory )
    {
        if ( trajectory.rx_ < trajectory.radius_ )
        {
            return true;
        }
        if ( trajectory.rx_ > container.lx_ - trajectory.radius_ )
        {
            return true;
        }

        if ( trajectory.ry_ < trajectory.radius_ )
        {
            return true;
        }
        if ( trajectory.ry_ > container.ly_ - trajectory.radius_ )
        {
            return true;
        }

        return false;
    }

    //calculates the container intersections
    inline int getContainerIntersections ()
    {
        double x_tr = trajectory_->rx_;
        double y_tr = trajectory_->ry_;
        double r_tr = trajectory_->radius_;

        int n = 0;
        for ( auto nr_seg = 0; nr_seg < 4; ++nr_seg )
        {
            double x1_obs = container_->segments->x_[nr_seg];
            double y1_obs = container_->segments->y_[nr_seg];

            double x2_obs = container_->segments->x_[nr_seg + 1];
            double y2_obs = container_->segments->y_[nr_seg + 1];

            double l_obs_inv = container_->segments->length_inv_[nr_seg];

            double ehx = ( y2_obs - y1_obs ) * l_obs_inv;
            double ehy = ( x1_obs - x2_obs ) * l_obs_inv;

            double dx = x1_obs - x_tr;
            double dy = y1_obs - y_tr;

            double altitude = dx * ehx + dy * ehy;

            if ( std::abs ( altitude ) >= r_tr )
            {
                continue;
            }

            double base = std::sqrt ( ( r_tr - altitude ) * ( r_tr + altitude ) );

            //first point
            dx = x_tr + altitude * ehx - base * ehy;
            dy = y_tr + altitude * ehy + base * ehx;

            //Check if the intersection points are between the points of the segment
            if ( ( dx - x1_obs ) * ( dx - x2_obs ) + ( dy - y1_obs ) * ( dy - y2_obs ) < 0 )
            {
                cont_inters_[2 * n]   = dx;
                cont_inters_[2 * n + 1] = dy;
                cont_inters_seg_[n] = nr_seg;
                n++;
            }

            dx = x_tr + altitude * ehx + base * ehy;
            dy = y_tr + altitude * ehy - base * ehx;

            if ( ( dx - x1_obs ) * ( dx - x2_obs ) + ( dy - y1_obs ) * ( dy - y2_obs ) < 0 )
            {
                cont_inters_[2 * n]   = dx;
                cont_inters_[2 * n + 1] = dy;
                cont_inters_seg_[n] = nr_seg;
                n++;
            }
        }
        return n;
    }

    //---------------------------Obstacles Intersections-----------------------------//

    //checks if the trajectory intersects the obstacles (Circular Obstacles)
    bool checkIntersections ( Obstacles::Circles &obstacles )
    {
        double r_tr = trajectory_->radius_;

        double r_obs;
        double dx, dy;
        double d_sq;

        for ( auto i = 0; i < obstacles.n_tot_; i++ )
        {
            r_obs = obstacles.radius_[i];

            dx = trajectory_->rx_ - obstacles.x_[i];
            dy = trajectory_->ry_ - obstacles.y_[i];

            //compute the squared distance between the two circles
            d_sq = dx * dx + dy * dy;

            //check if the two circles are closer than the sum of their radius:
            if ( d_sq >= ( r_obs + r_tr ) * ( r_obs + r_tr ) or d_sq <= ( r_obs - r_tr ) * ( r_obs - r_tr ) )
            {
                continue;
            }
            else
            {
                return true;
            }
        }
        return false;
    }

    //checks if the trajectory intersects the obstacles (Square Obstacles)
    bool checkIntersections ( Obstacles::Squares &obstacles )
    {
        double x_tr = trajectory_->rx_;
        double y_tr = trajectory_->ry_;
        double r_tr = trajectory_->radius_;

        for ( auto i = 0; i < obstacles.n_tot_ ; i++ )
        {
            double esx, esy;

            double altitude, base;
            double x0_obs, y0_obs;
            double x1_obs, y1_obs;

            double dx, dy;

            for ( int nr_seg = 0; nr_seg < 4; ++nr_seg )
            {

                esx = obstacles.ve_x_[4 * i + nr_seg];
                esy = obstacles.ve_y_[4 * i + nr_seg];

                //ehx = esy;
                //ehy = -esx;

                x0_obs = obstacles.px_[5 * i + nr_seg];
                y0_obs = obstacles.py_[5 * i + nr_seg];

                x1_obs = obstacles.px_[5 * i + nr_seg + 1];
                y1_obs = obstacles.py_[5 * i + nr_seg + 1];

                dx = x1_obs - x_tr;
                dy = y1_obs - y_tr;

                altitude = dx * esy - dy * esx;

                //base = std::sqrt(r_tr_sq - altitude * altitude);
                base = std::sqrt ( ( r_tr + altitude ) * ( r_tr - altitude ) );

                //first point
                dx = x_tr + altitude * esy + base * esx;
                dy = y_tr - altitude * esx + base * esy;

                //Check if the intersection points P are between the points of the segment A and B:
                //take scalar product of the PA.PB it should be < 0;

                if ( ( dx - x0_obs ) * ( dx - x1_obs ) + ( dy - y0_obs ) * ( dy - y1_obs ) < 0 )
                {
                    return true;
                }

                dx = x_tr + altitude * esy - base * esx;
                dy = y_tr - altitude * esx - base * esy;

                if ( ( dx - x0_obs ) * ( dx - x1_obs ) + ( dy - y0_obs ) * ( dy - y1_obs ) < 0 )
                {
                    return true;
                }
            }
        }
        return false;
    }

    //checks if the trajectory intersects the obstacles (Crosses Obstacles)
    bool checkIntersections ( Obstacles::Crosses &obstacles )
    {
        double x_tr = trajectory_->rx_;
        double y_tr = trajectory_->ry_;
        double r_tr = trajectory_->radius_;

        for ( auto i = 0; i < obstacles.n_tot_ ; i++ )
        {
            double esx, esy;

            double altitude, base;
            double x0_obs, y0_obs;
            double x1_obs, y1_obs;

            double dx, dy;


            for ( int nr_seg = 0; nr_seg < 12; ++nr_seg )
            {

                esx = obstacles.ve_x_[12 * i + nr_seg];
                esy = obstacles.ve_y_[12 * i + nr_seg];

                //ehx = esy;
                //ehy = -esx;

                x0_obs = obstacles.px_[13 * i + nr_seg];
                y0_obs = obstacles.py_[13 * i + nr_seg];

                x1_obs = obstacles.px_[13 * i + nr_seg + 1];
                y1_obs = obstacles.py_[13 * i + nr_seg + 1];

                dx = x1_obs - x_tr;
                dy = y1_obs - y_tr;

                altitude = dx * esy - dy * esx;

                //base = std::sqrt(r_tr_sq - altitude * altitude);
                base = std::sqrt ( ( r_tr + altitude ) * ( r_tr - altitude ) );

                //first point
                dx = x_tr + altitude * esy + base * esx;
                dy = y_tr - altitude * esx + base * esy;

                //Check if the intersection points P are between the points of the segment A and B:
                //take scalar product of the PA.PB it should be < 0;

                if ( ( dx - x0_obs ) * ( dx - x1_obs ) + ( dy - y0_obs ) * ( dy - y1_obs ) < 0 )
                {
                    return true;
                }

                dx = x_tr + altitude * esy - base * esx;
                dy = y_tr - altitude * esx - base * esy;

                if ( ( dx - x0_obs ) * ( dx - x1_obs ) + ( dy - y0_obs ) * ( dy - y1_obs ) < 0 )
                {
                    return true;
                }
            }
        }
        return false;
    }

    //Collision calculations: There routines must only return collisions inside of the container.

    //Traj: Circle; Obstacle: Circle
    inline int getIntersections ( Obstacles::Circles &obstacles )
    {
        double x_tr = trajectory_->rx_;
        double y_tr = trajectory_->ry_;
        double r_tr = trajectory_->radius_;

        int n = 0;
        int n_pre = 0;
        bool found = false;

        int n_checkedcells = 0;

        for ( int curr_cell = 0; curr_cell < ncells_; curr_cell++ )
        {
            n_checkedcells++;
            int index_last = neighlist_->getNeighboursInStencilCell ( curr_cell );

            for ( auto i = 0; i < index_last; i++ )
            {
                int n_obs = neighlist_->neighbour_obstacles_[i];

                double x_obs = obstacles.x_[n_obs];
                double y_obs = obstacles.y_[n_obs];
                double r_obs = obstacles.radius_[n_obs];

                double T = r_tr + r_obs ;
                double t = r_obs - r_tr;

                double max_dist_sq = T*T;
                double min_dist_sq = t*t;

                double dx = x_tr - x_obs;
                double dy = y_tr - y_obs;

                //compute the squared distance between the two circles
                double d_sq = dx * dx + dy * dy;

                //check if the two circles are closer than the sum of their radius or the other lies inside of the other
                if ( d_sq >= max_dist_sq or d_sq <= min_dist_sq )
                {
                    continue;
                }

                double d = std::sqrt ( d_sq );
                double d_inv = 1.0 / d;

                //compute the distance from the center of the trajectory circle on the jonction line of the two centers
                double u = d_sq + t * T;
                u /= 2.0 * d;

                //compute the vertical offset from this line
                double v = std::sqrt ( ( r_obs - u ) * ( r_obs + u ) );

                //compute unit vector along, and perpendicular to the connecting line
                double e1x = d_inv * dx;
                double e1y = d_inv * dy;

                //compute the two intersection points
                dx = x_obs + u * e1x + v * -e1y;
                dy = y_obs + u * e1y + v * e1x;

                if ( container_->isInside ( dx, dy ) )
                {
                    collision_buffer_[2 * n]   = dx;
                    collision_buffer_[2 * n + 1] = dy;
                    collided_obstacles_[n] = n_obs;
                    n++;
                }

                dx -= 2 * v * -e1y;
                dy -= 2 * v * e1x;

                if ( container_->isInside ( dx, dy ) )
                {
                    collision_buffer_[2 * n]   = dx;
                    collision_buffer_[2 * n + 1] = dy;
                    collided_obstacles_[n] = n_obs;
                    n++;
                }
            }

            if ( n_checkedcells == 78 )
            {
                n_pre = n;
            }

            if ( n_checkedcells > 78 )
            {
                if ( found == false and ( n - n_pre ) > 3 )
                {
                    found = true;
                    if ( ncells_ > curr_cell + 50 )
                    {
                        ncells_ = curr_cell + 50;
                    }
                }
            }
        }
        return n;
    }

//Traj: Circle; Obstacle: Squares
    inline int getIntersections ( Obstacles::Squares &obstacles )
    {
        double x_tr = trajectory_->rx_;
        double y_tr = trajectory_->ry_;
        double r_tr = trajectory_->radius_;

        bool found = false;

        int n = 0;
        int n_pre = 0;
        int n_checkedcells = 0;

        for ( int curr_cell = 0; curr_cell < ncells_; curr_cell++ )
        {
            n_checkedcells++;
            int index_last = neighlist_->getNeighboursInStencilCell ( curr_cell );

            for ( auto i = 0; i < index_last ; i++ )
            {
                double esx, esy;

                double altitude, base;
                double x0_obs, y0_obs;
                double x1_obs, y1_obs;

                double dx, dy;

                auto nr_obst = neighlist_->neighbour_obstacles_[i];

                for ( auto nr_seg = 0; nr_seg < 4; ++nr_seg )
                {
                    esx = obstacles.ve_x_[4 * nr_obst + nr_seg];
                    esy = obstacles.ve_y_[4 * nr_obst + nr_seg];

                    //unit vector perpendicular to the segments
                    //ehx = esy;
                    //ehy = -esx;

                    x0_obs = obstacles.px_[5 * nr_obst + nr_seg];
                    y0_obs = obstacles.py_[5 * nr_obst + nr_seg];

                    x1_obs = obstacles.px_[5 * nr_obst + nr_seg + 1];
                    y1_obs = obstacles.py_[5 * nr_obst + nr_seg + 1];

                    dx = x1_obs - x_tr;
                    dy = y1_obs - y_tr;

                    altitude = dx * esy - dy * esx;

                    //base = std::sqrt(r_tr * r_tr - altitude * altitude);
                    base = std::sqrt ( ( r_tr - altitude ) * ( r_tr + altitude ) );

                    if ( base > r_tr )
                    {
                        continue;
                    }

                    //first point
                    dx = x_tr + altitude * esy + base * esx;
                    dy = y_tr - altitude * esx + base * esy;

                    //Check if the intersection points P are between the points of the segment A and B:
                    //take scalar product of the PA.PB it should be < 0;

                    if ( ( dx - x0_obs ) * ( dx - x1_obs ) + ( dy - y0_obs ) * ( dy - y1_obs ) < 0 and container_->isInside ( dx, dy ) )
                    {
                        collision_buffer_[2 * n    ] = dx;
                        collision_buffer_[2 * n + 1] = dy;
                        collided_obstacles_[n] = 4 * nr_obst + nr_seg;
                        n++;
                    }

                    dx = x_tr + altitude * esy - base * esx;
                    dy = y_tr - altitude * esx - base * esy;

                    if ( ( dx - x0_obs ) * ( dx - x1_obs ) + ( dy - y0_obs ) * ( dy - y1_obs ) < 0 and container_->isInside ( dx, dy ) )
                    {
                        collision_buffer_[2 * n    ] = dx;
                        collision_buffer_[2 * n + 1] = dy;
                        collided_obstacles_[n] = 4 * nr_obst + nr_seg;
                        n++;
                    }
                }
            }
            if ( n_checkedcells == 78 )
            {
                n_pre = n;
            }

            if ( n_checkedcells > 78 )
            {
                if ( found == false and ( n - n_pre ) > 3 )
                {
                    found = true;

                    if ( ncells_ > curr_cell + 50 )
                    {
                        ncells_ = curr_cell + 50;
                    }
                }
            }
        }
        return n;
    }

//Traj: Circle; Obstacle: Crosses
    inline int getIntersections ( Obstacles::Crosses &obstacles )
    {
        double x_tr = trajectory_->rx_;
        double y_tr = trajectory_->ry_;
        double r_tr = trajectory_->radius_;

        bool found = false;

        int n = 0;
        int n_pre = 0;
        int n_checkedcells = 0;

        for ( int curr_cell = 0; curr_cell < ncells_; curr_cell++ )
        {
            n_checkedcells ++;
            int index_last = neighlist_->getNeighboursInStencilCell ( curr_cell );
            for ( auto i = 0; i < index_last ; i++ )
            {
                double esx, esy;

                double altitude, base;
                double x0_obs, y0_obs;
                double x1_obs, y1_obs;

                double dx, dy;

                auto nr_obst = neighlist_->neighbour_obstacles_[i];

                for ( auto nr_seg = 0; nr_seg < 12; ++nr_seg )
                {

                    esx = obstacles.ve_x_[12 * nr_obst + nr_seg];
                    esy = obstacles.ve_y_[12 * nr_obst + nr_seg];

                    //unit vector perpendicular to the segments
                    //ehx = esy;
                    //ehy = -esx;

                    x0_obs = obstacles.px_[13 * nr_obst + nr_seg];
                    y0_obs = obstacles.py_[13 * nr_obst + nr_seg];

                    x1_obs = obstacles.px_[13 * nr_obst + nr_seg + 1];
                    y1_obs = obstacles.py_[13 * nr_obst + nr_seg + 1];

                    dx = x1_obs - x_tr;
                    dy = y1_obs - y_tr;

                    altitude = dx * esy - dy * esx;

                    //base = std::sqrt(r_tr * r_tr - altitude * altitude);
                    base = std::sqrt ( ( r_tr + altitude ) * ( r_tr - altitude ) );

                    if ( base > r_tr )
                    {
                        continue;
                    }

                    //first point
                    dx = x_tr + altitude * esy + base * esx;
                    dy = y_tr - altitude * esx + base * esy;

                    //Check if the intersection points P are between the points of the segment A and B:
                    //take scalar product of the PA.PB it should be < 0;

                    if ( ( dx - x0_obs ) * ( dx - x1_obs ) + ( dy - y0_obs ) * ( dy - y1_obs ) < 0 and container_->isInside ( dx, dy ) )
                    {
                        collision_buffer_[2 * n    ] = dx;
                        collision_buffer_[2 * n + 1] = dy;
                        collided_obstacles_[n] = 12 * nr_obst + nr_seg;
                        n++;
                    }

                    dx = x_tr + altitude * esy - base * esx;
                    dy = y_tr - altitude * esx - base * esy;

                    if ( ( dx - x0_obs ) * ( dx - x1_obs ) + ( dy - y0_obs ) * ( dy - y1_obs ) < 0 and container_->isInside ( dx, dy ) )
                    {
                        collision_buffer_[2 * n    ] = dx;
                        collision_buffer_[2 * n + 1] = dy;
                        collided_obstacles_[n] = 12 * nr_obst + nr_seg;
                        n++;
                    }
                }
            }

            if ( n_checkedcells == 78 )
            {
                n_pre = n;
            }

            if ( n_checkedcells > 78 )
            {
                if ( found == false and ( n - n_pre ) > 3 )
                {
                    found = true;
                    if ( ncells_ > curr_cell + 50 )
                    {
                        ncells_ = curr_cell + 50;
                    }
                }
            }
        }
        return n;
    }

//Traj: Circle;
    inline void reEnter ( int side )
    {
        switch ( side )
        {
        case 0 :
        {
            trajectory_->ry_ += container_->ly_;
            trajectory_->cy_begin_   += container_->ly_;
            break;
        }
        case 1 :
        {
            trajectory_->rx_ -= container_->lx_;
            trajectory_->cx_begin_   -= container_->lx_;
            break;
        }
        case 2 :
        {
            trajectory_->ry_ -= container_->ly_;
            trajectory_->cy_begin_   -= container_->ly_;
            break;
        }
        case 3 :
        {
            trajectory_->rx_ += container_->lx_;
            trajectory_->cx_begin_   += container_->lx_;
            break;
        }
        }
    }

//mirror
//Traj: Circle; Obstacle: Circle
    inline void mirror ( Obstacles::Circles &obstacles, int n_obs, double cx, double cy, double collision_angle )
    {
        double x_obs = obstacles.x_[n_obs];
        double y_obs = obstacles.y_[n_obs];

        double epx = ( cy - obstacles.y_[n_obs] );
        double epy = ( obstacles.x_[n_obs] - cx );

        //double norm = 1.0 /obstacles.radius_[n_obs];
        double norm = 1.0/std::sqrt ( epx*epx+epy*epy );
        epx *= norm;
        epy *= norm;

        double dx = x_obs - trajectory_->rx_;
        double dy = y_obs - trajectory_->ry_;

        double h = epx * dx + epy * dy;

        dx = 2 * h * epx;
        dy = 2 * h * epy;

        trajectory_->rx_old_ = trajectory_->rx_;
        trajectory_->ry_old_ = trajectory_->ry_;

        trajectory_->rx_unfolded_old_ = trajectory_->rx_unfolded_;
        trajectory_->ry_unfolded_old_ = trajectory_->ry_unfolded_;

        trajectory_->rx_ += dx;
        trajectory_->ry_ += dy;

        trajectory_->rx_unfolded_ += dx;
        trajectory_->ry_unfolded_ += dy;

        trajectory_->cx_begin_ = cx;
        trajectory_->cy_begin_ = cy;

        trajectory_->theta_begin_old_ = trajectory_->theta_begin_;
        trajectory_->theta_begin_ = atan2 ( cy - trajectory_->ry_, cx - trajectory_->rx_ );

        trajectory_->kahan_time_addition ( collision_angle * trajectory_->radius_ );
    };

//Traj: Circle; Obstacle: Square
    inline void mirror ( Obstacles::Squares &obstacles, int n_encoded, double cx, double cy, double collision_angle )
    {
        int nr_obst = n_encoded / 4;
        int nr_seg = n_encoded % 4;

        //unit vector along the segment
        double epx = obstacles.ve_x_[4 * nr_obst + nr_seg];
        double epy = obstacles.ve_y_[4 * nr_obst + nr_seg];

        double dx = cx - trajectory_->rx_;
        double dy = cy - trajectory_->ry_;

        double b = epx * dx + epy * dy;

        dx = 2 * b * epx;
        dy = 2 * b * epy;

        trajectory_->rx_old_ = trajectory_->rx_;
        trajectory_->ry_old_ = trajectory_->ry_;

        trajectory_->rx_unfolded_old_ = trajectory_->rx_unfolded_;
        trajectory_->ry_unfolded_old_ = trajectory_->ry_unfolded_;

        trajectory_->rx_ += dx;
        trajectory_->ry_ += dy;

        trajectory_->rx_unfolded_ += dx;
        trajectory_->ry_unfolded_ += dy;

        trajectory_->cx_begin_ = cx;
        trajectory_->cy_begin_ = cy;

        trajectory_->theta_begin_old_ = trajectory_->theta_begin_;
        trajectory_->theta_begin_ = atan2 ( cy - trajectory_->ry_, cx - trajectory_->rx_ );

        //double diff = -trajectory.current_time_;
        trajectory_->kahan_time_addition ( collision_angle * trajectory_->radius_ );

        //if (trajectory.current_time_  > 1e5)
        //{
        //  deb_file << std::setprecision(18) << trajectory .current_time_ + diff << ' ' << collision_angle * trajectory.radius_ <<std::endl;
        //}
    };

//Traj: Circle; Obstacle: Crosses
    inline void mirror ( Obstacles::Crosses &obstacles, int n_encoded, double cx, double cy, double collision_angle )
    {

        int nr_obst = n_encoded / 12;
        int nr_seg = n_encoded % 12;

        //unity vector along the segment
        double epx = obstacles.ve_x_[12 * nr_obst + nr_seg];
        double epy = obstacles.ve_y_[12 * nr_obst + nr_seg];

        double dx = cx - trajectory_->rx_;
        double dy = cy - trajectory_->ry_;

        double h = epx * dx + epy * dy;

        dx = 2 * h * epx;
        dy = 2 * h * epy;

        trajectory_->rx_old_ = trajectory_->rx_;
        trajectory_->ry_old_ = trajectory_->ry_;

        trajectory_->rx_unfolded_old_ = trajectory_->rx_unfolded_;
        trajectory_->ry_unfolded_old_ = trajectory_->ry_unfolded_;

        trajectory_->rx_ += dx;
        trajectory_->ry_ += dy;

        trajectory_->rx_unfolded_ += dx;
        trajectory_->ry_unfolded_ += dy;

        trajectory_->cx_begin_ = cx;
        trajectory_->cy_begin_ = cy;

        trajectory_->theta_begin_old_ = trajectory_->theta_begin_;
        trajectory_->theta_begin_ = atan2 ( cy - trajectory_->ry_, cx - trajectory_->rx_ );

        trajectory_->kahan_time_addition ( collision_angle * trajectory_->radius_ );
    };

//Traj: Circle; Obstacle: Circle
    inline bool isNotValidIntersection ( Obstacles::Circles &obstacles, int n_obs, double intersection_x, double intersection_y )
    {
        double vx = trajectory_->ry_ - intersection_y;
        double vy = intersection_x - trajectory_->rx_;

        double ox = obstacles.x_[n_obs] - intersection_x;
        double oy = obstacles.y_[n_obs] - intersection_y;

        if ( vx * ox + vy * oy > 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

//Traj: Circle; Obstacle: Square
    inline bool isNotValidIntersection ( Obstacles::Squares &obstacles, int n_encoded, double intersection_x, double intersection_y )
    {
        auto nr_obst = n_encoded / 4;
        int edge = n_encoded % 4;

        //velocity vector:
        double vx = trajectory_->ry_ - intersection_y;
        double vy = intersection_x - trajectory_->rx_;

        //segment normal that is directed outside of the shape
        double nx =  obstacles.ve_y_[4 * nr_obst + edge];
        double ny = -obstacles.ve_x_[4 * nr_obst + edge];

        if ( vx * nx + vy * ny < 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

//Traj: Circle; Obstacle: Cross
    inline bool isNotValidIntersection ( Obstacles::Crosses &obstacles, int n_encoded, double intersection_x, double intersection_y )
    {
        auto nr_obst = n_encoded / 12;
        int edge = n_encoded % 12;

        //velocity vector:
        double vx = trajectory_->ry_ - intersection_y;
        double vy = intersection_x - trajectory_->rx_;

        //segment normal that is directed outside of the shape
        double nx =  obstacles.ve_y_[12 * nr_obst + edge];
        double ny = -obstacles.ve_x_[12 * nr_obst + edge];

        if ( vx * nx + vy * ny < 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }
};


