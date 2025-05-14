#pragma once

#include <iostream>

#include "globals.h"

namespace Obstacles
{
class Crosses
{
public:

    std::vector<double> x_;    // center position
    std::vector<double> y_;

    std::vector<double> ve_x_; // unit vector parrallel to edge
    std::vector<double> ve_y_;

    std::vector<double> px_;   // vertices of the squares
    std::vector<double> py_;

    double l_; // Length of the arm from the center
    double d_; // Widtht of the arm

    int n_tot_;
    int n_obstacles_;
    int n_ghosts_;

    double angle_;
    std::vector<double> angles_;

    double max_radius_;
    static constexpr int sides_ = 12;
    static constexpr int type_ = 2;


    Crosses()
    {
        //li_ = NULL;
        l_ = 0;
        d_ = 0;

        n_tot_ = 0;
        n_obstacles_ = 0;
        n_ghosts_ = 0;

        angle_ = 0;
        max_radius_ = 0.0;

    }

    ///l = width of the square d is the width of the arm from edge to edge.
    ///l_= is half the width of the cross.(intern)
    Crosses ( double l, double d, double angle, int n_obstacles, int n_ghosts )
    {
        l_ =  0.5 * l;
        d_ = d;

        max_radius_ = std::sqrt ( l_ * l_ + 0.25 * d_ * d_ );

        n_obstacles_ = n_obstacles;

        if ( n_obstacles > 200000000 )
        {
            std::cout << "too many obsatcles ... ... ... aborting" << std::endl;
            std::terminate();
        }

        n_ghosts_ = n_ghosts;
        n_tot_ = n_ghosts + n_obstacles;

        x_.resize ( n_tot_ );
        y_.resize ( n_tot_ );

        px_.resize ( 13 * n_tot_ );
        py_.resize ( 13 * n_tot_ );

        ve_x_.resize ( 12 * n_tot_ );
        ve_y_.resize ( 12 * n_tot_ );

        //li_ = new double [12] {d_, l_, l_, d_, l_, l_, d_, l_, l_, d_, l_, l_};

        // if unique angle is used
        angle_ = angle;

        // if each obstacles has a random angle
        if(angle_ < 0)
        {
            angles_.resize ( n_tot_ );
        }

        double mem_size = sizeof ( double ) * 20 * n_tot_;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Obstacles (Crosses)---" << std::endl;
            std::cout << "l_cross " << 2.0 * l_ << std::endl;
            std::cout << "d_cross " << d_ << std::endl;
            std::cout << "n Obstacles " << n_obstacles_ << std::endl;
            std::cout << "n Ghosts " << n_ghosts_ << std::endl;
            printf ( "total heapsize of the Obstacles is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    bool isInside ( int *obstacles_list, int last, double x, double y )
    {
        for ( auto i = 0; i < last; i++ )
        {
            int n = obstacles_list[i];

            double ap_x = x - px_[13 * n];
            double ap_y = y - py_[13 * n];

            double sp1 = ( ve_x_[12 * n] * ap_x + ve_y_[12 * n] * ap_y ) * d_;
            double sp2 = - ( ve_x_[12 * n + 11] * ap_x + ve_y_[12 * n + 11] * ap_y ) * 2 * l_;

            ap_x = x - px_[13 * n + 3];
            ap_y = y - py_[13 * n + 3];

            double sp3 = ( ve_x_[12 * n + 3] * ap_x + ve_y_[12 * n + 3] * ap_y ) * d_;
            double sp4 = - ( ve_x_[12 * n + 10] * ap_x + ve_y_[12 * n + 10] * ap_y ) * l_ * 2;

            if ( ( sp1 > 0 and sp1 < d_ * d_ and sp2 > 0 and sp2 < 4 * l_ * l_ )
                    or ( sp3 > 0 and sp3 < d_ * d_ and sp4 > 0 and sp4 < 4 * l_ * l_ ) )
            {
                return true;
            }
        }
        return false;
    }

    bool isInside ( int n, double x, double y )
    {
        double ap_x = x - px_[13 * n];
        double ap_y = y - py_[13 * n];

        double sp1 = ( ve_x_[12 * n] * ap_x + ve_y_[12 * n] * ap_y ) * d_;
        double sp2 = - ( ve_x_[12 * n + 11] * ap_x + ve_y_[12 * n + 11] * ap_y ) * 2 * l_;

        ap_x = x - px_[13 * n + 3];
        ap_y = y - py_[13 * n + 3];

        double sp3 = ( ve_x_[12 * n + 3] * ap_x + ve_y_[12 * n + 3] * ap_y ) * d_;
        double sp4 = - ( ve_x_[12 * n + 10] * ap_x + ve_y_[12 * n + 10] * ap_y ) * l_ * 2;

        if ( ( sp1 > 0 and sp1 < d_ * d_ and sp2 > 0 and sp2 < 4 * l_ * l_ )
                or ( sp3 > 0 and sp3 < d_ * d_ and sp4 > 0 and sp4 < 4 * l_ * l_ ) )
        {
            return true;
        }

        return false;
    }

    inline void copy_i_to_j_ ( int i, int j )
    {
        x_[j] = x_[i];
        y_[j] = y_[i];

        for ( int m = 0; m < 13; m++ )
        {
            px_[13 * j + m] = px_[13 * i + m];
            py_[13 * j + m] = py_[13 * i + m];
        }

        for ( int m = 0; m < 13; m++ )
        {
            ve_x_[12 * j + m] = ve_x_[12 * i + m];
            ve_y_[12 * j + m] = ve_y_[12 * i + m];
        }

        if ( angle_ < 0 )
        {
            angles_[j] = angles_[i];
        }
    }

    void set ( int n, double alpha, double x, double y )
    {
        x_[n] = x;
        y_[n] = y;

        if(angle_ < 0)
        {
            angles_[n] = alpha;
        }
        double d = 0.5 * d_;

        if ( alpha == 0.0 )
        {
            px_[13 * n] = x_[n] + l_;
            py_[13 * n] = y_[n] - d;

            px_[13 * n + 1] = x_[n] + l_;
            py_[13 * n + 1] = y_[n] + d;

            px_[13 * n + 2] = x_[n] + d;
            py_[13 * n + 2] = y_[n] + d;

            px_[13 * n + 3] = x_[n] + d;
            py_[13 * n + 3] = y_[n] + l_;

            px_[13 * n + 4] = x_[n] - d;
            py_[13 * n + 4] = y_[n] + l_;

            px_[13 * n + 5] = x_[n] - d;
            py_[13 * n + 5] = y_[n] + d;

            px_[13 * n + 6] = x_[n] - l_;
            py_[13 * n + 6] = y_[n] + d;

            px_[13 * n + 7] = x_[n] - l_;
            py_[13 * n + 7] = y_[n] - d;

            px_[13 * n + 8] = x_[n] - d;
            py_[13 * n + 8] = y_[n] - d;

            px_[13 * n + 9] = x_[n] - d;
            py_[13 * n + 9] = y_[n] - l_;

            px_[13 * n + 10] = x_[n] + d;
            py_[13 * n + 10] = y_[n] - l_;

            px_[13 * n + 11] = x_[n] + d;
            py_[13 * n + 11] = y_[n] - d;

            px_[13 * n + 12] = px_[13 * n];
            py_[13 * n + 12] = py_[13 * n];

            ve_x_[12 * n] = 0.0;
            ve_y_[12 * n] = 1.0;

            ve_x_[12 * n + 1] = -1.0;
            ve_y_[12 * n + 1] =  0.0;

            ve_x_[12 * n + 2] = 0.0;
            ve_y_[12 * n + 2] = 1.0;

            ve_x_[12 * n + 3] = -1.0;
            ve_y_[12 * n + 3] =  0.0;

            ve_x_[12 * n + 4] =  0.0;
            ve_y_[12 * n + 4] = -1.0;

            ve_x_[12 * n + 5] = -1.0;
            ve_y_[12 * n + 5] =  0.0;

            ve_x_[12 * n + 6] =  0.0;
            ve_y_[12 * n + 6] = -1.0;

            ve_x_[12 * n + 7] = 1.0;
            ve_y_[12 * n + 7] = 0.0;

            ve_x_[12 * n + 8] = 0.0;
            ve_y_[12 * n + 8] = -1.0;

            ve_x_[12 * n + 9] = 1.0;
            ve_y_[12 * n + 9] = 0.0;

            ve_x_[12 * n + 10] = 0.0;
            ve_y_[12 * n + 10] = 1.0;

            ve_x_[12 * n + 11] = 1.0;
            ve_y_[12 * n + 11] = 0.0;
        }
        else
        {
            px_[13 * n] = l_;
            py_[13 * n] = -d;

            px_[13 * n + 1] = l_;
            py_[13 * n + 1] = d;

            px_[13 * n + 2] = d;
            py_[13 * n + 2] = d;

            px_[13 * n + 3] = d;
            py_[13 * n + 3] = l_;

            px_[13 * n + 4] = -d;
            py_[13 * n + 4] = l_;

            px_[13 * n + 5] = -d;
            py_[13 * n + 5] = d;

            px_[13 * n + 6] = -l_;
            py_[13 * n + 6] =  d;

            px_[13 * n + 7] = -l_;
            py_[13 * n + 7] = -d;

            px_[13 * n + 8] = -d;
            py_[13 * n + 8] = -d;

            px_[13 * n + 9] = -d;
            py_[13 * n + 9] = -l_;

            px_[13 * n + 10] = d;
            py_[13 * n + 10] = -l_;

            px_[13 * n + 11] = d;
            py_[13 * n + 11] = -d;

            ve_x_[12 * n] = 0.0;
            ve_y_[12 * n] = 1.0;

            ve_x_[12 * n + 1] = -1.0;
            ve_y_[12 * n + 1] =  0.0;

            ve_x_[12 * n + 2] = 0.0;
            ve_y_[12 * n + 2] = 1.0;

            ve_x_[12 * n + 3] = -1.0;
            ve_y_[12 * n + 3] =  0.0;

            ve_x_[12 * n + 4] =  0.0;
            ve_y_[12 * n + 4] = -1.0;

            ve_x_[12 * n + 5] = -1.0;
            ve_y_[12 * n + 5] =  0.0;

            ve_x_[12 * n + 6] =  0.0;
            ve_y_[12 * n + 6] = -1.0;

            ve_x_[12 * n + 7] = 1.0;
            ve_y_[12 * n + 7] = 0.0;

            ve_x_[12 * n + 8] =  0.0;
            ve_y_[12 * n + 8] = -1.0;

            ve_x_[12 * n + 9] = 1.0;
            ve_y_[12 * n + 9] = 0.0;

            ve_x_[12 * n + 10] = 0.0;
            ve_y_[12 * n + 10] = 1.0;

            ve_x_[12 * n + 11] = 1.0;
            ve_y_[12 * n + 11] = 0.0;


            double c = std::cos ( alpha );
            double s = std::sin ( alpha );

            for ( auto i = 0; i < 12; i++ )
            {
                double px_new = px_[13 * n + i] * c - py_[13 * n + i] * s;
                double py_new = px_[13 * n + i] * s + py_[13 * n + i] * c;

                px_[13 * n + i] = px_new + x_[n];
                py_[13 * n + i] = py_new + y_[n];

                px_new = ve_x_[12 * n + i] * c - ve_y_[12 * n + i] * s;
                py_new = ve_x_[12 * n + i] * s + ve_y_[12 * n + i] * c;

                ve_x_[12 * n + i] = px_new;
                ve_y_[12 * n + i] = py_new;

            }

            px_[13 * n + 12] = px_[13 * n];
            py_[13 * n + 12] = py_[13 * n];
        }
    }

    // between nr_obs1 and all < nr_obs1
    inline int overlapp ( int nr_obs1, double shift_x, double shift_y )
    {
        int index1 =  13 * nr_obs1;

        for ( int i = 0 ; i < 12; i++ )
        {
            int side1 = index1 + i;

            double P1x = px_[side1] + shift_x;
            double P1y = py_[side1] + shift_y;

            double P2x = px_[side1 + 1] + shift_x;
            double P2y = py_[side1 + 1] + shift_y;

            for ( int j = 0 ; j < nr_obs1 ; j++ )
            {
                int index2 = j * 13;

                //std::cout << std::sqrt(dx * dx + dy *dy) << " d" << std::endl;

                for ( int k = 0; k < 12 ; k++ )
                {
                    int side2 = index2 + k;

                    double P3x = px_[side2];
                    double P3y = py_[side2];

                    double P4x = px_[side2 + 1];
                    double P4y = py_[side2 + 1];

                    double Ax = P2x - P1x;
                    double Ay = P2y - P1y;

                    double Bx = P3x - P4x;
                    double By = P3y - P4y;

                    double Cx = P1x - P3x;
                    double Cy = P1y - P3y;

                    double den = Ay * Bx - Ax * By;
                    double num1 = By * Cx - Bx * Cy;
                    double num2 = Ax * Cy - Ay * Cx;

                    if ( den == 0.0 )
                    {
                        continue;
                    }

                    double len1 = num1 / den;
                    double len2 = num2 / den;

                    if ( len1 > 1.0 or len1 < 0.0 or len2 > 1.0 or len2 < 0.0 )
                    {
                        continue;
                    }
                    else
                    {
                        return -1;
                    }
                }
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
        obstaclesFile << 12 << std::endl;
        obstaclesFile.precision ( 8 );
        
        for ( int i = 0; i < n_tot_; ++i )
        {
            for ( int j = 0; j < 12; j++ )
            {
                obstaclesFile << i << ' ' << px_[13*i+j] << ' ' << py_[13*i+j] << std::endl;
            }
        }
        obstaclesFile.close();
    }

};//Class
};//NameSpace
