#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

#include "globals.h"
#include <vector>
#include <sys/stat.h>

namespace Obstacles
{
class Squares
{
public:

    std::vector<double>x_;    // center position
    std::vector<double>y_;

    std::vector<double>ve_x_; // unit vector parrallel to edge
    std::vector<double>ve_y_;

    std::vector<double>px_;   // vertices of the squares
    std::vector<double>py_;
    
    std::vector<double> angles_;

    double l_;

    int n_tot_;
    int n_obstacles_;


    //-

    double l_sq_;
    double l_inv_;

    double max_radius_;
    int n_ghosts_;

    double angle_;
    double area_;
    static constexpr int sides_ = 4;
    static constexpr int type_ = 1;

    // for debuging
    //int * obs_ID;

    Squares ( double l, double angle, int n_obstacles, int n_ghosts )
    {
        l_ = l;
        area_ = l_*l_;
        l_sq_ = l * l;
        l_inv_ = 1.0 / l;

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

        px_.resize ( 5 * n_tot_ );
        py_.resize ( 5 * n_tot_ );

        ve_x_.resize ( 4 * n_tot_ );
        ve_y_.resize ( 4 * n_tot_ );

        // if unique angle is used
        angle_ = angle;

        if (angle < 0)
        {
            // if each obstacles has a random angle
            angles_.resize(n_tot_);
        }
        
        max_radius_ = l_ / std::sqrt ( 2.0 );

        double mem_size = sizeof ( double ) * 21 * n_tot_;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Obstacles (Squares)---" << std::endl;
            std::cout << "l_edge " << l_ << std::endl;
            std::cout << "angle " << angle << std::endl;
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

            double ap_x = x - px_[5 * n];
            double ap_y = y - py_[5 * n];

            double sp1 = ( ve_x_[4 * n] * ap_x + ve_y_[4 * n] * ap_y ) * l_;
            double sp2 = - ( ve_x_[4 * n + 3] * ap_x + ve_y_[4 * n + 3] * ap_y ) * l_;

            if ( sp1 > 0 and sp1 < l_ * l_ and sp2 > 0 and sp2 < l_sq_ )
            {
                return true;
            }
        }
        return false;
    }

    bool isInside ( int n, double x, double y )
    {
        double ap_x = x - px_[5 * n];
        double ap_y = y - py_[5 * n];

        double sp1 = ( ve_x_[4 * n] * ap_x + ve_y_[4 * n] * ap_y ) * l_;
        double sp2 = - ( ve_x_[4 * n + 3] * ap_x + ve_y_[4 * n + 3] * ap_y ) * l_;

        if ( sp1 > 0 and sp1 < l_ * l_ and sp2 > 0 and sp2 < l_sq_ )
        {
            return true;
        }

        return false;
    }



    void set ( int n, double alpha, double x, double y )
    {
        x_[n] = x;
        y_[n] = y;

        if(angle_ < 0)
        {
            angles_[n] = alpha;
        }
        
        if ( alpha == 0.0 )
        {
            double lh = 0.5 * l_;

            px_[5 * n] = x_[n] + lh;
            py_[5 * n] = y_[n] + lh;

            px_[5 * n + 1] = x_[n] - lh;
            py_[5 * n + 1] = y_[n] + lh;

            px_[5 * n + 2] = x_[n] - lh;
            py_[5 * n + 2] = y_[n] - lh;

            px_[5 * n + 3] = x_[n] + lh;
            py_[5 * n + 3] = y_[n] - lh;

            px_[5 * n + 4] = px_[5 * n];
            py_[5 * n + 4] = py_[5 * n];

            ve_x_[4 * n] = -1.0;
            ve_y_[4 * n] =  0.0;

            ve_x_[4 * n + 1] =  0.0;
            ve_y_[4 * n + 1] = -1.0;

            ve_x_[4 * n + 2] = 1.0;
            ve_y_[4 * n + 2] = 0.0;

            ve_x_[4 * n + 3] = 0.0;
            ve_y_[4 * n + 3] = 1.0;
        }
        else
        {
            double lh = 0.5 * l_;

            px_[5 * n] = lh;
            py_[5 * n] = lh;

            px_[5 * n + 1] = -lh;
            py_[5 * n + 1] = lh;

            px_[5 * n + 2] = -lh;
            py_[5 * n + 2] = -lh;

            px_[5 * n + 3] = lh;
            py_[5 * n + 3] = -lh;

            double c = std::cos ( alpha );
            double s = std::sin ( alpha );

            for ( auto i = 0; i < 4; i++ )
            {
                double px_new = px_[5 * n + i] * c - py_[5 * n + i] * s;
                double py_new = px_[5 * n + i] * s + py_[5 * n + i] * c;

                px_[5 * n + i] = px_new + x_[n];
                py_[5 * n + i] = py_new + y_[n];
            }

            px_[5 * n + 4] = px_[5 * n];
            py_[5 * n + 4] = py_[5 * n];

            ve_x_[4 * n] = ( px_[5 * n + 1] - px_[5 * n] ) * l_inv_;
            ve_y_[4 * n] = ( py_[5 * n + 1] - py_[5 * n] ) * l_inv_;

            ve_x_[4 * n + 1] = -ve_y_[4 * n];
            ve_y_[4 * n + 1] =  ve_x_[4 * n];

            ve_x_[4 * n + 2] = -ve_x_[4 * n];
            ve_y_[4 * n + 2] = -ve_y_[4 * n];

            ve_x_[4 * n + 3] =  ve_y_[4 * n];
            ve_y_[4 * n + 3] = -ve_x_[4 * n];

        }
    }
    void shift ( int n, double dx, double dy )
    {
        px_[5*n] += dx;
        px_[5*n+1] += dx;
        px_[5*n+2] += dx;
        px_[5*n+3] += dx;
        px_[5*n+4] += dx;

        py_[5*n] += dy;
        py_[5*n+1] += dy;
        py_[5*n+2] += dy;
        py_[5*n+3] += dy;
        py_[5*n+4] += dy;

        x_[n] += dx;
        y_[n] += dy;
    }


    // between nr_obs1 and all < nr_obs1
    inline int overlapp ( int nr_obs1, double shift_x, double shift_y )
    {
        int index1 =  5 * nr_obs1;

        for ( int i = 0 ; i < 4; i++ )
        {
            int side1 = index1 + i;

            double P1x = px_[side1] + shift_x;
            double P1y = py_[side1] + shift_y;

            double P2x = px_[side1 + 1] + shift_x;
            double P2y = py_[side1 + 1] + shift_y;

            for ( int j = 0 ; j < nr_obs1 ; j++ )
            {
                int index2 = j * 5;

                //std::cout << std::sqrt(dx * dx + dy *dy) << " d" << std::endl;

                for ( int k = 0; k < 4 ; k++ )
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

    inline int overlapp ( int nr_obs1, int nr_obs2, double shift_x, double shift_y )
    {
        int index1 =  5 * nr_obs1;

        for ( int i = 0 ; i < 4; i++ )
        {
            int side1 = index1 + i;

            double P1x = px_[side1] + shift_x;
            double P1y = py_[side1] + shift_y;

            double P2x = px_[side1 + 1] + shift_x;
            double P2y = py_[side1 + 1] + shift_y;


            int index2 = nr_obs2 * 5;


            for ( int k = 0; k < 4 ; k++ )
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
        return 1;
    }


    void writeObstacles ( std::string filename,  bool append )
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
        obstaclesFile << 4 << std::endl;
        obstaclesFile.precision ( 8 );

        for ( int i = 0; i < n_tot_; ++i )
        {
            for ( int j = 0; j < 4; j++ )
            {
                obstaclesFile << i << ' ' << px_[5 * i + j] << ' ' << py_[5 * i + j] << std::endl;
            }
        }
        obstaclesFile.close();
    }


    void writeObstaclepositions ( std::string filename )
    {
        std::ofstream obstaclesFile;
        obstaclesFile.open ( filename );

        obstaclesFile.precision ( 8 );

        for ( int i = 0; i < n_tot_; ++i )
        {
            obstaclesFile << x_[i] << ' ' << y_[i] << std::endl;
        }
        obstaclesFile.close();
    }


};
};

