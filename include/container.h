#pragma once

#include <fstream>
#include <iostream>
#include <random>

#include "obstacles.h"
#include "globals.h"
#include "shapes_segments.h"

namespace Neighbourlist
{
class Linear;
}

class Container
{
public:

    double lx_;
    double ly_;

    //double surface_;

    Shapes::Segment *segments;

    //contructors and destructors
    Container ( double lx, double ly )
    {
        lx_ = lx;
        ly_ = ly;

        segments = new Shapes::Segment ( 4 );
        segments->setVertex ( 0, 0, 0 );
        segments->setVertex ( lx, 0, 1 );
        segments->setVertex ( lx, ly, 2 );
        segments->setVertex ( 0, ly, 3 );
        segments->setVertex ( 0, 0, 4 );
        segments->initPath();

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Container---" << std::endl;
            std::cout << "lx " << lx_ << std::endl;
            std::cout << "ly " << ly_ << std::endl;
        }
    }

    ~Container()
    {
        delete segments;
    }

    //functions to inline

    inline bool isOutside ( double x, double y )
    {
        if ( x <  0.0 )
        {
            return true;
        }
        if ( x >= lx_ )
        {
            return true;
        }
        if ( y <  0.0 )
        {
            return true;
        }
        if ( y >= ly_ )
        {
            return true;
        }
        return false;
    }

    inline bool isInside ( double x, double y )
    {
        if ( x <  0.0 )
        {
            return false;
        }
        if ( x >= lx_ )
        {
            return false;
        }
        if ( y <  0.0 )
        {
            return false;
        }
        if ( y >= ly_ )
        {
            return false;
        }
        return true;
    }

    inline int checkBoundary ( double &x, double &y )
    {
        int flag = 0;

        if ( x <  0.0 )
        {
            flag -= 1;
        }

        if ( x >= lx_ )
        {
            flag += 1;
        }

        if ( y <  0.0 )
        {
            flag -= 3;
        }

        if ( y >= ly_ )
        {
            flag += 3;
        }

        return flag;
    }

    inline void applyBoundary ( double &x, double &y )
    {
        if ( x < 0 )
        {
            x += lx_;
        }
        if ( x >= lx_ )
        {
            x -= lx_;
        }
        if ( y < 0 )
        {
            y += ly_;
        }
        if ( y >= ly_ )
        {
            y -= ly_;
        }
    }

    inline void applyBoundary ( double &x, double &y, int flag )
    {
        if ( flag == 0 )
        {
            return;
        }

        if ( flag == -1 )
        {
            x += lx_;
        }
        if ( flag == 1 )
        {
            x -= lx_;
        }
        if ( flag == -3 )
        {
            y += ly_;
        }
        if ( flag == 3 )
        {
            y -= ly_;
        }
        if ( flag == 2 )
        {
            x += lx_;
            y -= ly_;
        }
        if ( flag == 4 )
        {
            x -= lx_;
            y -= ly_;
        }
        if ( flag == -2 )
        {
            x -= lx_;
            y += ly_;
        }
        if ( flag == -4 )
        {
            x += lx_;
            y += ly_;
        }
    }

    inline void remapInside ( double &x, double &y )
    {
        x = std::fmod ( x, lx_ );
        y = std::fmod ( y, ly_ );
    }

    void minimalImage ( double &dx, double &dy )
    {
        if ( dx <= -0.5 * lx_ )
        {
            dx += lx_;
        }
        if ( dx >   0.5 * lx_ )
        {
            dx -= lx_;
        }

        if ( dy <= -0.5 * ly_ )
        {
            dy += ly_;
        }
        if ( dy >   0.5 * ly_ )
        {
            dy -= ly_;
        }
    }

    void setRandomPositions ( Obstacles::Circles &obstacles )
    {
        //initializing random number generator
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        //initializing random position in the box and setting the radius
        for ( int i = 0; i < obstacles.n_obstacles_; ++i )
        {
            obstacles.x_[i] = uni_dist ( Globalrand::e1 ) * lx_;
            obstacles.y_[i] = uni_dist ( Globalrand::e1 ) * ly_;
        }
    }

    void setRandomPositions ( Obstacles::Squares &obstacles )
    {
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        //initializing random position in the box and setting the radius
        if ( obstacles.angle_ < 0 ) // Random orientation
        {
            for ( int i = 0; i < obstacles.n_obstacles_; ++i )
            {
                obstacles.set ( i, uni_dist ( Globalrand::e1 ) * M_PI * 0.5, uni_dist ( Globalrand::e1 ) * lx_, uni_dist ( Globalrand::e1 ) * ly_ );
            }
        }
        else
        {
            for ( int i = 0; i < obstacles.n_obstacles_; ++i )
            {
                obstacles.set ( i, obstacles.angle_, uni_dist ( Globalrand::e1 ) * lx_, uni_dist ( Globalrand::e1 ) * ly_ );
            }
        }
    }

    void setRandomPositions ( Obstacles::Crosses &obstacles )
    {
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        //initializing random position in the box and setting the radius
        if ( obstacles.angle_ < 0 ) // Random orientation
        {
            for ( int i = 0; i < obstacles.n_obstacles_; ++i )
            {
                obstacles.set ( i, uni_dist ( Globalrand::e1 ) * M_PI * 0.5, uni_dist ( Globalrand::e1) * lx_, uni_dist ( Globalrand::e1 ) * ly_ );
            }
        }
        else
        {
            for ( int i = 0; i < obstacles.n_obstacles_; ++i )
            {
                obstacles.set ( i, obstacles.angle_, uni_dist ( Globalrand::e1 ) * lx_, uni_dist ( Globalrand::e1 ) * ly_ );
            }
        }
    }

    void setRandomPositionsNoOverlapp ( Obstacles::Circles &obstacles )
    {
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        int n = 0;

        while ( n < obstacles.n_obstacles_ )
        {
            //adding
            obstacles.set ( n,  uni_dist ( Globalrand::e1 ) * lx_, uni_dist ( Globalrand::e1 ) * ly_ );
            //check overlapp
            if ( obstacles.overlapp ( n, 0.0, 0.0 ) < 0 ) //reject
            {
                continue;
            }
            else // we have to check the ghost image:
            {
                int flag = 0;

                if ( obstacles.x_[n] > lx_ - 2.0 * obstacles.max_radius_ ) // right border
                {
                    flag += 1;
                }

                if ( obstacles.x_[n] < 2.0 * obstacles.max_radius_ )     // left border
                {
                    flag -= 1;
                }

                if ( obstacles.y_[n] > ly_ - 2.0 * obstacles.max_radius_ ) // upper border
                {
                    flag += 3;
                }
                if ( obstacles.y_[n] < 2.0 * obstacles.max_radius_ )     // lowerborder
                {
                    flag -= 3;
                }

                switch ( flag )
                {
                case 1:
                {
                    if ( obstacles.overlapp ( n, -lx_, 0.0 ) < 0 ) //reject
                    {
                        continue;
                    }
                    break;
                }
                case -1:
                {
                    if ( obstacles.overlapp ( n, lx_, 0.0 ) < 0 ) //reject
                    {
                        continue;
                    }

                    break;
                }

                case 3:
                {
                    if ( obstacles.overlapp ( n, 0.0, -ly_ ) < 0 ) //reject
                    {
                        continue;
                    }

                    break;
                }
                case -3:
                {
                    if ( obstacles.overlapp ( n, 0.0, +ly_ ) < 0 ) //reject
                    {
                        continue;
                    }

                    break;
                }

                case 2:
                {
                    if ( obstacles.overlapp ( n, +lx_, -ly_ ) < 0 or obstacles.overlapp ( n, 0.0, -ly_ ) < 0 or obstacles.overlapp ( n, lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;
                }

                case 4:
                {
                    if ( obstacles.overlapp ( n, -lx_, -ly_ ) < 0 or obstacles.overlapp ( n, 0.0, -ly_ ) < 0 or obstacles.overlapp ( n, -lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;
                }

                case -2:
                {
                    if ( obstacles.overlapp ( n, -lx_, +ly_ ) < 0 or obstacles.overlapp ( n, 0.0, +ly_ ) < 0 or obstacles.overlapp ( n, -lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;

                }

                case -4:
                {
                    if ( obstacles.overlapp ( n, +lx_, +ly_ ) < 0 or obstacles.overlapp ( n, 0.0, +ly_ ) < 0 or obstacles.overlapp ( n, lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;
                }
                }
                n++;
            }
        }
        obstacles.n_obstacles_ = n;
    }

    void setRandomPositionsNoOverlappMontecarlo ( Obstacles::Squares &obstacles, int sqrt_n )
    {
        double d_x = lx_ / sqrt_n;
        
        // we should have resized the container in such a way that nx*nx is exacliy a Squares

        if ( obstacles.angle_ < 0 ) // Random orientation
        {
            if ( d_x < obstacles.max_radius_ )
            {
                std::cout << "DENSITY TO HIGH FOR NO OVERLAPP" << std::endl;
                return;
            }
        }
        else
        {
            if ( d_x < obstacles.l_ )
            {
                std::cout << "DENSITY TO HIGH FOR NO OVERLAPP" << std::endl;
                return;
            }
        }

        double d_offset = 0.5 * d_x;

        for ( int i = 0; i < sqrt_n; i++ )
        {
            for ( int j = 0; j < sqrt_n; j++ )
            {
                if ( i + sqrt_n * j < obstacles.n_obstacles_ )
                {
                    obstacles.set ( i + sqrt_n * j, 0, d_offset + i * d_x , d_offset + j * d_x);
                }
            }
        }
    }

    void setRandomPositionsNoOverlapp ( Obstacles::Squares &obstacles )
    {
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        int n = 0;

        while ( n < obstacles.n_obstacles_ )
        {
            //adding
            if ( obstacles.angle_ < 0 ) // Random orientation
            {
                obstacles.set ( n, uni_dist ( Globalrand::e1 ) * M_PI * 0.25, uni_dist ( Globalrand::e1 ) * lx_, uni_dist ( Globalrand::e1 ) * ly_ );
            }
            else
            {
                obstacles.set ( n, obstacles.angle_, uni_dist ( Globalrand::e1 ) * lx_, uni_dist ( Globalrand::e1 ) * ly_ );
            }
            //check overlapp
            if ( obstacles.overlapp ( n, 0.0, 0.0 ) < 0 ) //reject
            {
                continue;
            }
            else // we have to check the ghost image:
            {
                int flag = 0;

                if ( obstacles.x_[n] > lx_ - 2.0 * obstacles.max_radius_ ) // right border
                {
                    flag += 1;
                }

                if ( obstacles.x_[n] < 2.0 * obstacles.max_radius_ )     // left border
                {
                    flag -= 1;
                }

                if ( obstacles.y_[n] > ly_ - 2.0 * obstacles.max_radius_ ) // upper border
                {
                    flag += 3;
                }
                if ( obstacles.y_[n] < 2.0 * obstacles.max_radius_ )     // lowerborder
                {
                    flag -= 3;
                }

                switch ( flag )
                {
                case 1:
                {
                    if ( obstacles.overlapp ( n, -lx_, 0.0 ) < 0 ) //reject
                    {
                        continue;
                    }
                    break;
                }
                case -1:
                {
                    if ( obstacles.overlapp ( n, lx_, 0.0 ) < 0 ) //reject
                    {
                        continue;
                    }

                    break;
                }

                case 3:
                {
                    if ( obstacles.overlapp ( n, 0.0, -ly_ ) < 0 ) //reject
                    {
                        continue;
                    }

                    break;
                }
                case -3:
                {
                    if ( obstacles.overlapp ( n, 0.0, +ly_ ) < 0 ) //reject
                    {
                        continue;
                    }

                    break;
                }

                case 2:
                {
                    if ( obstacles.overlapp ( n, +lx_, -ly_ ) < 0 or obstacles.overlapp ( n, 0.0, -ly_ ) < 0 or obstacles.overlapp ( n, lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;
                }

                case 4:
                {
                    if ( obstacles.overlapp ( n, -lx_, -ly_ ) < 0 or obstacles.overlapp ( n, 0.0, -ly_ ) < 0 or obstacles.overlapp ( n, -lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;
                }

                case -2:
                {
                    if ( obstacles.overlapp ( n, -lx_, +ly_ ) < 0 or obstacles.overlapp ( n, 0.0, +ly_ ) < 0 or obstacles.overlapp ( n, -lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;

                }

                case -4:
                {
                    if ( obstacles.overlapp ( n, +lx_, +ly_ ) < 0 or obstacles.overlapp ( n, 0.0, +ly_ ) < 0 or obstacles.overlapp ( n, lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;
                }
                }
                n++;

            }
        }
        obstacles.n_obstacles_ = n;
    }

    void setRandomPositionsNoOverlapp ( Obstacles::Crosses &obstacles )
    {
        //initializing random number generator
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        int n = 0;

        while ( n < obstacles.n_obstacles_ )
        {
            //adding
            if ( obstacles.angle_ < 0 ) // Random orientation
            {
                obstacles.set ( n, uni_dist ( Globalrand::e1 ) * M_PI * 0.25, uni_dist ( Globalrand::e1 ) * lx_, uni_dist ( Globalrand::e1 ) * ly_ );
            }
            else
            {
                obstacles.set ( n, obstacles.angle_, uni_dist ( Globalrand::e1 ) * lx_, uni_dist ( Globalrand::e1 ) * ly_ );
            }
            //check overlapp
            if ( obstacles.overlapp ( n, 0.0, 0.0 ) < 0 ) //reject
            {
                continue;
            }
            else // we have to check the ghost image:
            {
                int flag = 0;

                if ( obstacles.x_[n] > lx_ - 2.0 * obstacles.max_radius_ ) // right border
                {
                    flag += 1;
                }

                if ( obstacles.x_[n] < 2.0 * obstacles.max_radius_ )     // left border
                {
                    flag -= 1;
                }

                if ( obstacles.y_[n] > ly_ - 2.0 * obstacles.max_radius_ ) // upper border
                {
                    flag += 3;
                }
                if ( obstacles.y_[n] < 2.0 * obstacles.max_radius_ )     // lowerborder
                {
                    flag -= 3;
                }

                switch ( flag )
                {
                case 1:
                {
                    if ( obstacles.overlapp ( n, -lx_, 0.0 ) < 0 ) //reject
                    {
                        continue;
                    }
                    break;
                }
                case -1:
                {
                    if ( obstacles.overlapp ( n, lx_, 0.0 ) < 0 ) //reject
                    {
                        continue;
                    }

                    break;
                }

                case 3:
                {
                    if ( obstacles.overlapp ( n, 0.0, -ly_ ) < 0 ) //reject
                    {
                        continue;
                    }

                    break;
                }
                case -3:
                {
                    if ( obstacles.overlapp ( n, 0.0, +ly_ ) < 0 ) //reject
                    {
                        continue;
                    }

                    break;
                }

                case 2:
                {
                    if ( obstacles.overlapp ( n, +lx_, -ly_ ) < 0 or obstacles.overlapp ( n, 0.0, -ly_ ) < 0 or obstacles.overlapp ( n, lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;
                }

                case 4:
                {
                    if ( obstacles.overlapp ( n, -lx_, -ly_ ) < 0 or obstacles.overlapp ( n, 0.0, -ly_ ) < 0 or obstacles.overlapp ( n, -lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;
                }

                case -2:
                {
                    if ( obstacles.overlapp ( n, -lx_, +ly_ ) < 0 or obstacles.overlapp ( n, 0.0, +ly_ ) < 0 or obstacles.overlapp ( n, -lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;

                }

                case -4:
                {
                    if ( obstacles.overlapp ( n, +lx_, +ly_ ) < 0 or obstacles.overlapp ( n, 0.0, +ly_ ) < 0 or obstacles.overlapp ( n, lx_, 0.0 ) < 0 )
                    {
                        continue;
                    }
                    break;
                }
                }
                n++;
            }
        }
        obstacles.n_obstacles_ = n;
    }
};

