#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "cell_list.h"
#include "trajectory.h"
#include <vector>
#include "globals.h"

/****************************************************************************
 * int64_t run from −9.223.372.036.854.775.808 to 9.223.372.036.854.775.808
 * this means as we need the squared radius ie. r_sq << int64_max
 * -> r << 2^31.5 -> r << 30370005000 (3.0e9)
 ****************************************************************************/

namespace Stencil
{

class Circular
{
public:

    //std::ofstream debugfile;
    
    const int64_t nx_;
    const int64_t ny_;

    const double cell_l_;

    int n_big_;
    int max_lenght_;
    double r_red_sq;

    int64_t *cells_y_;
    int64_t *cells_x_;
    
public:

    Circular ( CellList &celllist, Container &container, Trajectory::Circular &trajectory ) :
        nx_ ( celllist.nx_ ),
        ny_ ( celllist.ny_ ),
        cell_l_ ( celllist.l_ )
    {
        
        //----
        //debugfile.open("debug.txt");
        //----
        
        r_red_sq = ( trajectory.radius_ / cell_l_ );
        r_red_sq *= r_red_sq;
        n_big_ = nx_ > ny_ ? nx_ : ny_;
        double l_big_ = container.lx_ > container.ly_ ? container.lx_ : container.ly_;
        int64_t N_r = ( trajectory.radius_ ) / cell_l_ + 1;

        if ( trajectory.radius_ < l_big_ )
        {
            max_lenght_ = 8 * N_r * 7 + 14 * 10 + 8 * 4; // ... + > ( 8 * 4 corners + 12 * 7 begin and end)
        }
        else if ( trajectory.radius_ < 2 * l_big_ )
        {
            max_lenght_ = 4 * 7 * ( n_big_ + 12 );
        }
        else
        {
            max_lenght_ = 2 * 7 * ( n_big_ + 12 );
        }

        max_lenght_ *= 1.10; // adding some space ... as the approximation is a little thight

        cells_x_ = new int64_t[max_lenght_];
        cells_y_ = new int64_t[max_lenght_];

        if ( GlobalFlags::verbose )
        {
            double mem_size = sizeof ( int64_t ) * max_lenght_ * 2;
            std::cout << std::endl << "---Stencil(Circular)---" << std::endl;
            printf ( "the stencil hast a maximum of %d cells\n", max_lenght_ );
            printf ( "total heapsize of the celllist is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    ~Circular()
    {
        delete[]cells_x_;
        delete[]cells_y_;
    }

///Bresenham's circle drawing, for reference see https://doi.org/10.1145/359423.359432

    inline void bresenham_step ( int64_t &Xtemp, int64_t &Ytemp, int64_t &Xs, int64_t &Ys, double &delta, int &Q, int D, int &m1x, int &m1y, int &m2x, int &m2y, int &m3x, int &m3y, int &flag )
    {
        if ( Ys < 1 ) //quadrant crossing
        {
            Q--;
            delta = delta - 4 * Xs;

            std::swap ( Xs, Ys );
            Xs *= -1;

            m1x = m3x;
            m1y = m3y;

            std::swap ( m2x, m2y );

            m2x *= D;
            m2y *= -D;

            std::swap ( m3x, m3y );

            m3x *= D;
            m3y *= -D;
        }

        if ( delta <= 0 )
        {
            double del = 2 * delta + 2 * Ys - 1;

            if ( del <= 0 )
            {
                flag = 1;
                Xtemp += m1x;
                Ytemp += m1y;
            }
            else
            {
                flag = 2;
                Xtemp += m2x;
                Ytemp += m2y;
            }
        }
        else
        {
            double del = 2 * delta - 2 * Xs - 1;

            if ( del <= 0 )
            {
                flag = 2;
                Xtemp += m2x;
                Ytemp += m2y;
            }
            else
            {
                flag = 3;
                Xtemp += m3x;
                Ytemp += m3y;
            }
        }

        switch ( flag )
        {
        case 1:
        {
            ++Xs;
            delta += 2 * Xs + 1;
            break;
        }

        case 2:
        {
            ++Xs;
            --Ys;
            delta = delta + 2 * Xs - 2 * Ys + 2;
            break;
        }

        case 3:
        {
            --Ys;
            delta = delta - 2 * Ys + 1;
            break;
        }
        }

    }

    inline int getCells ( double xs, double ys, double rx, double ry ) //full circle
    {
        //std::cout << "full" << std::endl;
        int m1x, m1y, m2x, m2y, m3x, m3y;
        int64_t Xs, Ys, Xt, Yt;
        double delta;

        //cell offset:
        int64_t offsetX = std::round ( rx/cell_l_ );
        int64_t offsetY = std::round ( ry/cell_l_ );

        setPoint ( xs, ys, rx, ry, Xs, Ys );

        int64_t Xtemp = Xs;
        int64_t Ytemp = Ys;

        int Q;
        int flag = 0;
        // now we shift the beginning point cw wise:
        normalize_cw ( Xs, Ys, Q, m1x, m1y, m2x, m2y, m3x, m3y );

        //initialize the error
        delta = ( Xs+1 ) * ( Xs+1 ) + ( Ys-1 ) * ( Ys-1 ) - r_red_sq;

        //Shift the beginning by 6 in the opposite direction to get the beginning cell
        bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q, 1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
        bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q, 1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
        bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q, 1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
        bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q, 1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
        bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q, 1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
        bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q, 1, m1x, m1y, m2x, m2y, m3x, m3y, flag );

        //set the new beggining point
        Xs = Xtemp;
        Ys = Ytemp;

        //normalise to ccw:
        normalize_ccw ( Xs, Ys, Q, m1x, m1y, m2x, m2y, m3x, m3y );

        Xt = Xs;
        Yt = Ys;

        //write cells
        Q = 3;
        delta = 2 * ( Xs - Ys + 1 ); //reset radius to exact value at the point to avoid init spike

        int count = 0;
        flag = 0;

        while ( Q >= 0 or Xt > Xs  or  Yt < Ys )
        {
            writeCells ( Xtemp, Ytemp, count, m3x, m3y, flag );
            bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta,Q, -1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
        }

        //add offset:
        for ( int i = 0; i < count; ++i )
        {
            cells_x_[i] += offsetX;
            cells_y_[i] += offsetY;
        }
        return count;
    }

    inline int getCells ( double xs, double ys, double xt, double yt, double rx, double ry )
    {
        //Compute normalized
        int qs, qt;
        int m1x, m1y, m2x, m2y, m3x, m3y;

        bool complete = false;

        int64_t Xs, Ys, Xt, Yt;
        double delta;

        //cell offset:
        int64_t offsetX = std::round ( rx/cell_l_ );
        int64_t offsetY = std::round ( ry/cell_l_ );

        setPoint ( xs, ys, rx, ry, Xs, Ys );
        setPoint ( xt, yt, rx, ry, Xt, Yt );

        //not normalized starting and ending points

        int64_t Xtemp = Xs;
        int64_t Ytemp = Ys;

        int64_t Xend = Xt;
        int64_t Yend = Yt;

        //calculate original number of quadrant crossings
        normalize_ccw ( Xs, Ys, qs, m1x, m1y, m2x, m2y, m3x, m3y );
        normalize_ccw ( Xt, Yt, qt );

        int Q_reverse = 0;
        int Q = ( qt - qs ) % 4;
        double vp = ( xs-rx ) * ( yt-ry ) - ( ys-ry ) * ( xt-rx );
        
        if ( Q < 0 )
        {
            Q += 4;
        }

        if ( Q == 0 and vp < 0.0 ) // vector product
        {
            Q = 3;
        }
        else
        {
            --Q;
        }
        // now we shift the beginning point cw wise:
        Xs = Xtemp;
        Ys = Ytemp;

        normalize_cw ( Xs, Ys, qs, m1x, m1y, m2x, m2y, m3x, m3y );
        //initialize the error
        delta = ( Xs+1 ) * ( Xs+1 ) + ( Ys-1 ) * ( Ys-1 ) - r_red_sq;
        
        int flag ( 0 );
        for ( int i = 0; i < 6 ; ++i )
        {
            bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q_reverse, 1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
        }
        
        //now check if the end is too close
        if ( Q > 1 )
        {
            if ( std::abs ( Xtemp-Xend ) < 7 and std::abs ( Ytemp-Yend ) < 7 )
            {
                
                //std::cout << Q << "  " << Xtemp << ' ' << Ytemp << "  " << Xend << ' ' << Yend << std::endl;
                complete = true;
            }
        }

        int count = 0;
        if ( complete == false ) // draw
        {
            //std::cout << "not complete" << std::endl;
            Xs = Xtemp;
            Ys = Ytemp;

            //normalise to ccw:
            normalize_ccw ( Xs, Ys, qs, m1x, m1y, m2x, m2y, m3x, m3y );

            //recalculate Q! 
            Q = ( qt - qs ) % 4;
            if ( Q < 0 )
            {
                Q += 4;
            }
            
            if ( Q == 0 and vp < 0.0 ) // vector product
            {
                Q = 3;
            }
            else
            {
                --Q;
            }
            //write cells
            int flag ( 0 );
            delta = ( Xs+1 ) * ( Xs+1 ) + ( Ys-1 ) * ( Ys-1 ) - r_red_sq;
            while ( Q >= 0 or Xt > Xs  or  Yt < Ys )
            {
                //debugfile << (Xtemp + offsetX) * cell_l_ << ' ' << (Ytemp + offsetY) * cell_l_ << std::endl;
                writeCells ( Xtemp, Ytemp, count, m3x, m3y, flag );
                bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q, -1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
            }

            for ( int i = 0; i< 6; ++i )
            {
                writeCells ( Xtemp, Ytemp, count, m3x, m3y, flag );
                bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q, -1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
            }
        }
        else
        {
            Q=3;
            Xs = Xtemp;
            Ys = Ytemp;

            //normalise to ccw:
            normalize_ccw ( Xs, Ys, qs, m1x, m1y, m2x, m2y, m3x, m3y );

            Xt = Xs;
            Yt = Ys;

            //write cells
            delta = 2 * ( Xs - Ys + 1 );
            int flag ( 0 );
            while ( Q >= 0 or Xt > Xs  or  Yt < Ys )
            {
                writeCells ( Xtemp, Ytemp, count, m3x, m3y, flag );
                bresenham_step ( Xtemp, Ytemp, Xs, Ys, delta, Q, -1, m1x, m1y, m2x, m2y, m3x, m3y, flag );
            }
        }

        //add offset:
        for ( int i = 0; i < count; ++i )
        {
            cells_x_[i] += offsetX;
            cells_y_[i] += offsetY;
        }
        return count;
    }

private:

//normalise to first quadrant

    void normalize_ccw ( int64_t &X, int64_t &Y, int &q )
    {
        if ( X < 0 )
        {
            if ( Y < 0 )
            {
                std::swap ( X, Y );
                X = -X;
                Y = -Y;
                q = 3;
            }
            else
            {
                X = -X;
                q = 2;
            }
        }
        else
        {
            if ( Y < 0 )
            {
                Y = -Y;
                q = 0;
            }
            else
            {
                std::swap ( X, Y );
                q = 1;
            }
        }
    }

    void normalize_ccw ( int64_t &X, int64_t &Y, int &q, int &m1x, int &m1y, int &m2x, int &m2y, int &m3x, int &m3y )
    {
        if ( X < 0 )
        {
            if ( Y < 0 )   // III
            {
                X = -X;
                Y = -Y;
                std::swap ( X, Y );

                q = 3;
                m1x = 0;
                m1y = -1;
                m2x = 1;
                m2y = -1;
                m3x = 1;
                m3y = 0;
            }
            else     // II
            {
                X = -X;
                q = 2;
                m1x = -1;
                m1y = 0;
                m2x = -1;
                m2y = -1;
                m3x = 0;
                m3y = -1;
            }
        }
        else
        {
            if ( Y < 0 )   // IV
            {
                Y = -Y;
                q = 0;
                m1x = 1;
                m1y = 0;
                m2x = 1;
                m2y = 1;
                m3x = 0;
                m3y = 1;
            }
            else     // I
            {
                std::swap ( X, Y );
                q = 1;
                m1x = 0;
                m1y = 1;
                m2x = -1;
                m2y = 1;
                m3x = -1;
                m3y = 0;
            }
        }
    }
    void normalize_ccw ( double &x, double &y )
    {
        if ( x < 0.0 )
        {
            if ( y < 0.0 )
            {
                std::swap ( x, y );
                x = -x;
                y = -y;

            }
            else
            {
                x = -x;
            }
        }
        else
        {
            if ( y < 0.0 )
            {
                y = -y;
            }
            else
            {
                std::swap ( x, y );
            }
        }
    }

    void normalize_cw ( int64_t &X, int64_t &Y, int &q )
    {
        if ( X < 0 )
        {
            if ( Y < 0 )
            {
                X = -X;
                Y = -Y;

                q = 2;
            }
            else
            {
                X = -X;
                std::swap ( X,Y );

                q = 3;
            }
        }
        else
        {
            if ( Y < 0 )
            {
                Y = -Y;
                std::swap ( X,Y );
                q = 1;
            }
            else
            {
                q = 0;
            }
        }
    }

    void normalize_cw ( int64_t &X, int64_t &Y, int &q, int &m1x, int &m1y, int &m2x, int &m2y, int &m3x, int &m3y )
    {
        if ( X < 0 )
        {
            if ( Y < 0 )   // III
            {
                X = -X;
                Y = -Y;

                q = 2;

                m1x = -1;
                m1y = 0;
                m2x = -1;
                m2y = 1;
                m3x = 0;
                m3y = 1;
            }
            else     // II
            {
                X = -X;
                std::swap ( X,Y );

                q = 3;

                m1x = 0;
                m1y = 1;
                m2x = 1;
                m2y = 1;
                m3x = 1;
                m3y = 0;
            }
        }
        else
        {
            if ( Y < 0 )   // IV
            {
                Y = -Y;
                std::swap ( X,Y );

                q = 1;

                m1x = 0;
                m1y = -1;
                m2x = -1;
                m2y = -1;
                m3x = -1;
                m3y = 0;
            }
            else     // I
            {
                q = 0;
                m1x = 1;
                m1y = 0;
                m2x = 1;
                m2y = -1;
                m3x = 0;
                m3y = -1;
            }
        }
    }

    void writeColumn ( int count, int64_t cell_x, int64_t cell_y )
    {
        cells_x_[count    ] = cell_x;
        cells_x_[count + 1] = cell_x;
        cells_x_[count + 2] = cell_x;
        cells_x_[count + 3] = cell_x;
        cells_x_[count + 4] = cell_x;
        cells_x_[count + 5] = cell_x;
        cells_x_[count + 6] = cell_x;

        cells_y_[count    ] = cell_y - 3;
        cells_y_[count + 1] = cell_y - 2;
        cells_y_[count + 2] = cell_y - 1;
        cells_y_[count + 3] = cell_y ;
        cells_y_[count + 4] = cell_y + 1;
        cells_y_[count + 5] = cell_y + 2;
        cells_y_[count + 6] = cell_y + 3;
    }

    inline void writeRow ( int count, int64_t cell_x, int64_t cell_y )
    {
        cells_y_[count    ] = cell_y;
        cells_y_[count + 1] = cell_y;
        cells_y_[count + 2] = cell_y;
        cells_y_[count + 3] = cell_y;
        cells_y_[count + 4] = cell_y;
        cells_y_[count + 5] = cell_y;
        cells_y_[count + 6] = cell_y;

        cells_x_[count    ] = cell_x - 3;
        cells_x_[count + 1] = cell_x - 2;
        cells_x_[count + 2] = cell_x - 1;
        cells_x_[count + 3] = cell_x ;
        cells_x_[count + 4] = cell_x + 1;
        cells_x_[count + 5] = cell_x + 2;
        cells_x_[count + 6] = cell_x + 3;
    }

    inline void writeRow_limit ( int &count, int64_t cell_x, int64_t cell_y )
    {
        if ( cell_x < 0 )
        {
            if ( cell_y < 0 )
            {
                cell_x -= 3;
                while ( cell_x <= cell_y )
                {
                    cells_x_[count] = cell_x;
                    cells_y_[count++] = cell_y;
                    ++cell_x;
                }
            }
            else
            {
                cell_x -= 3;
                while ( cell_x <= -cell_y )
                {
                    cells_x_[count] = cell_x;
                    cells_y_[count++] = cell_y;
                    ++cell_x;
                }
            }
        }
        else
        {
            if ( cell_y < 0 )
            {
                cell_x += 3;
                while ( cell_x >= -cell_y )
                {
                    cells_x_[count] = cell_x ;
                    cells_y_[count++] = cell_y;
                    --cell_x;
                }
            }
            else
            {
                cell_x += 3;
                while ( cell_x >= cell_y )
                {
                    cells_x_[count] = cell_x;
                    cells_y_[count++] = cell_y;
                    --cell_x;
                }
            }
        }
    }

    inline void writeColumn_limit ( int &count, int64_t cell_x, int64_t cell_y )
    {
        if ( cell_x < 0 )
        {
            if ( cell_y < 0 )
            {
                cell_y -= 3;
                while ( cell_y < cell_x )
                {
                    cells_x_[count] = cell_x;
                    cells_y_[count++] = cell_y;
                    ++cell_y;
                }
            }
            else
            {
                cell_y += 3;
                while ( cell_y > -cell_x )
                {
                    cells_x_[count] = cell_x;
                    cells_y_[count++] = cell_y;
                    --cell_y;
                }
            }
        }
        else
        {
            if ( cell_y < 0 )
            {
                cell_y -= 3;
                while ( cell_y < -cell_x )
                {
                    cells_x_[count] = cell_x ;
                    cells_y_[count++] = cell_y;
                    ++cell_y;
                }
            }
            else
            {
                cell_y += 3;
                while ( cell_y > cell_x )
                {
                    cells_x_[count] = cell_x;
                    cells_y_[count++] = cell_y;
                    --cell_y;
                }
            }
        }
    }

    inline void writeCorner1 ( int &count, int64_t cell_x, int64_t cell_y )
    {
        if ( cell_x != cell_y and - cell_x != cell_y )
        {
            return;
        }
        if ( cell_x < 0 )
        {
            if ( cell_y < 0 )
            {
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y;

                cells_x_[count] = cell_x - 1;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y - 1;

                cells_x_[count] = cell_x - 2;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y - 2;

                cells_x_[count] = cell_x - 1;
                cells_y_[count++] = cell_y - 1;
            }
            else
            {
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y;

                cells_x_[count] = cell_x - 1;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y + 1;

                cells_x_[count] = cell_x - 2;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y + 2;

                cells_x_[count] = cell_x - 1;
                cells_y_[count++] = cell_y + 1;
            }
        }
        else
        {
            if ( cell_y < 0 )
            {
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y;

                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y - 1;

                cells_x_[count] = cell_x + 2;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y - 2;

                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y - 1;
            }
            else
            {
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y;

                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y + 1;

                cells_x_[count] = cell_x + 2;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y + 2;

                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y + 1;
            }
        }
    }

    inline void writeCorner2 ( int &count, int64_t cell_x, int64_t cell_y )
    {
        if ( cell_x < 0 ) //count +10
        {
            if ( cell_y < 0 )
            {
                cells_x_[count] = cell_x;     // diag
                cells_y_[count++] = cell_y;

                cells_x_[count] = cell_x - 1; //1st line
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y - 1;

                cells_x_[count] = cell_x - 2; // 2nd line
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y - 2;
                cells_x_[count] = cell_x - 1;
                cells_y_[count++] = cell_y - 1;

                cells_x_[count] = cell_x - 3; // 3rd line
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x - 2;
                cells_y_[count++] = cell_y - 1;
                cells_x_[count] = cell_x - 1;
                cells_y_[count++] = cell_y - 2;
                cells_x_[count] = cell_x ;
                cells_y_[count++] = cell_y - 3 ;
            }
            else
            {
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y;

                cells_x_[count] = cell_x - 1;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y + 1;

                cells_x_[count] = cell_x - 2;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y + 2;
                cells_x_[count] = cell_x - 1;
                cells_y_[count++] = cell_y + 1;

                cells_x_[count] = cell_x - 3; // 3rd line
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x - 2;
                cells_y_[count++] = cell_y + 1;
                cells_x_[count] = cell_x - 1;
                cells_y_[count++] = cell_y + 2;
                cells_x_[count] = cell_x ;
                cells_y_[count++] = cell_y + 3;
            }
        }
        else
        {
            if ( cell_y < 0 )
            {
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y;

                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y - 1;

                cells_x_[count] = cell_x + 2;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y - 2;
                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y - 1;

                cells_x_[count] = cell_x + 3; // 3rd line
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x + 2;
                cells_y_[count++] = cell_y - 1;
                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y - 2;
                cells_x_[count] = cell_x ;
                cells_y_[count++] = cell_y - 3;
            }
            else
            {
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y;

                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y + 1;

                cells_x_[count] = cell_x + 2;
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x;
                cells_y_[count++] = cell_y + 2;
                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y + 1;

                cells_x_[count] = cell_x + 3; // 3rd line
                cells_y_[count++] = cell_y;
                cells_x_[count] = cell_x + 2;
                cells_y_[count++] = cell_y + 1;
                cells_x_[count] = cell_x + 1;
                cells_y_[count++] = cell_y + 2;
                cells_x_[count] = cell_x ;
                cells_y_[count++] = cell_y + 3;
            }
        }
    }

    void writeCells ( int64_t cellx, int64_t celly, int &count, int m3x, int m3y, int flag )
    {        
        int64_t x_comp = std::abs ( cellx );
        int64_t y_comp = std::abs ( celly );

        if ( x_comp > y_comp + 5 )
        {
            writeRow ( count, cellx, celly );
            count += 7;
            return;
        }

        if ( x_comp < y_comp - 5 )
        {
            writeColumn ( count, cellx, celly );
            count += 7;
            return;
        }

        //-----------------------------------------

        if ( x_comp > y_comp + 1 )
        {
            writeRow_limit ( count, cellx, celly );
            return;
        }
        if ( x_comp < y_comp - 1 )
        {
            writeColumn_limit ( count, cellx, celly );
            return;
        }
        //------------------------------------------

        if ( x_comp == y_comp + 1 )
        {
            writeRow_limit ( count, cellx, celly );
            if ( flag == 2 )
            {
                writeCorner1 ( count, cellx - m3x, celly - m3y );
            }
            return;
        }

        if ( x_comp == y_comp - 1 )
        {
            writeColumn_limit ( count, cellx, celly );
            if ( flag == 2 )
            {
                writeCorner1 ( count, cellx - m3x, celly - m3y );
            }
            return;
        }

        if ( x_comp == y_comp )
        {
            writeCorner2 ( count, cellx, celly );
            return;
        }
    }

    ///for reference see
    ///Fundamental Algorithms for Computer Graphics
    ///NATO Advanced Study Institute directed by J.E. Bresenham, R.A. Earnshaw, M.L.V.
    ///https://link.springer.com/book/10.1007%2F978-3-642-84574-1
    ///page 211

    void setPoint ( double cx, double cy, double rx, double ry, int64_t &X_tt, int64_t &Y_tt )
    {
        int64_t X = std::floor ( ( cx - rx ) / cell_l_ );
        int64_t Y = std::floor ( ( cy - ry ) / cell_l_ );

        X_tt = X;
        Y_tt = Y;

        double error_00 =  X * X + Y * Y - r_red_sq;
        //std::cout << "error00 " << error_00 << std::endl;

        double error_tt = error_00 + 2 * X + 1 ;
        //std::cout << "error01 " << error_tt << std::endl;

        if ( std::abs ( error_tt ) < std::abs ( error_00 ) )
        {
            error_00 = error_tt;
            X_tt = X + 1;
            Y_tt = Y;
        }
        error_tt = error_tt + 2 * Y + 1;
        //std::cout << "error11 " << error_tt << std::endl;
        if ( std::abs ( error_tt ) < std::abs ( error_00 ) )
        {
            error_00 = error_tt;
            X_tt = X + 1;
            Y_tt = Y + 1;
        }
        error_tt = error_tt - 2 * X - 1;
        //std::cout << "error01 " <<error_tt << std::endl;
        if ( std::abs ( error_tt ) < std::abs ( error_00 ) )
        {
            error_00 = error_tt;
            X_tt = X;
            Y_tt = Y + 1;
        }
    }
};
};
