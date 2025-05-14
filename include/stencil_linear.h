#pragma once

#include <iostream>
#include <iomanip>

#include "cell_list.h"
#include "trajectory.h"
#include "globals.h"

namespace Stencil
{

class Linear
{
public:
    CellList *celllist_;
    int *cells_y_;
    int *cells_x_;
    /*
          int *test_x_;
          int *test_y_;*/

    const int nx_, ny_;
    const double cell_l_;
    int n_big_;
    int max_lenght_;

public:

    Linear ( CellList &celllist ) :
        nx_ ( celllist.nx_ ),
        ny_ ( celllist.ny_ ),
        cell_l_ ( celllist.l_ )
    {
        celllist_ = &celllist;
        n_big_ = nx_ > ny_ ? nx_ : ny_;
        max_lenght_ = 5 * ( celllist_->nx_ > celllist_->ny_ ? celllist_->nx_ : celllist_->ny_ ) + 25;
        cells_x_ = new int[max_lenght_];
        cells_y_ = new int[max_lenght_];

//         test_x_ = new int[max_lenght_];
//         test_y_ = new int[max_lenght_];

        if ( GlobalFlags::verbose )
        {
            double mem_size = sizeof ( int ) * max_lenght_ * 10 + 4;
            std::cout << std::endl << "---Stencil(Linear)---" << std::endl;
            printf ( "the stencil hast a maximum of %d cells\n", max_lenght_ );
            printf ( "total heapsize of the celllist is %.2f kB \n", mem_size / 1000.0 );
        }

    }

    ~Linear()
    {
        delete[] cells_x_;
        delete[] cells_y_;
    }

    int getCells ( int cell_x1, int cell_y1, int cell_x2, int cell_y2, double vx, double vy ) //Bresenham's line drawing
    {
        if ( vy > 0 )
        {
            if ( vx > 0 )   // I
            {
                if ( vx > vy )
                {
                    //std::cout << "lineForwardUp" << std::endl;
                    return lineForwardUp ( cell_x1, cell_y1, cell_x2, cell_y2 );
                }
                else
                {
                    //std::cout << "lineForwardUp_T" << std::endl;
                    return lineForwardUp_T ( cell_x1, cell_y1, cell_x2, cell_y2 );
                }
            }
            else     // II
            {
                //cell_x1 -= 1;
                //cell_x2 -= 1;

                if ( -vx > vy )
                {
                    //std::cout << "lineBackwardUp" << std::endl;
                    return lineBackwardUp ( cell_x1, cell_y1, cell_x2, cell_y2 );
                }
                else
                {
                    //std::cout << "lineBackwardUp_T" << std::endl;
                    return lineBackwardUp_T ( cell_x1, cell_y1, cell_x2, cell_y2 );
                }
            }
        }
        else     // (vy < 0)
        {
            //cell_y1 -= 1;
            //cell_y2 -= 1;

            if ( vx > 0 )   // IV
            {
                if ( vx > -vy )
                {
                    //std::cout << "lineForwardDown" << std::endl;
                    return lineForwardDown ( cell_x1, cell_y1, cell_x2, cell_y2 );
                }
                else
                {
                    //std::cout << "lineForwardDown_T" << std::endl;
                    return lineForwardDown_T ( cell_x1, cell_y1, cell_x2, cell_y2 );
                }
            }
            else     // III
            {

                //cell_x1 -= 1;
                //cell_x2 -= 1;

                if ( -vx > -vy )
                {
                    //std::cout << "lineBackwardDown" << std::endl;
                    return lineBackwardDown ( cell_x1, cell_y1, cell_x2, cell_y2 );
                }
                else
                {
                    //std::cout << "lineBackwardDown_T" << std::endl;
                    return lineBackwardDown_T ( cell_x1, cell_y1, cell_x2, cell_y2 );
                }
            }
        }
        return 0;
    }

    int getCellsWindtree ( int cell_x1, int cell_y1, double vx, double vy )
    {
        //std::cout << std::setprecision(16) << vx << ' ' << vy << std::endl;
        if ( vx == 1.0 )
        {
            return writeNarrowYRowForward ( cell_x1, cell_y1 );
        }
        if ( vx == -1.0 )
        {
            return writeNarrowYRowBackward ( cell_x1, cell_y1 );
        }
        if ( vy == 1.0 )
        {
            return writeNarrowXRowUpward ( cell_x1, cell_y1 );
        }
        if ( vy == -1.0 )
        {
            return writeNarrowXRowDownward ( cell_x1, cell_y1 );
        }
        std::cout << "Direction Problem !! " << std::endl;
        return 0;
    }

private:

    //Bresenham Quadrants
    inline int lineForwardUp ( int x0, int y0, int x1, int y1 )
    {
        const int dx = x1 - x0;
        const int dy = y1 - y0;
        const int dx2 = 2 * dx;
        const int dy2 = 2 * dy;
        int err = dy2 - dx;

        const int y_lim = ny_ + 2;

        int count = 0;
        //write one row behind
        writeYRow ( count, x0 - 1, y0 );
        count += 5;

        for ( ;; )
        {
            writeYRow ( count, x0, y0 );

            x0++;
            count += 5;

            if ( x0 >= nx_ )
            {
                return count;
            }

            if ( err > 0 )
            {
                y0++;
                err -= dx2;
                if ( y0 > y_lim )
                {
                    return count;
                }
            }
            err += dy2;
        }
        return count;
    }

    inline int lineForwardDown ( int x0, int y0, int x1, int y1 )
    {
        const int dx = x1 - x0;
        const int dy = y0 - y1;
        const int dx2 = 2 * dx;
        const int dy2 = 2 * dy;

        int err = dy2 - dx;

        int count = 0;
        //write one row behind
        writeYRow ( count, x0 - 1, y0 );
        count += 5;

        for ( ;; )
        {
            writeYRow ( count, x0, y0 );

            x0++;
            count += 5;

            if ( x0 >= nx_ )
            {
                return count;
            }

            if ( err > 0 )
            {
                y0--;
                err -= dx2;
                if ( y0 < -2 )
                {
                    return count;
                }
            }
            err += dy2;
        }
        return count;
    }

    inline int lineForwardUp_T ( int x0, int y0, int x1, int y1 )
    {
        const int dx = x1 - x0;
        const int dy = y1 - y0;
        const int dx2 = 2 * dx;
        const int dy2 = 2 * dy;
        int err = dx2 - dy;
        const int x_lim = nx_ + 2;

        int count = 0;

        //write one row behind
        writeXRow ( count, x0, y0 - 1 );
        count += 5;

        for ( ;; )
        {
            writeXRow ( count, x0, y0 );

            y0++;
            count += 5;

            if ( y0 >= ny_ )
            {
                return count;
            }

            if ( err > 0 )
            {
                x0++;
                err -= dy2;
                if ( x0 > x_lim )
                {
                    return count;
                }
            }
            err += dx2;
        }
        return count;
    }

    inline int lineForwardDown_T ( int x0, int y0, int x1, int y1 )
    {
        const int dx = x1 - x0;
        const int dy = y0 - y1;
        const int dx2 = 2 * dx;
        const int dy2 = 2 * dy;
        int err = dx2 - dy;

        int count = 0;
        //write one row behind
        writeXRow ( count, x0, y0 + 1 );
        count += 5;

        for ( ;; )
        {
            writeXRow ( count, x0, y0 );

            y0--;
            count += 5;

            if ( y0 < 0 )
            {
                return count;
            }

            if ( err > 0 )
            {
                x0++;
                err -= dy2;
                if ( x0 >= nx_ + 2 )
                {
                    return count;
                }
            }
            err += dx2;
        }
        return count;
    }

    inline int lineBackwardUp ( int x0, int y0, int x1, int y1 )
    {
        const int dx = x0 - x1;
        const int dy = y1 - y0;
        const int dx2 = 2 * dx;
        const int dy2 = 2 * dy;
        int err = dy2 - dx;

        const int y_lim = ny_ + 2;
        int count = 0;

        //write one row behind
        writeYRow ( count, x0 + 1, y0 );
        count += 5;

        for ( ;; )
        {
            writeYRow ( count, x0, y0 );
            count += 5;

            x0--;
            if ( x0 < 0 )
            {
                return count;
            }

            if ( err > 0 )
            {
                y0++;
                err -= dx2;
                if ( y0 > y_lim )
                {
                    return count;
                }
            }
            err += dy2;
        }
        return count;
    }

    inline int lineBackwardUp_T ( int x0, int y0, int x1, int y1 )
    {
        const int dx = x0 - x1;
        const int dy = y1 - y0;
        const int dx2 = 2 * dx;
        const int dy2 = 2 * dy;
        int err = dx2 - dy;

        int count = 0;
        //write one row behind
        writeXRow ( count, x0, y0 - 1 );
        count += 5;

        for ( ;; )
        {
            writeXRow ( count, x0, y0 );

            y0++;
            count += 5;

            if ( y0 > ny_ )
            {
                return count;
            }

            if ( err > 0 )
            {
                x0--;
                err -= dy2;
                if ( x0 < -2 )
                {
                    return count;
                }
            }
            err += dx2;
        }
        return count;
    }

    inline int lineBackwardDown ( int x0, int y0, int x1, int y1 )
    {
        const int dx = x0 - x1;
        const int dy = y0 - y1;
        const int dx2 = 2 * dx;
        const int dy2 = 2 * dy;
        int err = dy2 - dx;

        int count = 0;
        //write one row behind
        writeYRow ( count, x0 + 1, y0 );
        count += 5;

        for ( ;; )
        {
            writeYRow ( count, x0, y0 );

            x0--;
            count += 5;

            if ( x0 < 0 )
            {
                return count;
            }

            if ( err > 0 )
            {
                y0--;
                err -= dx2;
                if ( y0 < -2 )
                {
                    return count;
                }
            }
            err += dy2;
        }
        return count;
    }

    inline int lineBackwardDown_T ( int x0, int y0, int x1, int y1 )
    {
        const int dx = x0 - x1;
        const int dy = y0 - y1;
        const int dx2 = 2 * dx;
        const int dy2 = 2 * dy;
        int err = dx2 - dy;

        int count = 0;
        //write one row behind
        writeXRow ( count, x0, y0 + 1 );
        count += 5;

        for ( ;; )
        {
            writeXRow ( count, x0, y0 );

            y0--;
            count += 5;

            if ( y0 < 0 )
            {
                return count;
            }

            if ( err > 0 )
            {
                x0--;
                err -= dy2;
                if ( x0 < -2 )
                {
                    return count;
                }
            }
            err += dx2;
        }
        return count;
    }

    //cell_nr = ix + iy * nx_;
    inline void writeYRow ( int count, int cell_x, int cell_y )
    {
        //std::cout << "count " << count << std::endl;
        cells_x_[count    ] = cell_x;
        cells_x_[count + 1] = cell_x;
        cells_x_[count + 2] = cell_x;
        cells_x_[count + 3] = cell_x;
        cells_x_[count + 4] = cell_x;

        cells_y_[count    ] = cell_y - 2;
        cells_y_[count + 1] = cell_y - 1;
        cells_y_[count + 2] = cell_y    ;
        cells_y_[count + 3] = cell_y + 1;
        cells_y_[count + 4] = cell_y + 2;
    }

    inline void writeXRow ( int count, int cell_x, int cell_y )
    {
        //std::cout << "count " << count << std::endl;
        cells_x_[count    ] = cell_x - 2;
        cells_x_[count + 1] = cell_x - 1;
        cells_x_[count + 2] = cell_x    ;
        cells_x_[count + 3] = cell_x + 1;
        cells_x_[count + 4] = cell_x + 2;

        cells_y_[count    ] = cell_y;
        cells_y_[count + 1] = cell_y;
        cells_y_[count + 2] = cell_y;
        cells_y_[count + 3] = cell_y;
        cells_y_[count + 4] = cell_y;
    }

    inline int writeNarrowYRowForward ( int cell_x, int cell_y )
    {
        int i = 0;
        cell_x--;
        while ( cell_x < nx_ )
        {
            cells_x_[i] = cell_x;
            cells_y_[i] = cell_y - 1;
            i++;
            cells_x_[i] = cell_x;
            cells_y_[i] = cell_y    ;
            i++;
            cells_x_[i] = cell_x;
            cells_y_[i] = cell_y + 1;
            i++;
            cell_x++;
        }
        return i;
    }

    inline int writeNarrowYRowBackward ( int cell_x, int cell_y )
    {
        int i = 0;
        cell_x ++;
        while ( cell_x >= 0 )
        {
            cells_x_[i] = cell_x;
            cells_y_[i] = cell_y - 1;
            i++;
            cells_x_[i] = cell_x;
            cells_y_[i] = cell_y    ;
            i++;
            cells_x_[i] = cell_x;
            cells_y_[i] = cell_y + 1;
            i++;
            cell_x--;
        }
        return i;
    }

    inline int writeNarrowXRowUpward ( int cell_x, int cell_y )
    {
        int i = 0;
        cell_y--;
        while ( cell_y < ny_ )
        {
            cells_x_[i] = cell_x - 1;
            cells_y_[i] = cell_y;
            i++;
            cells_x_[i] = cell_x;
            cells_y_[i] = cell_y;
            i++;
            cells_x_[i] = cell_x + 1;
            cells_y_[i] = cell_y;
            i++;
            cell_y++;
        }
        return i;
    }

    inline int writeNarrowXRowDownward ( int cell_x, int cell_y )
    {
        int i = 0;
        cell_y++;
        while ( cell_y >= 0 )
        {
            cells_x_[i] = cell_x - 1;
            cells_y_[i] = cell_y;
            i++;
            cells_x_[i] = cell_x;
            cells_y_[i] = cell_y;
            i++;
            cells_x_[i] = cell_x + 1;
            cells_y_[i] = cell_y;
            i++;
            cell_y--;
        }
        return i;
    }
};
};



