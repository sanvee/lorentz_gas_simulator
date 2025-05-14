#pragma once

#include "cell_list.h"
#include "stencils.h"

#include <stdint.h>


class Testing;
namespace Neighbourlist
{
class Simple
{
    //takes care of
    //-size of the cells in the box
    //-adjust cutoff distance
    //-creation of the stencil
    //-returns the obstacles to be considered for collision search

public:

    int *neighbour_obstacles_;

    //private:
    CellList *celllist_;
    Stencil::FullPeriodicTrimmedStencil *stencil_;
    Stencil::FullPeriodicTrimmedStencil *placement_stencil_;

    double cutoff_;
    double target_cell_length_;

    int celllist_nx_;
    int celllist_ny_;
    int last_cell_x_;
    int last_cell_y_;
    int last_list_offset_;
    int nlist_length_;

public:

    Simple()
    {
        neighbour_obstacles_ = NULL;
        celllist_ = NULL;
        stencil_ = NULL;
        placement_stencil_ = NULL;
    }

    ~Simple()
    {
        delete[] neighbour_obstacles_;
    }

    Simple ( int lenght, CellList &celllist, Stencil::FullPeriodicTrimmedStencil &stencil, Stencil::FullPeriodicTrimmedStencil &placement_stencil )
    {
        nlist_length_ = lenght;
        celllist_ = &celllist;
        stencil_ = &stencil;
        placement_stencil_ = &placement_stencil;

        neighbour_obstacles_ = new int[nlist_length_];

        celllist_nx_ = celllist_->nx_;
        celllist_ny_ = celllist_->ny_;

        last_cell_x_ = -1;
        last_cell_y_ = -1;
        last_list_offset_ = -1;

        double mem_size = sizeof ( int ) * nlist_length_;
        //GlobalCounters::memsize += mem_size;
        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Neighbourlist (simple)---" << std::endl;
            printf ( "neighbourlist lenght %i obstacles \n", nlist_length_ );
            printf ( "total heapsize of the neighbourlist is %.2f kB \n", mem_size / 1000.0 );
        }
    }

public:

    inline int getNeighbours ( double rx, double ry )
    {
        int curr_cell_x = std::floor ( rx / celllist_->l_ );
        int curr_cell_y = std::floor ( ry / celllist_->l_ );

        //Shrink to cellist dimension
        curr_cell_x = curr_cell_x % celllist_nx_;
        curr_cell_y = curr_cell_y % celllist_ny_;

        //remap it to real cell
        celllist_->remapCells ( curr_cell_x, curr_cell_y );

        //no need to change anything if we have not moved out of the last cell.
        if ( curr_cell_x == last_cell_x_ and curr_cell_y == last_cell_y_ )
        {
            return last_list_offset_;
        }

        last_cell_x_ = curr_cell_x;
        last_cell_y_ = curr_cell_y;

        int neighbour_cell_x, neighbour_cell_y;
        int list_offset = 0;

        for ( auto i = 0; i < stencil_->n_stencil_; ++i )
        {
            neighbour_cell_x = curr_cell_x + stencil_->cell_offsets_x_[i];
            neighbour_cell_y = curr_cell_y + stencil_->cell_offsets_y_[i];

            celllist_->remapCells ( neighbour_cell_x, neighbour_cell_y );

            int obs_num = celllist_->head_cells_[neighbour_cell_x + neighbour_cell_y * celllist_nx_];

            while ( obs_num >= 0 )
            {
                neighbour_obstacles_[list_offset] = obs_num;
                obs_num = celllist_->list_cells_[obs_num];
                ++list_offset;
            }
        }

        last_list_offset_ = list_offset;
        return list_offset;
    }

    int getPossibleCoveringNeighbours ( double x, double y )
    {
        int curr_cell_x = std::floor ( x / celllist_->l_ );
        int curr_cell_y = std::floor ( y / celllist_->l_ );

        int neighbour_cell_x, neighbour_cell_y;
        int list_offset = 0;

        for ( auto i = 0; i < placement_stencil_->n_stencil_; ++i )
        {
            neighbour_cell_x = curr_cell_x + placement_stencil_->cell_offsets_x_[i];
            neighbour_cell_y = curr_cell_y + placement_stencil_->cell_offsets_y_[i];

            celllist_->remapCells ( neighbour_cell_x, neighbour_cell_y );

            int obs_num = celllist_->head_cells_[neighbour_cell_x + neighbour_cell_y * celllist_nx_];

            while ( obs_num >= 0 )
            {
                neighbour_obstacles_[list_offset] = obs_num;
                obs_num = celllist_->list_cells_[obs_num];
                ++list_offset;
            }
        }
        return list_offset;
    }
};
};


