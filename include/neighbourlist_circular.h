#pragma once

#include "trajectory.h"
#include "simulationbox.h"
#include "cell_list.h"
#include "stencils.h"

#include <stdint.h>
#include <cmath>

class Testing;

namespace Neighbourlist
{

class Circular
{
    //takes care of
    //-size of the cells in the box
    //-adjust cutoff distance
    //-creation of the stencil
    //-returns the obstacles to be considered for collision search
    //-implicit cast of cell number from int64_t in stencil to int in n_list !!! There should not be overflow!

public:

    int *neighbour_obstacles_;
    int nlist_length_;
    int celllist_nx_;
    int celllist_ny_;

    int ncells_;
    double cell_l_;
    int max_cells_;

    //private:
    CellList *celllist_;
    Stencil::Circular *stencil_;
    Stencil::FullPeriodicTrimmedStencil *placement_stencil_;

    //-----------------------------
    //int debug_counter = 0;
    //std::ofstream obstaclesFile;
    //-----------------------------

public:

    Circular()
    {
        neighbour_obstacles_ = NULL;
        celllist_ = NULL;
        stencil_ = NULL;
        placement_stencil_ = NULL;
    }

    ~Circular()
    {
        delete[] neighbour_obstacles_;
    }

    Circular ( int lenght, CellList &celllist, Stencil::Circular &stencil, Stencil::FullPeriodicTrimmedStencil &placement_stencil )
    {
        nlist_length_ = lenght;
        neighbour_obstacles_ = new int [nlist_length_];

        //Cell list stuff
        celllist_ = &celllist;
        cell_l_ = celllist_->l_;
        celllist_nx_ = celllist_->nx_;
        celllist_ny_ = celllist_->ny_;

        //Stencils
        stencil_ = &stencil;
        placement_stencil_ = &placement_stencil;

        double mem_size = sizeof ( int ) * nlist_length_;

        //GlobalCounters::memsize += mem_size;
        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Neighbourlist (circular)---" << std::endl;
            printf ( "neighbourlist lenght %d obstacles \n", nlist_length_ );
            printf ( "total heapsize of the neighbourlist is %.2f kB \n", mem_size / 1000.0 );
        }
    }

public:

    inline int getNeighbours ( double cx, double cy, double ex, double ey, double rx, double ry )
    {
        int ncells = stencil_->getCells ( cx, cy, ex, ey, rx, ry );
        //std::cout << " " << ncells << " " << std::endl;
        int list_offset = 0;
        for ( int i = 0; i < ncells; ++i )
        {
            auto cell_x = stencil_->cells_x_[i];
            if ( cell_x < 0 or cell_x >= celllist_nx_ )
            {
                continue;
            }
            auto cell_y = stencil_->cells_y_[i];
            if ( cell_y < 0 or cell_y >= celllist_ny_ )
            {
                continue;
            }

            int obs_num = celllist_->head_cells_[cell_x + cell_y * celllist_nx_];

            while ( obs_num >= 0 )
            {
                neighbour_obstacles_[list_offset] = obs_num;
                obs_num = celllist_->list_cells_[obs_num];
                ++list_offset;
            }
        }
        return list_offset;
    }

    inline int getNeighbours ( double cx, double cy, double rx, double ry )
    {

        int ncells = stencil_->getCells ( cx, cy, rx, ry );
        //std::cout << " " << ncells << " " << std::endl;
        //      ncells_ = ncells;
        int list_offset = 0;
        for ( int i = 0; i < ncells; ++i )
        {
            int cell_x = stencil_->cells_x_[i];
            if ( cell_x < 0 or cell_x >= celllist_nx_ )
            {
                continue;
            }
            int cell_y = stencil_->cells_y_[i];
            if ( cell_y < 0 or cell_y >= celllist_ny_ )
            {
                continue;
            }

            int obs_num = celllist_->head_cells_[cell_x + cell_y * celllist_nx_];

            while ( obs_num >= 0 )
            {
                neighbour_obstacles_[list_offset] = obs_num;
                obs_num = celllist_->list_cells_[obs_num];
                ++list_offset;
            }
        }
        return list_offset;
    }

    inline int getNeighbourCells ( double cx, double cy, double ex, double ey, double rx, double ry )
    {
        return stencil_->getCells ( cx, cy, ex, ey, rx, ry );
    }

    inline int getNeighbourCells ( double cx, double cy, double rx, double ry )
    {
        //int ncells = stencil_->getCells ( cx, cy, rx, ry );
        //std::cout << ncells << std::endl;
        //return ncells;
        return stencil_->getCells ( cx, cy, rx, ry );
    }

    int getNeighboursInStencilCell ( int curr_cell )
    {
        int list_offset = 0;
        int cell_x = stencil_->cells_x_[curr_cell];
        if ( cell_x < 0 or cell_x >= celllist_nx_ )
        {
            return 0;
        }
        int cell_y = stencil_->cells_y_[curr_cell];
        if ( cell_y < 0 or cell_y >= celllist_ny_ )
        {
            return 0;
        }

        int obs_num = celllist_->head_cells_[cell_x + cell_y * celllist_nx_];

        while ( obs_num >= 0 )
        {
            neighbour_obstacles_[list_offset] = obs_num;
            obs_num = celllist_->list_cells_[obs_num];
            ++list_offset;
        }

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

            int obs_num = celllist_->head_cells_[ neighbour_cell_x + neighbour_cell_y * celllist_->nx_];

            while ( obs_num >= 0 )
            {
                neighbour_obstacles_[list_offset] = obs_num;
                obs_num = celllist_->list_cells_[obs_num];
                ++list_offset;
            }
        }
        return list_offset;
    }



    /*
        inline int getNeighbours ( double cx, double cy, double rx, double ry ) //, Obstacles::Squares &obstacles )
        {
            int ncells = stencil_->getCells ( cx, cy, rx, ry );
            ncells_ = ncells;
            int list_offset = 0;
            for ( int i = 0; i < ncells; ++i )
            {
                int64_t cell_x = stencil_->cells_x_[i];
                if ( cell_x < 0 or cell_x >= celllist_nx_ )
                {
                    continue;
                }
                int64_t cell_y = stencil_->cells_y_[i];
                if ( cell_y < 0 or cell_y >= celllist_ny_ )
                {
                    continue;
                }

                int64_t cell_nr = cell_x + cell_y * celllist_nx_;
                int obs_nr = celllist_->head_cells_[cell_nr];

                while ( obs_nr >= 0 )
                {
                    neighbour_obstacles_[list_offset] = obs_nr;
                    obs_nr = celllist_->list_cells_[obs_nr];
                    ++list_offset;
                }
            }
            return list_offset;
        }

        inline int getNeighbourCells ( double cx, double cy, double ex, double ey, double rx, double ry )
        {
            int ncells = stencil_->getCells ( cx, cy, ex, ey, rx, ry );
            return ncells;
        }

        inline int getNeighbourCells ( double cx, double cy, double rx, double ry )
        {
            int ncells = stencil_->getCells ( cx, cy, rx, ry );
            return ncells;
        }

        inline int getNeighboursInStencilCell ( int nr )
        {
            int list_offset = 0;
            int64_t cell_x = stencil_->cells_x_[nr];
            if ( cell_x < 0 or cell_x >= celllist_nx_ )
            {
                return 0;
            }
            int64_t cell_y = stencil_->cells_y_[nr];
            if ( cell_y < 0 or cell_y >= celllist_ny_ )
            {
                return 0;
            }

            int64_t cell_nr = cell_x + cell_y * celllist_nx_;
            int obs_nr = celllist_->head_cells_[cell_nr];

            while ( obs_nr >= 0 )
            {
                neighbour_obstacles_[list_offset] = obs_nr;
                obs_nr = celllist_->list_cells_[obs_nr];
                ++list_offset;
            }
            return list_offset;
        }

        int getPossibleCoveringNeighbours ( double x, double y )
        {
            int curr_cell_x = std::floor ( x / celllist_->l_ );
            int curr_cell_y = std::floor ( y / celllist_->l_ );

            int neighbour_cell_x, neighbour_cell_y;
            int obs_num;
            int list_offset = 0;
            for ( auto i = 0; i < placement_stencil_->n_stencil_; ++i )
            {
                neighbour_cell_x = curr_cell_x + placement_stencil_->cell_offsets_x_[i];
                neighbour_cell_y = curr_cell_y + placement_stencil_->cell_offsets_y_[i];

                celllist_->remapCells ( neighbour_cell_x, neighbour_cell_y );

                obs_num = celllist_->head_cells_[ neighbour_cell_x + neighbour_cell_y * celllist_->nx_];

                while ( obs_num >= 0 )
                {
                    neighbour_obstacles_[list_offset] = obs_num;
                    obs_num = celllist_->list_cells_[obs_num];
                    ++list_offset;
                }
            }
            return list_offset;
        }
    */


};
};


