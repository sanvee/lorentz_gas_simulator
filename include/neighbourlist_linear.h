#pragma once

#include "trajectory.h"
#include "cell_list.h"
#include "stencils.h"

#include <stdint.h>

class Testing;
namespace Neighbourlist
{

class Linear
{
    //takes care of
    //-size of the cells in the box
    //-adjust cutoff distance
    //-creation of the stencil
    //-returns the obstacles to be considered for collision search

public:

    int *neighbour_obstacles_;
    int nlist_length_;
    int celllist_nx_;
    int celllist_ny_;

    double shift_len_;

    double cell_l_;

    //private:
    CellList *celllist_;
    Stencil::Linear *stencil_;
    Stencil::FullPeriodicTrimmedStencil *placement_stencil_;

public:

    Linear()
    {
        neighbour_obstacles_ = NULL;
        celllist_ = NULL;
        stencil_ = NULL;
        placement_stencil_ = NULL;
    }

    ~Linear()
    {
        delete[] neighbour_obstacles_;
    }

    Linear ( int lenght, CellList &celllist, Stencil::Linear &stencil, Stencil::FullPeriodicTrimmedStencil &placement_stencil )
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

        shift_len_ = 2 * ( ( celllist_->nx_ + celllist_->ny_ ) * celllist_->l_ );

        double mem_size = sizeof ( int ) * nlist_length_;

        //GlobalCounters::memsize += mem_size;
        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Neighbourlist (linear)---" << std::endl;
            printf ( "neighbourlist lenght %d obstacles \n", nlist_length_ );
            printf ( "total heapsize of the neighbourlist is %.2f kB \n", mem_size / 1000.0 );
        }
    }
    /*
        template<class T_SimulationBox>
        Linear ( T_SimulationBox &simbox, bool windtree )
        {
            double target_cell_length = 0;
            if ( windtree )
            {
                target_cell_length = simbox.obstacles_->max_radius_ * 2.001;
            }
            else
            {
                target_cell_length = simbox.obstacles_->max_radius_ * 2.901;
            }

            celllist_ = simbox.celllist_;
            max_cells_ = celllist_->n_cells_;

            stencil_ = new Stencil::Linear ( *celllist_ );

            nlist_length_ = std::ceil ( simbox.obstacles_->n_tot_ * celllist_->l_ * celllist_->l_ / ( simbox.container_->lx_ * simbox.container_->ly_ ) * 2 ) * stencil_->max_lenght_;

            neighbour_obstacles_ = new int [nlist_length_];

            placement_stencil_ = new Stencil::FullPeriodicTrimmedStencil ( *celllist_, simbox.obstacles_->max_radius_ );
            double mem_size = sizeof ( int ) * nlist_length_;

            shift_len_ = 2 * ( simbox.container_->lx_ + simbox.container_->ly_ );

            cell_l_ = celllist_->l_;

            celllist_nx_ = celllist_->nx_;
            celllist_ny_ = celllist_->ny_;

            //GlobalCounters::memsize += mem_size;
            if ( GlobalFlags::verbose )
            {
                std::cout << std::endl << "---Neighbourlist (linear)---" << std::endl;
                printf ( "neighbourlist lenght %d obstacles \n", nlist_length_ );
                printf ( "total heapsize of the neighbourlist is %.2f kB \n", mem_size / 1000.0 );
            }
        }
    */
public:



    inline int getNeighbours ( int cell_x1, int cell_y1, double vx, double vy )
    {
        int cell_x2, cell_y2;

        celllist_->getCell ( ( cell_x1 + 0.5 ) * cell_l_ + shift_len_ * vx, ( cell_y1 + 0.5 ) * cell_l_ + shift_len_ * vy, cell_x2, cell_y2 );
        int ncells = stencil_->getCells ( cell_x1, cell_y1, cell_x2, cell_y2, vx, vy );
        int list_offset = 0;

        for ( auto i = 0; i < ncells; i++ )
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

    inline int getNeighboursInStencilCell ( int nr )
    {
        int list_offset = 0;
        auto cell_x = stencil_->cells_x_[nr];
        if ( cell_x < 0 or cell_x >= celllist_nx_ )
        {
            return 0;
        }
        auto cell_y = stencil_->cells_y_[nr];
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


    inline int getNeighbourCells ( int cell_x1, int cell_y1, double vx, double vy )
    {
        int cell_x2, cell_y2;
        celllist_->getCell ( ( cell_x1 + 0.5 ) * cell_l_ + shift_len_ * vx, ( cell_y1 + 0.5 ) * cell_l_ + shift_len_ * vy, cell_x2, cell_y2 );

        return stencil_->getCells ( cell_x1, cell_y1, cell_x2, cell_y2, vx, vy );
    }


    inline int getNeighbourCellsWindtree ( int cell_x1, int cell_y1, double vx, double vy )
    {
        return stencil_->getCellsWindtree ( cell_x1, cell_y1, vx, vy );
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

