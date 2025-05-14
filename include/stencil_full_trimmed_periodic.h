#pragma once

#include <iostream>
#include <cmath>
#include <vector> 
#include "globals.h"
#include "cell_list.h"


namespace Stencil
{

class FullPeriodicTrimmedStencil
{
private:
    // dimensions of the stencil in both dimensions ( from sxl_ to sxh )
    int s_h_;
    int s_l_;

    int s_max_; //holds the maximun of the possible number of cells in the stencil.
public:
    int n_stencil_;

    double cutoff_;
    double cutoff_sq_;

    std::vector<int> cell_offsets_x_;
    std::vector<int> cell_offsets_y_;

    FullPeriodicTrimmedStencil ( CellList &celllist, double cutoff )
    {
        cutoff_ = cutoff;
        cutoff_sq_ = cutoff * cutoff;
        int nx = celllist.nx_;
        int ny = celllist.ny_;
        

        //maximal Stencil size per dimension.

        s_h_ = std::ceil( cutoff_ * celllist.l_inv_ );
        
        s_l_ = -s_h_;


        // we need to trimm the stencil so it does not overlapp itself
        // trimm in both direction if nx or ny is even then we ignore the lower cells.

        if ( 2 * s_h_ + 1 >= nx )
        {
            s_h_ = nx / 2 ;
            s_l_ = -nx / 2 + ( nx + 1 ) % 2;
        }

        if ( 2 * s_h_ + 1 >= ny )
        {
            s_h_ = ny / 2 ;
            s_l_ = -ny / 2 + ( ny + 1 ) % 2;
        }

        s_max_ = ( s_h_ - s_l_ + 1 ) * ( s_h_ - s_l_ + 1 );

        cell_offsets_x_.reserve(s_max_);
        cell_offsets_y_.reserve(s_max_);
        
        n_stencil_ = 0;

        for ( int i = s_l_; i <= s_h_; ++i )
        {
            for ( int j = s_l_; j <= s_h_; ++j )
            {
                if ( celllist.cellDistanceSquared ( i,j ) < cutoff_sq_ )
                {
                    cell_offsets_x_.push_back(i);
                    cell_offsets_y_.push_back(j);
                }
            }
        }
        
        cell_offsets_x_.shrink_to_fit();
        cell_offsets_y_.shrink_to_fit();

        n_stencil_ = cell_offsets_x_.size();
        
        double mem_size = sizeof ( int ) * n_stencil_ * 2.0;

        if ( GlobalFlags::verbose )
        {
          std::cout << std::endl << "---Stencil(simple)---" << std::endl;
          printf ( "stencil has %d cells \n", n_stencil_);
          printf ( "total heapsize of the stencil is %.2f kB \n", mem_size / 1000.0 );
        }
    };
};
};
