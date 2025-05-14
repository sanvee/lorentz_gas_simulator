#pragma once

#include "container.h"
#include "obstacles.h"
#include "globals.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

namespace SimulationBox // sets up a Periodic box with ghost obstacles
{

class PeriodicBoxCrossObstacles
{
public:

    int nr_ghosts_approx_;
    int nr_total_approx_;
    Container *container_;
    Obstacles::Crosses *obstacles_;
    CellList *celllist_;

public:

    ~PeriodicBoxCrossObstacles()
    {
        delete obstacles_;
        delete container_;
        if ( celllist_ != NULL )
        {
            delete celllist_;
        }
    };

    PeriodicBoxCrossObstacles()
    {
        nr_ghosts_approx_ = 0;
        nr_total_approx_ = 0;

        container_ = NULL;
        obstacles_ = NULL;
        celllist_  = NULL;
    };

    PeriodicBoxCrossObstacles ( double lx, double ly, double target_cell_length,  int n_obstacles, double l, double d, double alpha, bool overlapp, bool fill_celllist )
    {
        double max_radius = std::sqrt ( 0.25 * l * l + 0.25 * d * d );

        // approximate upper limit of the number of ghost Obstacles)
        nr_ghosts_approx_ = std::ceil ( ( ( n_obstacles ) / ( lx * ly ) ) * ( 2 * max_radius * lx + 2 * max_radius * ly + 4 * max_radius * max_radius ) * 2 ) + 100;
        nr_total_approx_ = nr_ghosts_approx_ + n_obstacles;

        obstacles_ = new Obstacles::Crosses ( l, d, alpha, n_obstacles, nr_ghosts_approx_ );

        if ( GlobalFlags::verbose )
        {
            if ( overlapp )
            {
                std::cout << "Overlapping Obstacles" << std::endl;
            }
            else
            {
                std::cout << "NON Overlapping Obstacles" << std::endl;
            }
        }

        container_ = new Container ( lx, ly );

        if ( overlapp )
        {
            container_->setRandomPositions ( *obstacles_ );
        }
        else
        {
            container_->setRandomPositionsNoOverlapp ( *obstacles_ );
        }

        setGhosts ( obstacles_ );

        if ( fill_celllist )
        {
            celllist_ = new CellList ( obstacles_->n_tot_, container_->lx_, container_->ly_, target_cell_length );
            celllist_->fill ( *obstacles_ );
        }
        else
        {
            celllist_ = NULL;
        }
    }

    inline void setGhosts ( Obstacles::Crosses *obstacles )
    {
        int n_obstacles = obstacles->n_obstacles_;
        double lx = container_->lx_;
        double ly = container_->ly_;
        double max_radius = obstacles->max_radius_;

        //adding the ghost atoms:
        int flag;
        int count = n_obstacles;
        double alpha = 0;
        for ( int i = 0; i < n_obstacles; i++ )
        {
            if ( obstacles->angle_ < 0 )
            {
                alpha = obstacles->angles_[i];
            }
            else
            {
                alpha = obstacles->angle_;
            }

            if ( count >= nr_total_approx_ - 3 )
            {
                std::cout << "WARNING: GHOST LIST FULL IN CLASS SimulationBox::PeriodicBoxBinaryObstacles" << std::endl;
                std::cout << "not all Ghost atoms might be listed" << std::endl;
                return;
            }

            flag = 0;

            if ( obstacles->x_[i] > lx - max_radius ) // right border
            {
                flag += 1;
            }

            if ( obstacles->x_[i] < max_radius )     // left border
            {
                flag -= 1;
            }

            if ( obstacles->y_[i] > ly - max_radius ) // upper border
            {
                flag += 3;
            }
            if ( obstacles->y_[i] < max_radius )     // lowerborder
            {
                flag -= 3;
            }

            switch ( flag )
            {
            case 1:
            {
                obstacles->set ( count++, alpha, obstacles->x_[i] - lx, obstacles->y_[i] );
                break;
            }
            case -1:
            {
                obstacles->set ( count++, alpha, obstacles->x_[i] + lx, obstacles->y_[i] );
                break;
            }

            case 3:
            {
                obstacles->set ( count++, alpha, obstacles->x_[i], obstacles->y_[i] - ly );
                break;
            }
            case -3:
            {
                obstacles->set ( count++, alpha, obstacles->x_[i], obstacles->y_[i] + ly );
                break;
            }

            case 2:
            {
                obstacles->set ( count++, alpha, obstacles->x_[i] + lx, obstacles->y_[i] - ly );
                obstacles->set ( count++, alpha, obstacles->x_[i], obstacles->y_[i] - ly );
                obstacles->set ( count++, alpha, obstacles->x_[i] + lx, obstacles->y_[i] );
                break;
            }

            case 4:
            {
                obstacles->set ( count++, alpha, obstacles->x_[i] - lx, obstacles->y_[i] - ly );
                obstacles->set ( count++, alpha, obstacles->x_[i], obstacles->y_[i] - ly );
                obstacles->set ( count++, alpha, obstacles->x_[i] - lx, obstacles->y_[i] );
                break;
            }

            case -2:
            {
                obstacles->set ( count++, alpha, obstacles->x_[i] - lx, obstacles->y_[i] + ly );
                obstacles->set ( count++, alpha, obstacles->x_[i], obstacles->y_[i] + ly );
                obstacles->set ( count++, alpha, obstacles->x_[i] - lx, obstacles->y_[i] );
                break;

            }

            case -4:
            {
                obstacles->set ( count++, alpha, obstacles->x_[i] + lx, obstacles->y_[i] + ly );
                obstacles->set ( count++, alpha, obstacles->x_[i], obstacles->y_[i] + ly );
                obstacles->set ( count++, alpha, obstacles->x_[i] + lx, obstacles->y_[i] );
                break;
            }
            }
        }
        // Updating the number of Obstacles:
        obstacles->n_tot_ = count;
        obstacles->n_ghosts_ = count - n_obstacles;
    }

    void writeObstacles ( std::string filename )
    {
        std::ofstream obstaclesFile;
        obstaclesFile.open ( filename );
        obstaclesFile.precision ( 5 );
        obstaclesFile << 12 << std::endl;
        for ( int i = 0; i < obstacles_->n_tot_; ++i )
        {
            for ( int j = 0; j < 12; j++ )
            {
                obstaclesFile << i << ' ' << obstacles_->px_[13*i+j] << ' ' << obstacles_->py_[13*i+j] << std::endl;
            }
        }
        obstaclesFile.close();
    }
};
};
