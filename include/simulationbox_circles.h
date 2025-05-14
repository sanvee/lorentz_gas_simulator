#pragma once

#include "container.h"
#include "obstacles.h"
#include "cell_list.h"
#include "globals.h"

#include <iostream>
#include <cmath>

namespace SimulationBox // sets up a Periodic box with ghost obstacles
{

class PeriodicBoxCircularObstacles
{
public:

    int nr_ghosts_approx_;
    int nr_total_approx_;
    Container *container_;
    Obstacles::Circles *obstacles_;
    CellList *celllist_;

public:

    ~PeriodicBoxCircularObstacles()
    {
        delete obstacles_;
        delete container_;

        if ( celllist_ != NULL )
        {
            delete celllist_;
        }
    };

    PeriodicBoxCircularObstacles()
    {
        nr_ghosts_approx_ = 0;
        nr_total_approx_ = 0;
        //obstacle_density_ = 0.0;

        container_ = NULL;
        obstacles_ = NULL;
        celllist_ = NULL;
    };
    /*
        PeriodicBoxBinaryCircularObstacles ( double lx, double ly, int n_obstacles, int n_ghosts, int nr_ghosts_approx, double r1, double r2, CellList::FloorList &celllist )
        {
            celllist_ = &celllist;
            // approximate upper limit of the number of ghost Obstacles)
            nr_ghosts_approx_ = nr_ghosts_approx;
            nr_total_approx_ = nr_ghosts_approx_ + n_obstacles;

            obstacles_ = new Obstacles::Circles ( n_obstacles, nr_ghosts_approx_ );
            container_ = new Container ( lx, ly );
            obstacles_->n_ghosts_ = n_ghosts;
            obstacles_->n_tot_ = n_obstacles + n_ghosts;
            obstacles_->setSizesBinary ( r1, r2 );
            //obstacle_density_ = obstacles_->total_surface_ / container_->surface_;
        }

        PeriodicBoxBinaryCircularObstacles ( double lx, double ly, double x, double y, double r )
        {
            // approximate upper limit of the number of ghost Obstacles)
            nr_ghosts_approx_ = 4;
            nr_total_approx_ = nr_ghosts_approx_ + 1;

            container_ = new Container ( lx, ly );
            obstacles_ = new Obstacles::Circles ( 1, nr_ghosts_approx_ );
            obstacles_->setOne ( x, y, r, 1 );
            //obstacle_density_ = obstacles_->total_surface_ / container_->surface_;
            setGhosts();
        }
    */
    PeriodicBoxCircularObstacles ( double lx, double ly, double target_cell_length, int n1, double r1, bool overlapp, bool fill_celllist )
    {
        int n_obstacles = n1;
        double max_radius = r1;

        // approximate upper limit of the number of ghost Obstacles)
        nr_ghosts_approx_ = std::ceil ( ( ( n_obstacles ) / ( lx * ly ) ) * ( 2 * max_radius * lx + 2 * max_radius * ly + 4 * max_radius * max_radius ) * 2 ) + 100;
        nr_total_approx_ = nr_ghosts_approx_ + n_obstacles;

        obstacles_ = new Obstacles::Circles ( n_obstacles, nr_ghosts_approx_ );
        obstacles_->setSizes (r1);

        container_ = new Container ( lx, ly );

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

    inline void setGhosts ( Obstacles::Circles *obstacles )
    {
        int n_obstacles = obstacles->n_obstacles_;
        double lx = container_->lx_;
        double ly = container_->ly_;
        double max_radius = obstacles->max_radius_;

        //adding the ghost atoms:
        int flag;
        int count = n_obstacles;

        for ( int i = 0; i < n_obstacles; i++ )
        {
            if ( count >= nr_total_approx_ - 3 )
            {
                std::cout << "WARNING: GHOST LIST FULL IN CLASS SimulationBox::PeriodicBoxBinaryObstacles" << std::endl;
                std::cout << "not all Ghost atoms might be listed" << std::endl;
                return;
            }

            flag = 0;

            if ( obstacles->x_[i] > lx - max_radius )   // right border
            {
                flag += 1;
            }

            if ( obstacles->x_[i] < max_radius )     // left border
            {
                flag -= 1;
            }

            if ( obstacles->y_[i] > ly - max_radius )   // upper border
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
                obstacles->setOne ( obstacles->x_[i] - lx, obstacles->y_[i], obstacles->radius_[i], count++ );
                break;
            }
            case -1:
            {
                obstacles->setOne ( obstacles->x_[i] + lx, obstacles->y_[i], obstacles->radius_[i], count++ );
                break;
            }

            case 3:
            {
                obstacles->setOne ( obstacles->x_[i], obstacles->y_[i] - ly, obstacles->radius_[i], count++ );
                break;
            }

            case -3:
            {
                obstacles->setOne ( obstacles->x_[i], obstacles->y_[i] + ly, obstacles->radius_[i], count++ );
                break;
            }

            case 2:
            {
                obstacles->setOne ( obstacles->x_[i] + lx, obstacles->y_[i] - ly, obstacles->radius_[i], count++ );
                obstacles->setOne ( obstacles->x_[i], obstacles->y_[i] - ly, obstacles->radius_[i], count++ );
                obstacles->setOne ( obstacles->x_[i] + lx, obstacles->y_[i], obstacles->radius_[i], count++ );
                break;
            }

            case 4:
            {
                obstacles->setOne ( obstacles->x_[i] - lx, obstacles->y_[i] - ly, obstacles->radius_[i], count++ );
                obstacles->setOne ( obstacles->x_[i], obstacles->y_[i] - ly, obstacles->radius_[i], count++ );
                obstacles->setOne ( obstacles->x_[i] - lx, obstacles->y_[i], obstacles->radius_[i], count++ );
                break;
            }

            case -2:
            {
                obstacles->setOne ( obstacles->x_[i] - lx, obstacles->y_[i] + ly, obstacles->radius_[i], count++ );
                obstacles->setOne ( obstacles->x_[i], obstacles->y_[i] + ly, obstacles->radius_[i], count++ );
                obstacles->setOne ( obstacles->x_[i] - lx, obstacles->y_[i], obstacles->radius_[i], count++ );
                break;
            }

            case -4:
            {
                obstacles->setOne ( obstacles->x_[i] + lx, obstacles->y_[i] + ly, obstacles->radius_[i], count++ );
                obstacles->setOne ( obstacles->x_[i], obstacles->y_[i] + ly, obstacles->radius_[i], count++ );
                obstacles->setOne ( obstacles->x_[i] + lx, obstacles->y_[i], obstacles->radius_[i], count++ );
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
        obstaclesFile.precision ( 8 );
        obstaclesFile << obstacles_->type_ << std::endl;
        obstaclesFile << obstacles_->n_tot_ << std::endl;
        obstaclesFile << 1 << std::endl;
        for ( int i = 0; i < obstacles_->n_tot_; ++i )
        {
            obstaclesFile << i << ' ' << obstacles_->x_[i] << ' ' << obstacles_->y_[i] << ' ' << obstacles_->radius_[i] << std::endl;
        }
        obstaclesFile.close();
    }
};
};
