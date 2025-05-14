#pragma once

#include "container.h"
#include "config.h"
#include "obstacles.h"
#include "globals.h"
#include "cell_list.h"
#include "mc_obstacles.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sys/stat.h>

namespace SimulationBox // sets up a Periodic box with ghost obstacles
{
class PeriodicBoxSquareObstacles
{
public:

    int nr_ghosts_approx_;
    int nr_total_approx_;

    Container *container_;
    Obstacles::Squares *obstacles_;
    CellList *celllist_;

public:

    ~PeriodicBoxSquareObstacles()
    {
        delete obstacles_;
        delete container_;
        
        if ( celllist_ != NULL )
        {
            delete celllist_;
        }
    };

    PeriodicBoxSquareObstacles()
    {
        nr_ghosts_approx_ = 0;
        nr_total_approx_ = 0;

        container_ = NULL;
        obstacles_ = NULL;
        celllist_ = NULL;
    };

    PeriodicBoxSquareObstacles ( double lx, double ly, double target_cell_length, int n_obstacles, double l, double alpha, bool overlapp, bool mc, int sweeps, bool fill_celllist )
    {
        double max_radius = l / std::sqrt ( 2. );

        // approximate upper limit of the number of ghost Obstacles)
        nr_ghosts_approx_ = std::ceil ( ( ( n_obstacles ) / ( lx * ly ) ) * ( 2 * max_radius * lx + 2 * max_radius * ly + 4 * max_radius * max_radius ) * 2 ) + 100;
        nr_total_approx_ = nr_ghosts_approx_ + n_obstacles;

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
            obstacles_ = new Obstacles::Squares ( l, alpha, n_obstacles, nr_ghosts_approx_ );
            container_ = new Container ( lx, ly );
            container_->setRandomPositions ( *obstacles_ );
        }
        else
        {
            if ( mc )
            {
                double n_density = n_obstacles / (lx*ly);
                //make square container
                int sqrt_n = std::ceil ( sqrt ( n_obstacles ) );
                int N = sqrt_n*sqrt_n;
                double L =  std::sqrt ( N / n_density );
                obstacles_ = new Obstacles::Squares ( l, 0, N, nr_ghosts_approx_ );
                container_ = new Container ( L, L );

                if ( GlobalFlags::verbose )
                {
                    std::cout << "old  lx ly  " << lx << ' ' << ly << std::endl; 
                    std::cout << "rearanged L " << L  << std::endl;
                    std::cout << " old N      " << n_obstacles << std::endl;
                    std::cout << "rearanged N " << N << std::endl;
                    std::cout << "rearanged density " << N / ( L*L ) * l * l << std::endl;
                    std::cout << "doing " << sweeps << " mc sweeps" << std::endl;
                }

                container_->setRandomPositionsNoOverlappMontecarlo ( *obstacles_, sqrt_n );

                double rate = Montecarlo::canonicalObstacles ( *container_, *obstacles_, sweeps );
                if ( GlobalFlags::verbose )
                {
                    std::cout << "acc rate " << rate << std::endl;
                }

            }
            else
            {
                obstacles_ = new Obstacles::Squares ( l, alpha, n_obstacles, nr_ghosts_approx_ );
                container_ = new Container ( lx, ly );
                container_->setRandomPositionsNoOverlapp ( *obstacles_ );
            }
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

    inline void setGhosts ( Obstacles::Squares *obstacles )
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
            //obstacles->obs_ID[i] = i;
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
                obstacles->set ( count, alpha, obstacles->x_[i] - lx, obstacles->y_[i] );
                //obstacles->obs_ID[count] = i;
                count++;
                break;
            }
            case -1:
            {
                obstacles->set ( count, alpha, obstacles->x_[i] + lx, obstacles->y_[i] );
                //obstacles->obs_ID[count] = i;
                count++;
                break;
            }

            case 3:
            {
                obstacles->set ( count, alpha, obstacles->x_[i], obstacles->y_[i] - ly );
                //obstacles->obs_ID[count] = i;
                count++;
                break;
            }
            case -3:
            {
                obstacles->set ( count, alpha, obstacles->x_[i], obstacles->y_[i] + ly );
                //obstacles->obs_ID[count] = i;
                count++;
                break;
            }

            case 2:
            {
                obstacles->set ( count, alpha, obstacles->x_[i] + lx, obstacles->y_[i] - ly );
                //obstacles->obs_ID[count] = i;
                count++;
                obstacles->set ( count, alpha, obstacles->x_[i], obstacles->y_[i] - ly );
                //obstacles->obs_ID[count] = i;
                count++;
                obstacles->set ( count, alpha, obstacles->x_[i] + lx, obstacles->y_[i] );
                //obstacles->obs_ID[count] = i;
                count++;
                break;
            }

            case 4:
            {
                obstacles->set ( count, alpha, obstacles->x_[i] - lx, obstacles->y_[i] - ly );
                //obstacles->obs_ID[count] = i;
                count++;
                obstacles->set ( count, alpha, obstacles->x_[i], obstacles->y_[i] - ly );
                //obstacles->obs_ID[count] = i;
                count++;
                obstacles->set ( count, alpha, obstacles->x_[i] - lx, obstacles->y_[i] );
                //obstacles->obs_ID[count] = i;
                count++;
                break;
            }

            case -2:
            {
                obstacles->set ( count, alpha, obstacles->x_[i] - lx, obstacles->y_[i] + ly );
                //obstacles->obs_ID[count] = i;
                count++;
                obstacles->set ( count, alpha, obstacles->x_[i], obstacles->y_[i] + ly );
                //obstacles->obs_ID[count] = i;
                count++;
                obstacles->set ( count, alpha, obstacles->x_[i] - lx, obstacles->y_[i] );
                //obstacles->obs_ID[count] = i;
                count++;
                break;

            }

            case -4:
            {
                obstacles->set ( count, alpha, obstacles->x_[i] + lx, obstacles->y_[i] + ly );
                //obstacles->obs_ID[count] = i;
                count++;
                obstacles->set ( count, alpha, obstacles->x_[i], obstacles->y_[i] + ly );
                //obstacles->obs_ID[count] = i;
                count++;
                obstacles->set ( count, alpha, obstacles->x_[i] + lx, obstacles->y_[i] );
                //obstacles->obs_ID[count] = i;
                count++;
                break;
            }
            }
        }
        // Updating the number of Obstacles:
        obstacles->n_tot_ = count;
        obstacles->n_ghosts_ = count - n_obstacles;
    }

    void writeObstacles_config ( std::string folder, std::string filename )
    {
        mkdir ( folder.c_str(), 0755 );
        std::string sub_folder = folder + "/configurations";
        mkdir ( sub_folder.c_str(), 0755 );
        std::ofstream obstaclesFile;
        obstaclesFile.open ( sub_folder + '/' + filename );
        obstaclesFile.precision ( std::numeric_limits<long double>::digits10 + 1 );

        obstaclesFile << obstacles_->n_obstacles_ << std::endl;
        obstaclesFile << container_->lx_ << ' ' << container_->ly_ << std::endl;
        obstaclesFile << obstacles_->l_ << std::endl;

        for ( int i = 0; i < obstacles_->n_obstacles_; ++i )
        {
            obstaclesFile << obstacles_->x_[i] << ' ' << obstacles_->y_[i] << std::endl;
        }

        for ( int i = 0; i < obstacles_->n_obstacles_; ++i )
        {
            for ( int j = 0; j < 4; j++ )
            {
                obstaclesFile << i << ' ' << obstacles_->px_[5 * i + j] << ' ' << obstacles_->py_[5 * i + j] << std::endl;
            }
        }
        obstaclesFile.close();
    }

    void writeObstacles ( std::string filename )
    {
        std::ofstream obstaclesFile;
        obstaclesFile.open ( filename );
        obstaclesFile.precision ( 8 );
        obstaclesFile << obstacles_->type_ << std::endl;
        obstaclesFile << obstacles_->n_tot_ << std::endl;
        obstaclesFile << 4 << std::endl;
        
        for ( int i = 0; i < obstacles_->n_tot_; ++i )
        {
            for ( int j = 0; j < 4; j++ )
            {
                obstaclesFile << i << ' ' << obstacles_->px_[5*i+j] << ' ' << obstacles_->py_[5*i+j] << std::endl;
            }
        }
        obstaclesFile.close();
    }
};
};
