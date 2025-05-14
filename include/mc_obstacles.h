#pragma once

#include "obstacles.h"
#include "container.h"
#include "cell_list.h"

namespace Montecarlo
{
    double canonicalObstacles ( Container &container, Obstacles::Squares &obstacles, int sweeps )
    {
        std::random_device rd;
        std::mt19937_64 engine ( rd() );
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        double lx = container.lx_;
        double ly = container.ly_;
        double max_radius = obstacles.max_radius_;

        double *dx_tot = new double [obstacles.n_obstacles_]();
        double *dy_tot = new double [obstacles.n_obstacles_]();

        int n_obstacles = obstacles.n_obstacles_;
        double d_avg = std::sqrt ( ( container.lx_*container.ly_ ) /n_obstacles );
        //double l_shift_max = 20.*std::sqrt ( ( container.lx_*container.ly_ ) - n_obstacles * obstacles.area_ ) / n_obstacles;
        double l_shift_max = 0.25 * std::sqrt ( (( container.lx_*container.ly_ ) - n_obstacles * obstacles.area_)  / n_obstacles);
        if ( GlobalFlags::verbose )
        {
            std::cout << "shift max " << l_shift_max << std::endl;
        }
        if ( d_avg > 0.1 * lx )
        {
            d_avg = 0.1 * ly;
        }
        if ( d_avg < obstacles.max_radius_ )
        {
            d_avg = obstacles.max_radius_;
        }

        CellList *celllist = new CellList ( n_obstacles, container.lx_, container.ly_, d_avg * 2 );
        celllist->fill ( obstacles );
        int n_acc = 0;
        for ( int i = 0; i < sweeps; i++ )
        {
            for ( int p = 0; p < n_obstacles; p++ )
            {
                int obs_nr = uni_dist ( engine ) * n_obstacles;

                double dx = ( 2*uni_dist ( engine )-1 ) * l_shift_max;
                double dy = ( 2*uni_dist ( engine )-1 ) * l_shift_max;

                double new_x = obstacles.x_[obs_nr] + dx;
                double new_y = obstacles.y_[obs_nr] + dy;

                //std::cout << "n: " << new_x << ' ' << new_y << std::endl;
                container.applyBoundary ( new_x, new_y );

                double old_x = obstacles.x_[obs_nr];
                double old_y = obstacles.y_[obs_nr];

                double shift_x = new_x - old_x;
                double shift_y = new_y - old_y;

                obstacles.shift ( obs_nr, shift_x, shift_y );

                int cell_x, cell_y;
                celllist->getCell ( new_x, new_y, cell_x, cell_y );


                bool accept = true;

                for ( int cx = -1; cx <= 1; cx++ )
                {
                    for ( int cy = -1; cy <= 1; cy++ )
                    {
                        int curr_cell_x = cell_x + cx;
                        int curr_cell_y = cell_y + cy;

                        celllist->remapCells ( curr_cell_x, curr_cell_y );

                        int cell_nr = curr_cell_x + curr_cell_y * celllist->nx_;
                        int current_obstacle = celllist->head_cells_[cell_nr];

                        while ( current_obstacle >= 0 )
                        {
                            int o_flag = 1;
                            int flag = 0;

                            if ( current_obstacle == obs_nr )
                            {
                                current_obstacle = celllist->list_cells_[current_obstacle];
                                continue;
                            }

                            int overlapp = obstacles.overlapp ( obs_nr, current_obstacle, 0.0, 0.0 );

                            if ( overlapp == -1 )
                            {
                                o_flag = -1;
                            }
                            if ( obstacles.x_[obs_nr] > lx - 2 *  max_radius )   // right border
                            {
                                flag += 1;
                            }
                            if ( obstacles.x_[obs_nr] < 2 * max_radius )     // left border
                            {
                                flag -= 1;
                            }
                            if ( obstacles.y_[obs_nr] > ly - 2 * max_radius )   // upper border
                            {
                                flag += 3;
                            }
                            if ( obstacles.y_[obs_nr] < 2 * max_radius )     // lowerborder
                            {
                                flag -= 3;
                            }

                            switch ( flag )
                            {
                                case 1:
                                {
                                    overlapp = obstacles.overlapp ( obs_nr, current_obstacle, -lx, 0 );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }
                                    break;
                                }
                                case -1:
                                {
                                    overlapp = obstacles.overlapp ( obs_nr, current_obstacle, +lx, 0 );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    break;
                                }

                                case 3:
                                {
                                    overlapp = obstacles.overlapp ( obs_nr, current_obstacle, 0, -ly );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    break;
                                }
                                case -3:
                                {
                                    int overlapp = obstacles.overlapp ( obs_nr, current_obstacle, 0,  +ly );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    break;
                                }

                                case 2:
                                {
                                    overlapp =obstacles.overlapp ( obs_nr, current_obstacle,  +lx,  -ly );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    overlapp =obstacles.overlapp ( obs_nr, current_obstacle,0, -ly );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    overlapp =obstacles.overlapp ( obs_nr, current_obstacle,  +lx,0 );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    break;
                                }

                                case 4:
                                {
                                    overlapp =obstacles.overlapp ( obs_nr, current_obstacle, -lx, -ly );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    overlapp =obstacles.overlapp ( obs_nr, current_obstacle, 0, -ly );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    overlapp =obstacles.overlapp ( obs_nr, current_obstacle,  -lx, 0 );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    break;
                                }

                                case -2:
                                {
                                    overlapp =obstacles.overlapp ( obs_nr, current_obstacle, -lx, +ly );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    overlapp =obstacles.overlapp ( obs_nr, current_obstacle, 0, +ly );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    overlapp =obstacles.overlapp ( obs_nr, current_obstacle, -lx, 0 );

                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    break;
                                }

                                case -4:
                                {
                                    overlapp = obstacles.overlapp ( obs_nr, current_obstacle, +lx, +ly );
                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    overlapp = obstacles.overlapp ( obs_nr, current_obstacle, 0, +ly );
                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }

                                    overlapp = obstacles.overlapp ( obs_nr, current_obstacle, +lx, 0 );
                                    if ( overlapp == -1 )
                                    {
                                        o_flag = -1;
                                    }
                                    break; 
                                }
                            }
                            if ( o_flag < 0 )
                            {
                                accept = false;
                                break;
                            }
                            else
                            {
                                accept = true;
                            }

                            current_obstacle = celllist->list_cells_[current_obstacle];
                        }
                        if ( accept == false )
                        {
                            break;
                        }

                    }
                    if ( accept == false )
                    {
                        break;
                    }

                }

                if ( accept == false )
                {
                    obstacles.shift ( obs_nr, -shift_x, -shift_y );
                }
                else
                {
                    dx_tot[obs_nr] += dx;
                    dy_tot[obs_nr] += dy;
                    n_acc++;
                    if ( dx_tot[obs_nr] * dx_tot[obs_nr] + dy_tot[obs_nr] * dy_tot[obs_nr] >= 0.25 * celllist->l_ )
                    {
                        celllist->fill ( obstacles );
                        for ( int bla = 0; bla < obstacles.n_obstacles_; bla++ )
                        {
                            dx_tot[bla] = 0;
                            dy_tot[bla] = 0;
                        }
                    }
                }
            }
        }
        delete celllist;
        return ( ( double ) n_acc/ ( sweeps * obstacles.n_obstacles_ ) );
    }
}
