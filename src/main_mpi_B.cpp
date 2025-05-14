
/***********Includes********/
#include <mpi.h>
#include <iostream>
#include <sys/stat.h>
#include <chrono>
#include <cmath>
#include <string.h>
#include <stdint.h>

/******* My Includes *******/
#include "simulationbox.h"
#include "neighbourlist.h"
#include "trajectory.h"
#include "timescales.h"
#include "config.h"
#include "propagators.h"
#include "dynamic_correlations_mpi.h"
#include "globals.h"
/***************************/

//CCW trajectories assumed

int main ( int argc, char **argv )
{
    ////////////////
    //On all Nodes//
    ////////////////

    //gl_tr_file.open("traj.txt");
    //gl_obs_file.open("obs.txt");

    // Setting up MPI stuff
    int MY_MPI_STATUS = MPI_Init ( &argc, &argv );
    int world_rank;
    int world_size;

    MPI_Comm_rank ( MPI_COMM_WORLD, &world_rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &world_size );

    //read in configuration
    Setup::Configuration config ( argc, &argv );
    if ( !config.Check() )
    {
        std::cout << "Config file: " << config.configfile_ << " error" << std::endl;
        return -1;
    }

    //Measurement
    const int n_boxes =         config.n_boxes_;
    const int n_trajectories =  config.n_trajectory_;

    if ( world_rank == 0 )
    {
        GlobalFlags::verbose = true;
    }
    else
    {
        GlobalFlags::verbose = false;
    }

    //setting the timescale:
    TimeScale::LogScale  *logscale = new TimeScale::LogScale ( config.t_begin_, config.t_end_, config.n_bins_log, config.n_averages_ );

    //result datasize
    int32_t R_DATA_SIZE = logscale->time_.size();
    int32_t VH_DATA_SIZE = R_DATA_SIZE * config.vh_nbin_space_;

    //Creating empty results
    CorrelationFunctions_Mpi::vacmsd_results *results;
    if ( config.task_ == 0 )
    {
        results = new CorrelationFunctions_Mpi::vacmsd_results ( R_DATA_SIZE );
    }
    if ( config.task_ == 1 )
    {
        results = new CorrelationFunctions_Mpi::vacmsd_results ( R_DATA_SIZE, *logscale, config.vh_rmax_, config.vh_nbin_space_ );
    }
    if ( config.task_ == 2 )
    {
        if ( ( config.lx_ != config.ly_ ) or ( config.lx_ < 0 or config.ly_ < 0 ) )
        {
            std::cout << "Warning lx =! ly OR NOT DEFINED problem for FSQT" << std::endl;
            std::abort();
        }
        results = new CorrelationFunctions_Mpi::vacmsd_results ( R_DATA_SIZE, *logscale, config.q_, config.lx_, config.modes_ );
    }
    if ( config.task_ == 3 )
    {
        results = new CorrelationFunctions_Mpi::vacmsd_results ( R_DATA_SIZE, 3 );
    }

    //----------------------------------------
    int traj_per_box = std::ceil ( n_trajectories / static_cast<double> ( n_boxes ) );
    int boxes_per_rank = std::floor ( static_cast<double> ( n_boxes ) / world_size );
    int rest_traj_per_box = 0;
    if ( boxes_per_rank == 0 )
    {
        boxes_per_rank = 1;
    }
    if ( traj_per_box * boxes_per_rank * world_size < n_trajectories )
    {
        rest_traj_per_box = std::ceil ( ( n_trajectories - traj_per_box * boxes_per_rank * world_size ) / static_cast<double> ( world_size ) );
    }

    std::chrono::high_resolution_clock::time_point t1, t2, t3;

    if ( world_rank == 0 )
    {
        GlobalFlags::verbose = true;

        std::cout << "MPI_STATUS: " << MY_MPI_STATUS << std::endl << std::endl;

        std::cout << std::endl << "------QUICK INFO------" << std::endl;

        if ( config.task_ == 0 )
        {
            std::cout << "Computing Msd" << std::endl;
        }
        if ( config.task_ == 1 )
        {
            std::cout << "Computing Van-Hove" << std::endl;
        }
        if ( config.task_ == 2 )
        {
            std::cout << "Computing Fsqt" << std::endl;
        }
        if ( config.task_ == 3 )
        {
            std::cout << "Computing Moments" << std::endl;
        }


        std::cout << "B_FIELD: " << config.magnetic_field_ << std::endl;
        std::cout << "RHO    : " << config.reduced_density_ << std::endl;
        std::cout << "RHO_ALT: " << config.n_obstacles_/ ( config.lx_*config.ly_ ) << std::endl;

        if ( config.binary_ )
        {
            std::cout << "RHO 2: " << config.reduced_density_2_ << std::endl;
        }

        std::cout << "Mpi Worldsize: " << world_size << std::endl;
        std::cout << "traj_per_box " << traj_per_box << std::endl;
        std::cout << "boxes_per_rank " << boxes_per_rank << std::endl;
        std::cout << "rest_traj_per_box " << rest_traj_per_box << std::endl;
        std::cout << "Data Lenght: " << R_DATA_SIZE << std::endl;

        if ( config.task_ >= 1 )
        {
            std::cout << "Van Hove Data Lengh: " << VH_DATA_SIZE << std::endl;
        }

        //Start Timing:
        t1 = std::chrono::high_resolution_clock::now();
    }

    //Box characteristics:
    double lx = config.lx_;
    double ly = config.ly_;
    int n_obstacles = 0;

    //Trajectory charactreistics
    double r_trajectory = 1.0/config.magnetic_field_;

    //setting the work for each node:
    int count = 0;
    for ( int i = 0 ; i <= boxes_per_rank; i++ )
    {
        int n = traj_per_box;
        if ( i == boxes_per_rank )
        {
            if ( rest_traj_per_box == 0 )
            {
                continue;
            }
            else
            {
                n = rest_traj_per_box;
            }
        }
        count += n;

        //create trajectory
        Trajectory::Circular trajectory ( n, config.t_end_, r_trajectory );

        //now switch between obstacles types
        switch ( config.obstacle_type )
        {
        case 1: //----------circles----------
        {
            //setting n_obstacles
            int n1;
            double max_radius = config.r1_ > config.r2_ ? config.r1_ : config.r2_;

            if ( config.n_obstacles_ < 0 )
            {
                n1 = static_cast<int> ( config.reduced_density_ * lx * ly / ( M_PI * config.r1_ * config.r1_ ) );
            }
            else
            {
                n1 = config.n_obstacles_;
                lx = std::sqrt (( M_PI * config.r1_ * config.r1_ ) * n1 / config.reduced_density_);
                ly = lx;
            }

            //create simulation box
            double target_cell_length = max_radius * 2.901;
            SimulationBox::PeriodicBoxCircularObstacles simbox ( lx, ly, target_cell_length, n1, config.r1_, config.overlapping_, true );

            //for big radii use circular cellists instead of simple ones
            if ( ( r_trajectory / target_cell_length ) > 10.0 and !config.simple_nlist_ )
            {
                //Create Stencil:
                Stencil::Circular stencil ( *simbox.celllist_, *simbox.container_, trajectory );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //estimate nlist lenght:
                int nlist_length = std::ceil ( simbox.obstacles_->n_tot_ * std::pow ( simbox.celllist_->l_,2 ) / ( simbox.container_->lx_ * simbox.container_->ly_ ) * 2 ) * stencil.max_lenght_;

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }
                Neighbourlist::Circular neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                CircularPropagator propagator ( simbox, trajectory, neighlist, config.filter_circling_);
                propagator.setStartPositions ( simbox.obstacles_ );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            else
            {
                //Create Stencil:
                //simple stencil cutoff:
                double cutoff = trajectory.radius_ + simbox.obstacles_->max_radius_;
                Stencil::FullPeriodicTrimmedStencil stencil ( *simbox.celllist_, cutoff );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //now we need to approximate the amount of neighbour obstacles:
                //we take the first approximation of the number of obstacles in the stencil
                int nlist_length = std::ceil ( ( simbox.obstacles_->n_tot_ / static_cast<double> ( simbox.celllist_->n_cells_ ) + 9 ) * 1.5 * stencil.n_stencil_ );

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }

                Neighbourlist::Simple neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                SimplePropagator propagator ( simbox, trajectory, neighlist, config.filter_circling_ );

                //set start positions:
                propagator.setStartPositions ( simbox );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            break;
        }

        case 2: //----------squares----------
        {
            //setting n_obstacles
            if ( config.n_obstacles_ < 0 )
            {
                n_obstacles = static_cast<int> ( config.reduced_density_ * lx * ly / ( config.l_edge_ * config.l_edge_ ) );
            }
            else
            {
                n_obstacles = config.n_obstacles_;
                lx = ly = config.l_edge_ * std::sqrt ( n_obstacles/config.reduced_density_ );
            }

            //create simulation box
            double max_radius = config.l_edge_ / std::sqrt ( 2.0 );
            double target_cell_length = max_radius * 2.901;

            SimulationBox::PeriodicBoxSquareObstacles simbox ( lx, ly, target_cell_length, n_obstacles, config.l_edge_, config.theta_rot_, config.overlapping_, config.montcarlo_, config.mc_sweeps_, true );


            //for big radii use circular cellists instead of simple ones
            if ( ( r_trajectory / target_cell_length ) > 10.0 and !config.simple_nlist_ )
            {
                //Create Stencil:
                Stencil::Circular stencil ( *simbox.celllist_, *simbox.container_, trajectory );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //estimate nlist lenght:
                int nlist_length = std::ceil ( simbox.obstacles_->n_tot_ * std::pow ( simbox.celllist_->l_,2 ) / ( simbox.container_->lx_ * simbox.container_->ly_ ) * 2 ) * stencil.max_lenght_;

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }
                Neighbourlist::Circular neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                CircularPropagator propagator ( simbox, trajectory, neighlist, config.filter_circling_ );
                propagator.setStartPositions ( simbox.obstacles_ );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            else
            {
                //Create Stencil:
                //simple stencil cutoff:
                double cutoff = trajectory.radius_ + simbox.obstacles_->max_radius_;
                Stencil::FullPeriodicTrimmedStencil stencil ( *simbox.celllist_, cutoff );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //now we need to approximate the amount of neighbour obstacles:
                //we take the first approximation of the number of obstacles in the square sourrounding the stencil
                int nlist_length = std::ceil ( ( simbox.obstacles_->n_tot_ / static_cast<double> ( simbox.celllist_->n_cells_ ) + 9 ) * 1.5 * stencil.n_stencil_ );

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }

                Neighbourlist::Simple neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                SimplePropagator propagator ( simbox, trajectory, neighlist, config.filter_circling_ );

                //set start positions:
                propagator.setStartPositions ( simbox );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            break;//--------------------------
        }
        case 3: //----------crosses----------
        {
            //setting n_obstacles
            double cross_area = 2 * config.l_cross_ * config.d_cross_ - config.d_cross_ * config.d_cross_;

            if ( config.n_obstacles_ < 0 )
            {
                n_obstacles = static_cast<int> ( config.reduced_density_ * lx * ly / cross_area );
            }
            else
            {
                n_obstacles = config.n_obstacles_;
                lx = ly = std::sqrt ( n_obstacles * cross_area / config.reduced_density_ );
            }

            //create simulation box
            double max_radius = std::sqrt ( 0.25 * ( config.l_cross_ * config.l_cross_ + config.d_cross_ * config.d_cross_ ) );

            double target_cell_length = max_radius * 2.901;
            SimulationBox::PeriodicBoxCrossObstacles simbox ( lx, ly, target_cell_length, n_obstacles, config.l_cross_, config.d_cross_, config.theta_rot_, config.overlapping_, true );

            //for big radii use circular cellists instead of simple ones
            if ( ( r_trajectory / target_cell_length ) > 10.0 and !config.simple_nlist_ )
            {
                //Create Stencil:
                Stencil::Circular stencil ( *simbox.celllist_, *simbox.container_, trajectory );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //estimate nlist lenght:
                int nlist_length = std::ceil ( simbox.obstacles_->n_tot_ * std::pow ( simbox.celllist_->l_,2 ) / ( simbox.container_->lx_ * simbox.container_->ly_ ) * 2 ) * stencil.max_lenght_;

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }
                Neighbourlist::Circular neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                CircularPropagator propagator ( simbox, trajectory, neighlist, config.filter_circling_ );
                propagator.setStartPositions ( simbox.obstacles_ );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            else
            {
                //Create Stencil:
                //simple stencil cutoff:
                double cutoff = trajectory.radius_ + simbox.obstacles_->max_radius_;
                Stencil::FullPeriodicTrimmedStencil stencil ( *simbox.celllist_, cutoff );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //now we need to approximate the amount of neighbour obstacles:
                //we take the first approximation of the number of obstacles in the square sourrounding the stencil
                int nlist_length = std::ceil ( ( simbox.obstacles_->n_tot_ / static_cast<double> ( simbox.celllist_->n_cells_ ) + 9 ) * 1.5 * stencil.n_stencil_ );

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }

                Neighbourlist::Simple neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                SimplePropagator propagator ( simbox, trajectory, neighlist, config.filter_circling_ );

                //set start positions:
                propagator.setStartPositions ( simbox );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            break;
        }
        case -1: //----------Binary----------
        {
            int nobst;
            //setting n_obstacles1 (circles)
            if ( config.n_obstacles_ < 0 )
            {
                nobst = static_cast<int> ( config.reduced_density_ * lx * ly / ( M_PI * config.r1_ * config.r1_ ) );
            }
            else
            {
                nobst = config.n_obstacles_;
                lx = std::sqrt (( M_PI * config.r1_ * config.r1_ ) * nobst / config.reduced_density_);
                ly = lx;
            }
            
            SimulationBox::PeriodicBoxCircularObstacles simbox1 ( lx, ly, 1.0, nobst, config.r1_, config.overlapping_, false );
            
            //setting n_obstacles
            int n1;
            if ( config.n_obstacles_2_ <= 0 )
            {
                n1 = static_cast<int> ( config.reduced_density_2_ * lx * ly / ( M_PI * config.r2_ * config.r2_ ) );
            }
            else
            {
                std::cout << "You need to specify the density not the particle number of the secon obstacle type" << std::endl;
                return -1;
            }

            SimulationBox::PeriodicBoxCircularObstacles simbox2 ( lx, ly, 1, n1, config.r2_, config.overlapping_, false );
            
            //create foldername to dump configs-----------------------------------------
            char strbuff[15];
            std::snprintf ( strbuff,15,"D_%.4G", config.reduced_density_ );
            std::string density_str ( strbuff );

            std::snprintf ( strbuff,15,"T_%.4G", config.t_end_ );
            std::string tend_str ( strbuff );

            // write to file
            std::string foldername = config.foldername_;
            mkdir ( foldername.c_str(), 0755 );

            if ( config.extra_folder_ )
            {
                foldername += '/' + density_str + '_' + std::to_string ( static_cast<int> ( lx ) ) + 'x' + std::to_string ( static_cast<int> ( ly ) );
                mkdir ( foldername.c_str(), 0755 );
            }
            std::string filename = std::to_string ( world_rank ) + '_' + std::to_string ( i ) + ".txt";

            //---------------------------------------------------------------------------

            double max_radius = ( config.r2_ > config.r1_ ) ? config.r2_ : config.r1_;
            double target_cell_length = max_radius * 2.901;

            SimulationBox::BinaryBox < SimulationBox::PeriodicBoxCircularObstacles, SimulationBox::PeriodicBoxCircularObstacles, Obstacles::Circles, Obstacles::Circles> simbox ( &simbox1, &simbox2, simbox1.obstacles_, simbox2.obstacles_, target_cell_length );

            if ( config.dump_config_ == true )
            {
                simbox.writeObstacles_config ( foldername, filename );
            }

            //for big radii use circular cellists instead of simple ones
            if ( ( r_trajectory / target_cell_length ) > 10.0 and !config.simple_nlist_ )
            {
                //Create Stencil:
                Stencil::Circular stencil ( *simbox.celllist_, *simbox.container_, trajectory );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //estimate nlist lenght:
                int nlist_length = std::ceil ( simbox.obstacles_->n_tot_ * std::pow ( simbox.celllist_->l_,2 ) / ( simbox.container_->lx_ * simbox.container_->ly_ ) * 2 ) * stencil.max_lenght_;

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }

                Neighbourlist::Circular neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                CircularPropagatorBinary propagator ( simbox, trajectory, neighlist, config.filter_circling_ );
                propagator.setStartPositions ( simbox.obstacles_ );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            else
            {

                //Create Stencil:
                //simple stencil cutoff:
                double cutoff = trajectory.radius_ + simbox.obstacles_->max_radius_;
                Stencil::FullPeriodicTrimmedStencil stencil ( *simbox.celllist_, cutoff );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //now we need to approximate the amount of neighbour obstacles:
                //we take the first approximation of the number of obstacles in the square sourrounding the stencil
                int nlist_length = std::ceil ( ( simbox.obstacles_->n_tot_ / static_cast<double> ( simbox.celllist_->n_cells_ ) + 9 ) * 1.5 * stencil.n_stencil_ );

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }

                Neighbourlist::Simple neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                SimplePropagatorBinary propagator ( simbox, trajectory, neighlist, config.filter_circling_ );

                //set start positions:
                propagator.setStartPositions ( simbox.obstacles_ );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            break;
        }
        case -2: //----------Binary----------
        {
            double max_radius_squares = config.l_edge_ / sqrt ( 2 );
            //setting n_obstacles1 (Squares)
            if ( config.n_obstacles_ < 0 )
            {
                n_obstacles = static_cast<int> ( config.reduced_density_ * lx * ly / ( config.l_edge_ * config.l_edge_ ) );
            }
            else
            {
                n_obstacles = config.n_obstacles_;
                lx = ly = config.l_edge_ * std::sqrt ( n_obstacles/config.reduced_density_ );
            }
            SimulationBox::PeriodicBoxSquareObstacles simbox1 ( lx, ly, 1, n_obstacles, config.l_edge_, config.theta_rot_, config.overlapping_, config.montcarlo_, config.mc_sweeps_, false );


            //setting n_obstacles
            int n1;
            double max_radius_circles = config.r1_;

            if ( config.n_obstacles_2_ <= 0 )
            {
                n1 = static_cast<int> ( config.reduced_density_2_ * lx * ly / ( M_PI * config.r1_ * config.r1_ ) );
            }
            else
            {
                std::cout << "You need to specify the density not the particle number of the second obstacle type" << std::endl;
                return -1;
            }
            
            SimulationBox::PeriodicBoxCircularObstacles simbox2 ( lx, ly, 1, n1, config.r1_, config.overlapping_, true );

            //create foldername to dump configs-----------------------------------------
            char strbuff[15];
            std::snprintf ( strbuff,15,"D_%.4G", config.reduced_density_ );
            std::string density_str ( strbuff );

            std::snprintf ( strbuff,15,"T_%.4G", config.t_end_ );
            std::string tend_str ( strbuff );

            // write to file
            std::string foldername = config.foldername_;
            mkdir ( foldername.c_str(), 0755 );

            if ( config.extra_folder_ )
            {
                foldername += '/' + density_str + '_' + std::to_string ( static_cast<int> ( lx ) ) + 'x' + std::to_string ( static_cast<int> ( ly ) );
                mkdir ( foldername.c_str(), 0755 );
            }
            std::string filename = std::to_string ( world_rank ) + '_' + std::to_string ( i ) + ".txt";

            //---------------------------------------------------------------------------


            double max_radius = ( max_radius_circles > max_radius_squares ) ? max_radius_circles : max_radius_squares;
            double target_cell_length = max_radius * 2.901;

            SimulationBox::BinaryBox <SimulationBox::PeriodicBoxSquareObstacles, SimulationBox::PeriodicBoxCircularObstacles, Obstacles::Squares, Obstacles::Circles> simbox ( &simbox1, &simbox2, simbox1.obstacles_, simbox2.obstacles_, target_cell_length );

            if ( config.dump_config_ == true )
            {
                simbox.writeObstacles_config ( foldername, filename );
            }

            //for big radii use circular cellists instead of simple ones
            if ( ( r_trajectory / target_cell_length ) > 10.0 and !config.simple_nlist_ )
            {
                //Create Stencil:
                Stencil::Circular stencil ( *simbox.celllist_, *simbox.container_, trajectory );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //estimate nlist lenght:
                int nlist_length = std::ceil ( simbox.obstacles_->n_tot_ * std::pow ( simbox.celllist_->l_,2 ) / ( simbox.container_->lx_ * simbox.container_->ly_ ) * 2 ) * stencil.max_lenght_;

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }

                Neighbourlist::Circular neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                CircularPropagatorBinary propagator ( simbox, trajectory, neighlist, config.filter_circling_ );
                propagator.setStartPositions ( simbox.obstacles_ );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            else
            {

                //Create Stencil:
                //simple stencil cutoff:
                double cutoff = trajectory.radius_ + simbox.obstacles_->max_radius_;
                Stencil::FullPeriodicTrimmedStencil stencil ( *simbox.celllist_, cutoff );

                //Create placement stencil:
                Stencil::FullPeriodicTrimmedStencil placement_stencil ( *simbox.celllist_, simbox.obstacles_->max_radius_ );

                //initialize accelerators:
                //now we need to approximate the amount of neighbour obstacles:
                //we take the first approximation of the number of obstacles in the square sourrounding the stencil
                int nlist_length = std::ceil ( ( simbox.obstacles_->n_tot_ / static_cast<double> ( simbox.celllist_->n_cells_ ) + 9 ) * 1.5 * stencil.n_stencil_ );

                // cut to max value
                if ( nlist_length > simbox.obstacles_->n_tot_ )
                {
                    nlist_length = simbox.obstacles_->n_tot_;
                }

                Neighbourlist::Simple neighlist ( nlist_length, *simbox.celllist_, stencil, placement_stencil );

                //initialize propagator:
                SimplePropagatorBinary propagator ( simbox, trajectory, neighlist, config.filter_circling_);

                //set start positions:
                propagator.setStartPositions ( simbox.obstacles_ );

                //calculate
                CorrelationFunctions_Mpi::VacMsd_Mpi ( simbox, trajectory, propagator,  *results, *logscale );
            }
            break;
        }
        }

        if ( world_rank == 0 and i == 0 )
        {
            t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = ( t2 - t1 ) * ( boxes_per_rank + 1 ) /  60.;
            std::cout << std::endl << std::endl;
            std::cout << "projected duration: " << duration.count() << " minutes" << std::endl;
            GlobalFlags::verbose = false;
        }
    }

//------------------------------------------//
//-----------Results Aggregation------------//
//------------------------------------------//

    if ( world_rank == 0 )
    {
        MPI_Reduce ( MPI_IN_PLACE, results->msd_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( MPI_IN_PLACE, results->vaf_xx_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( MPI_IN_PLACE, results->vaf_yy_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( MPI_IN_PLACE, results->vaf_xy_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( MPI_IN_PLACE, results->vaf_yx_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( MPI_IN_PLACE, results->m4_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

        if ( results->task_ == 1 )
        {
            MPI_Reduce ( MPI_IN_PLACE, results->van_hove_self_count_, VH_DATA_SIZE, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD );

            MPI_Reduce ( MPI_IN_PLACE, results->m6_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce ( MPI_IN_PLACE, results->m8_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce ( MPI_IN_PLACE, results->m10_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        }
        if ( results->task_ == 2 )
        {
            MPI_Reduce ( MPI_IN_PLACE, results->fsqt_, R_DATA_SIZE * ( int ) results->q_.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        }

        if ( results->task_ == 3 )
        {
            MPI_Reduce ( MPI_IN_PLACE, results->m6_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce ( MPI_IN_PLACE, results->m8_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce ( MPI_IN_PLACE, results->m10_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        }

        MPI_Reduce ( MPI_IN_PLACE, results->normalisator_, R_DATA_SIZE, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD );

        t3 = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration = ( t3 - t1 ) / 60.;
        std::cout << "duration: " << duration.count() << " minutes with mpi world_size: " << world_size << std::endl;

        results->normalize();

        std::string otss;
        switch ( config.obstacle_type )
        {
        case 1 :
            otss = "crcl";
            break;
        case 2:
            otss = "sqr";
            break;
        case 3:
            otss = "crss";
            break;
        case -1 :
            otss = "crcl_bin";
            break;
        case -2:
            otss = "sqr_bin";
            break;
        case -3:
            otss = "crss_bin";
            break;
        }

        char strbuff[15];

        std::snprintf ( strbuff,15,"B_%.4G", config.magnetic_field_ );
        std::string bfield_str ( strbuff );

        std::snprintf ( strbuff,15,"D_%.4G", config.reduced_density_ );
        std::string density_str ( strbuff );

        std::snprintf ( strbuff,15,"T_%.4G", config.t_end_ );
        std::string tend_str ( strbuff );

        // write to file
        std::string foldername = config.foldername_;
        mkdir ( foldername.c_str(), 0755 );

        if ( config.extra_folder_ )
        {
            foldername += "/" + bfield_str + "_" + density_str + '_' + std::to_string ( static_cast<int> ( lx ) ) + 'x' + std::to_string ( static_cast<int> ( ly ) );
            mkdir ( foldername.c_str(), 0755 );
        }

        std::string filename = foldername + "/" + bfield_str + "_" + density_str + '_' + tend_str + "_N_" + std::to_string ( n_trajectories ) + '_' + config.suffix_  + otss;

        std::ofstream result_file;
        IO::makeFilenameUnique ( filename, ".txt" );
        result_file.open ( filename + ".txt" );

        //Print Header:

        result_file << config.magnetic_field_ << " B Field;" << " obstacles: " << config.obstacle_type_ << ' '
                    << "; overlapping: " << config.overlapping_ << "; rotation: " << config.theta_rot_
                    << "; is windtree: " << config.is_windtree_<< std::endl;
        result_file << config.reduced_density_ << " Reduced density" << std::endl;
        result_file << lx << " Container lx" << std::endl;
        result_file << ly << " Container ly" << std::endl;
        result_file << config.t_end_ << " Lifetime" << std::endl;
        result_file << config.n_averages_ << " Time Averages" << std::endl;
        result_file << n_trajectories << " Number of Trajectories" << std::endl;
        result_file << n_obstacles << " Number of Obstacles" << std::endl;
        result_file << "duration: " << duration.count() << " minutes with mpi world_size: " << world_size << std::endl;

        result_file.precision ( 17 );

        for ( auto i = 0; i < R_DATA_SIZE ; i++ )
        {
            result_file << logscale->time_[i] << ' ' << results->msd_[i] << ' ' <<  results->m4_[i] << ' ' << results->vaf_xx_[i] << ' ' << results->vaf_yy_[i] << ' ' << results->vaf_xy_[i] << ' ' << results->vaf_yx_[i] << std::endl;
        }
        result_file.close();

        // ********************van hove part*********************
        if ( results->task_ == 1 )
        {
            std::ofstream result_file_vh;
            result_file_vh.open ( filename + ".vhs" );

            //Print Header:
            result_file_vh << config.magnetic_field_ << " B Field;" << " obstacles: " << config.obstacle_type_ << ' '
                           << "; overlapping: " << config.overlapping_ << "; rotation: " << config.theta_rot_
                           << "; is windtree: " << config.is_windtree_<< std::endl;
            result_file_vh << config.reduced_density_ << " Reduced density" << std::endl;
            result_file_vh << lx << " Container lx" << std::endl;
            result_file_vh << ly << " Container ly" << std::endl;
            result_file_vh << config.t_end_ << " Lifetime" << std::endl;
            result_file_vh << config.n_averages_ << " Time Averages" << std::endl;
            result_file_vh << n_obstacles << " Number of Obstacles" << std::endl;
            result_file_vh << n_trajectories << " Number of Trajectories" << std::endl;
            result_file_vh << "duration: " << duration.count() << " minutes with mpi world_size: " << world_size << std::endl;

            result_file_vh.precision ( 17 );
            //write t
            for ( auto i = 0; i < R_DATA_SIZE ; i++ )
            {
                result_file_vh << logscale->time_[i] << ' ';
            }
            result_file_vh << std::endl;

            result_file_vh << results->rmin_ << ' ' << results->rmax_ << ' ' << results->nbingr_ <<std::endl;


            for ( auto i = 0; i < R_DATA_SIZE ; i++ )
            {
                for ( auto j = 0; j < config.vh_nbin_space_; j++ )
                {
                    result_file_vh << results->van_hove_self_[i * config.vh_nbin_space_ + j] << ' ';
                }
                result_file_vh << std::endl;
            }
            result_file_vh.close();

            //----moments file---//
            std::ofstream result_file_mom;
            result_file_mom.open ( filename + ".mom" );

            //Print Header:
            result_file_mom << config.magnetic_field_ << " B Field" << std::endl;
            result_file_mom << config.reduced_density_ << " Reduced density" << std::endl;
            result_file_mom << lx << " Container lx" << std::endl;
            result_file_mom << ly << " Container ly" << std::endl;
            result_file_mom << config.t_end_ << " Lifetime" << std::endl;
            result_file_mom << config.n_averages_ << " Time Averages" << std::endl;
            result_file_mom << n_trajectories << " Number of Trajectories" << std::endl;
            result_file_mom << n_obstacles << " Number of Obstacles" << std::endl;
            result_file_mom << "duration: " << duration.count() << " minutes with mpi world_size: " << world_size << std::endl;

            result_file_mom.precision ( 17 );

            for ( auto i = 0; i < R_DATA_SIZE ; i++ )
            {
                result_file_mom << logscale->time_[i] << ' ' << results->msd_[i] << ' ' <<  results->m4_[i] << ' ' << results->m6_[i] << ' ' << results->m8_[i] << ' ' << results->m10_[i] << std::endl;
            }
            result_file_mom.close();

        }

        // **********************fsqt_ part***********************
        if ( results->task_ == 2 )
        {
            std::ofstream result_file_fs;
            result_file_fs.open ( filename+ '_' + std::to_string ( results->qmin_ ).substr ( 0,5 ) + "to" + std::to_string ( results->qmax_ ).substr ( 0,5 ) +".fsqt" );

            //Print Header:
            result_file_fs << config.magnetic_field_ << " B Field;" << " obstacles: " << config.obstacle_type_ << ' '
                           << "; overlapping: " << config.overlapping_ << "; rotation: " << config.theta_rot_
                           << "; is windtree: " << config.is_windtree_<< std::endl;
            result_file_fs << config.reduced_density_ << " Reduced density" << std::endl;
            result_file_fs << lx << " Container lx" << std::endl;
            result_file_fs << ly << " Container ly" << std::endl;
            result_file_fs << config.t_end_ << " Lifetime" << std::endl;
            result_file_fs << config.n_averages_ << " Time Averages" << std::endl;
            result_file_fs << n_obstacles << " Number of Obstacles" << std::endl;
            result_file_fs << n_trajectories << " Number of Trajectories" << std::endl;
            result_file_fs << "duration: " << duration.count() << " minutes with mpi world_size: " << world_size << std::endl;

            result_file_fs.precision ( 10 );

            for ( auto j = 0; j < static_cast<int> ( results->q_.size() ); ++j )
            {
                result_file_fs << results->q_[j] << ' ' ;
            }
            result_file_fs << std::endl;

            for ( auto i = 0; i < R_DATA_SIZE ; i++ )
            {
                result_file_fs << logscale->time_[i] << ' ';
                for ( auto j = 0; j < static_cast<int> ( results->q_.size() ); ++j )
                {
                    result_file_fs << results->fsqt_[R_DATA_SIZE * j + i] << ' ';
                }
                result_file_fs  << std::endl;
            }
            result_file_fs.close();
        }
        // *********************moments_ part **********************//
        if ( results->task_ == 3 )
        {
            std::ofstream result_file_mom;
            result_file_mom.open ( filename + ".mom" );

            //Print Header:
            result_file_mom << config.magnetic_field_ << " B Field" << std::endl;
            result_file_mom << config.reduced_density_ << " Reduced density" << std::endl;
            result_file_mom << lx << " Container lx" << std::endl;
            result_file_mom << ly << " Container ly" << std::endl;
            result_file_mom << config.t_end_ << " Lifetime" << std::endl;
            result_file_mom << config.n_averages_ << " Time Averages" << std::endl;
            result_file_mom << n_trajectories << " Number of Trajectories" << std::endl;
            result_file_mom << n_obstacles << " Number of Obstacles" << std::endl;
            result_file_mom << "duration: " << duration.count() << " minutes with mpi world_size: " << world_size << std::endl;

            result_file_mom.precision ( 10 );

            for ( auto i = 0; i < R_DATA_SIZE ; i++ )
            {
                result_file_mom << logscale->time_[i] << ' ' << results->msd_[i] << ' ' <<  results->m4_[i] << ' ' << results->m6_[i] << ' ' << results->m8_[i] << ' ' << results->m10_[i] << std::endl;
            }
            result_file_mom.close();
        }

    }
    else
    {
        MPI_Reduce ( results->msd_, results->msd_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( results->vaf_xx_, results->vaf_xx_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( results->vaf_yy_, results->vaf_yy_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( results->vaf_xy_, results->vaf_xy_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( results->vaf_yx_, results->vaf_yx_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce ( results->m4_, results->m4_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

        if ( results->task_ == 1 )
        {
            MPI_Reduce ( results->van_hove_self_count_, results->van_hove_self_count_, VH_DATA_SIZE, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD );

            MPI_Reduce ( results->m6_, results->m6_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce ( results->m8_, results->m8_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce ( results->m10_, results->m10_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        }

        if ( results->task_ == 2 )
        {
            MPI_Reduce ( results->fsqt_, results->fsqt_, R_DATA_SIZE * ( int ) results->q_.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        }

        if ( results->task_ == 3 )
        {
            MPI_Reduce ( results->m6_, results->m6_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce ( results->m8_, results->m8_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce ( results->m10_, results->m10_, R_DATA_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        }

        MPI_Reduce ( results->normalisator_, results->normalisator_, R_DATA_SIZE, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD );
    }

    int rest = 0;
    if ( rest_traj_per_box != 0 )
    {
        rest = 1;
    }

    if ( world_rank == 0 )
    {
        std::cout << "averaged over " << count * world_size << " trajectories in " << ( boxes_per_rank + rest ) * world_size << " boxes" << std::endl;
    }

    delete logscale;
    delete results;

    MPI_Finalize();

    return 0;
}

