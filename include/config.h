#pragma once
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string.h>

namespace Setup
{
struct Configuration
{
    //
    double magnetic_field_ = -1;

    std::string suffix_ = "";
    std::string configfile_ = "";
    std::string foldername_ = "lorentz";
    bool extra_folder_ = true;
    bool dump_config_ = false;
    bool binary_ = false;
    bool filter_circling_ = true;

    int n_traj_test_ = 0;
    int n_box_test_ = 0;

    //Box
    double lx_ = -1;
    double ly_ = -1;

    bool simple_nlist_ = false;

    //Obstacles
    std::string obstacle_type_ = "";
    int obstacle_type = -1;
    double r1_ = -1;
    double r2_ = -1;
    double l_edge_ = -1;
    double l_cross_ = -1;
    double d_cross_ = -1;
    double theta_rot_ =-1;
    double reduced_density_ = -1.0;
    double reduced_density_2_ = -1.0;
    bool overlapping_ = true;
    bool montcarlo_ = false;
    int mc_sweeps_ = 0;
    int random_seed_ = -1;
    int n_obstacles_2_ = 0;

    //trajectory
    double windtree_angle_ = 0.0;
    bool is_windtree_ = false;



    //Measurement
    int task_ = 0;
    int n_boxes_ = -1;
    int n_trajectory_ = -1;
    int n_obstacles_ = -1;
    int n_averages_ = -1;

    bool is_lin_scale_ = false;

    double t_begin_ = -1;
    double t_end_ = -1;
    int n_bins_log = -1;

    double vh_rmax_ = -1;
    int vh_nbin_space_ = -1;

    std::vector<double> q_;
    bool modes_ = false;

    Configuration ( int &argc, char ***argv )
    {
        while ( true )
        {
            int res = getopt ( argc, *argv, "B:c:f:D:d:r:s:L:F:t:T:N:SR:WZn:b:k:C" );
            if ( res == -1 )
            {
                break;
            }
            switch ( res ) //overrrides config file
            {
            case 'B': //magnetic field
            {
                magnetic_field_ = atof ( optarg );
                break;
            }
            case 'D': //Reduced density
            {
                reduced_density_ = atof ( optarg );
                break;
            }
            case 'd': //Reduces density scaterers
            {
                reduced_density_2_ = atof ( optarg );
                break;
            }
            case 'f': //folder
            {
                foldername_ = optarg;
                break;
            }

            case 'c': // config filename
            {
                configfile_ = optarg;
                break;
            }
            case 's': //filename suffix
            {
                suffix_ = optarg;
                break;
            }
            case 'r': // VanHove Rmax
            {
                vh_rmax_ = atof ( optarg );
                break;
            }
            case 't': // t begin
            {
                t_begin_ = atof ( optarg );
                break;
            }
            case 'T': // t end
            {
                t_end_ = atof ( optarg );
                break;
            }
            case 'L': // box dimension
            {
                lx_ = ly_ = atof ( optarg );
                break;
            }

            case 'N':
            {
                n_obstacles_ = atoi ( optarg );
                break;
            }
            case 'F': // no separate folder
            {
                extra_folder_ = false;
                foldername_ = optarg;
                break;
            }

            case 'S': // force simple neighborlist
            {
                simple_nlist_ = true;
                break;
            }
            case 'R': //radom set random seed
            {
                random_seed_ = atoi ( optarg );
                break;
            }

            case 'W': //write configs
            {
                dump_config_ = true;
                break;
            }

            case 'Z': //Binary system
            {
                binary_ = true;
                break;
            }
            case 'n': 
            {
                n_obstacles_2_ = atoi ( optarg );
                break;
            }
            
            case 'b': 
            {
                n_box_test_ = atoi ( optarg );
                break;
            }
            case 'k': 
            {
                n_traj_test_ = atoi ( optarg );
                break;
            }
            case 'C': 
            {
                filter_circling_ = false;
                break;
            }

            default:
            {
                std::cout << "unknown argument: " << res << std::endl;
                exit ( -1 );
            }
            }
        }

        std::ifstream configFile;
        configFile.open ( configfile_, std::ifstream::in );

        if ( configFile.is_open() )
        {

            std::string line, name, value;
            while ( getline ( configFile, line ) )
            {
                line.erase ( std::remove_if ( line.begin(), line.end(), isspace ), line.end() );
                if ( line[0] == '#' || line.empty() )
                {
                    continue;
                }

                auto delimiterPos = line.find ( "=" );

                if ( delimiterPos == std::string::npos )
                {
                    name = line;
                    value = "";
                }
                else
                {
                    name = line.substr ( 0, delimiterPos );
                    value = line.substr ( delimiterPos + 1 );
                }

                if ( name.compare ( "obstacle_type" ) == 0 )
                {
                    obstacle_type_ = value;
                }

                /*
                if ( lx_ < 0.0 or ly_ < 0.0 )
                {
                    if ( name.compare ( "lx" ) == 0 )
                    {
                        lx_ = stod ( value );
                    }
                    if ( name.compare ( "ly" ) == 0 )
                    {
                        ly_ = stod ( value );
                    }
                }

                if(n_obstacles_ < 0 and (lx_ < 0.0 or ly_ < 0.0 ))
                {
                    if ( name.compare ( "n_obstacles" ) == 0 )
                    {
                        n_obstacles_= stoi ( value );
                    }
                }
                */

                if ( name.compare ( "l_edge" ) == 0 )
                {
                    l_edge_ = stod ( value );
                }

                if ( name.compare ( "l_cross" ) == 0 )
                {
                    l_cross_ = stod ( value );
                }
                if ( name.compare ( "d_cross" ) == 0 )
                {
                    d_cross_ = stod ( value );
                }
                
                if ( name.compare ( "binary" ) == 0 )
                {
                    binary_ = true;
                }

                if ( name.compare ( "rotation" ) == 0 )
                {
                    theta_rot_= stod ( value );
                }
                if ( name.compare ( "windtree_angle" ) == 0 )
                {
                    windtree_angle_= stod ( value );
                }
                if ( name.compare ( "windtree" ) == 0 )
                {
                    is_windtree_ = true;
                }
                if ( name.compare ( "no_overlap" ) == 0 )
                {
                    overlapping_ = false;
                }

                if ( name.compare ( "montecarlo" ) == 0 )
                {
                    montcarlo_ = true;
                }

                if ( name.compare ( "mc_sweeps" ) == 0 )
                {
                    mc_sweeps_ = stod ( value );
                }

                if ( name.compare ( "r1" ) == 0 )
                {
                    r1_ = stod ( value );
                }
                if ( name.compare ( "r2" ) == 0 )
                {
                    r2_ = stod ( value );
                }

                if ( name.compare ( "n_boxes" ) == 0 )
                {
                    n_boxes_ = stoi ( value );
                }
                if ( name.compare ( "n_trajectory" ) == 0 )
                {
                    n_trajectory_ = stoi ( value );
                }

                if ( name.compare ( "n_averages" ) == 0 )
                {
                    n_averages_ = stoi ( value );
                }

                if ( name.compare ( "lin_scale" ) == 0 )
                {
                    is_lin_scale_ = true;
                }

                if ( name.compare ( "van_hove" ) == 0 )
                {
                    task_ = 1;
                }

                if ( name.compare ( "fsqt_lambda" ) == 0 )
                {
                    task_ = 2;
                    modes_ = false;

                    char * pch;
                    char * cp = strdup ( value.c_str() );
                    pch = strtok ( cp,",;" );
                    while ( pch != NULL )
                    {
                        double q = std::stod ( pch );
                        q_.push_back ( q );
                        pch = strtok ( NULL,",;" );
                    }
                    free ( cp );
                }

                if ( name.compare ( "fsqt_modes" ) == 0 )
                {
                    task_ = 2;
                    modes_ = true;

                    char * pch;
                    char * cp = strdup ( value.c_str() );
                    pch = strtok ( cp,",;" );
                    while ( pch != NULL )
                    {
                        double q = std::stod ( pch );
                        q_.push_back ( q );
                        pch = strtok ( NULL,",;" );
                    }
                    free ( cp );
                }

                if ( name.compare ( "moments" ) == 0 )
                {
                    task_ = 3;
                }

                if ( name.compare ( "space_bins" ) == 0 )
                {
                    vh_nbin_space_ = stoi ( value );
                }
                if ( name.compare ( "n_bins_log" ) == 0 )
                {
                    n_bins_log = stoi ( value );
                }
            }

            //check obstacle type
            if ( binary_ )
            {   
                if ( obstacle_type_.compare ( "circles" ) == 0 )
                {
                    obstacle_type = -1;
                }
                if ( obstacle_type_.compare ( "squares" ) == 0 )
                {
                    obstacle_type = -2;
                }
                if ( obstacle_type_.compare ( "crosses" ) == 0 )
                {
                    obstacle_type = -3;
                }
            }
            else
            {
                if ( obstacle_type_.compare ( "circles" ) == 0 )
                {
                    obstacle_type = 1;
                }
                if ( obstacle_type_.compare ( "squares" ) == 0 )
                {
                    obstacle_type = 2;
                }
                if ( obstacle_type_.compare ( "crosses" ) == 0 )
                {
                    obstacle_type = 3;
                }
            }
        }
        else
        {
            std::cerr << "Couldn't open config file "<< configfile_ <<  " for reading.\n";
        }
    }

    bool Check()
    {
        if ( task_ >= 1 )
        {
            if ( vh_nbin_space_ < 0 )
            {
                vh_nbin_space_ = 1000;
            }
            if ( vh_rmax_ < 0.0 )
            {
                vh_rmax_ = 3 * sqrt ( sqrt ( 2 ) * t_end_ / reduced_density_ );
            }
        }
        if ( reduced_density_ < 0 )
        {
            std::cout << "Please specify density in config" << std::endl;
            return false;
        }
        if ( !is_lin_scale_ )
        {
            if ( t_begin_ <= 0 )
            {
                std::cout << "t_begin cannot be 0 in Logscale" << std::endl;
                return false;
            }
        }
        return true;
    }
};
}
