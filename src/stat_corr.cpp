#include <iostream>
#include <iomanip>
#include <experimental/filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <unistd.h>

#include "stat_corr.h"

int main ( int argc, char **argv )
{
    std::vector<std::string> filelist;
    std::string path_str;
    int n_max_files = 500;
    int nbin_gr = 400;
    int nbin_sq = 200;
    double lengthscale_gr = 25;
    double lengthscale_sq =25;
    int n_shell_avg_sq = 150;

    bool do_sq = true;

    while ( 1 )
    {
        int result = getopt ( argc, argv, "f:N:l:q:n:k:gha:" );
        if ( result == -1 ) break;
        switch ( result )
        {
            case '?':
                std::cout << "unkown parameter " << optarg << std::endl;
                break;
            case ':': // missing argument of a parameter
                fprintf ( stderr, "missing argument.\n" );
                break;
            case 'f':
                path_str = optarg;
                break;
            case 'N':
                n_max_files = atoi ( optarg );
                break;
            case 'n':
                nbin_gr = atoi ( optarg );
                break;
            case 'k':
                nbin_sq = atoi ( optarg );
                break;
            case 'l':
                lengthscale_gr = atof ( optarg );
                break;
            case 'q':
                lengthscale_sq = atof ( optarg );
                break;
                
            case 'a':
                n_shell_avg_sq = atof ( optarg );
                break;
                
            case 'g':
                do_sq = false;
                break;
            case 'h':
                std::cout << "options:" << std::endl;
                std::cout << "-f:\t <filepath>" << std::endl;
                std::cout << "-N:\t maximum numbers of files" << std::endl;
                std::cout << "-l:\t lenghtscale pair correlation function" << std::endl;
                std::cout << "-q:\t lenghtscale in q space structure factor" << std::endl;
                std::cout << "-g:\t compute gr only" << std::endl;
                std::cout << "-n:\t nbin gr" << std::endl;
                std::cout << "-k:\t nbin sq" << std::endl;
                std::cout << "-h:\t help" << std::endl;
                break;
        }
    }
    //----------REMOVE IF NOT TESTING MANUALLY ------------------------------
    path_str="/home/stylx/temp_results/stat_corr_thesis/no_overlapp/sq_0.5_mc_0/";
    //-----------------------------------------------------------------------
    
    std::cout << std::endl << "path: " << path_str << std::endl;

    if ( path_str.empty() )
    {
        std::cout << "please provide path" << std::endl;
        return -1;
    }

    for ( const auto &file : std::experimental::filesystem::directory_iterator ( path_str ) )
    {
        if ( file.path().extension() == ".txt" )
        {
            filelist.push_back ( file.path() );
        }
    }

    if ( filelist.size() == 0 )
    {
        std::cout << "no valid file in folder" << std::endl;
    }

    int n_systems = 0;

    std::vector<double> gr ( nbin_gr );
    std::vector<double> gr_avg ( nbin_gr, 0.0 );
    std::vector<double> bins_gr;

    std::vector<double> sq;
    std::vector<double> sq_avg;
    std::vector<double> bins_sq;

    if ( do_sq )
    {
        sq.resize ( nbin_sq );
        sq_avg.resize ( nbin_sq, 0.0 );
    }

    if ( int(filelist.size()) < n_max_files )
    {
        n_max_files = filelist.size();
    }

    for ( int i = 0; i < n_max_files ; i++ )
    {
        std::cout << "processing file " << filelist[i] << '('  << n_systems+1 << "/" << n_max_files << ')' << std::endl;
        std::ifstream file;
        file.open ( filelist[i] );

        if ( !file.is_open() )
        {
            std::cout << "Could not open file: " << filelist[i] << std::endl;
            return -1;
        }

        std::string line;
        std::stringstream ss;

        //------Open the file and begin analysis

        int n_particles;
        double lx, ly;
        int sides;

        std::getline ( file, line );

        ss << line;
        ss >> n_particles;

        ss.clear();
        ss.str ( std::string() );
        std::getline ( file, line );
        ss << line;
        ss >> lx >> ly;

        ss.clear();
        ss.str ( std::string() );
        std::getline ( file, line );
        ss << line;
        ss >> sides;

        std::vector<double>x ( n_particles );
        std::vector<double>y ( n_particles );

        for ( int n = 0; n < n_particles; n++ )
        {
            ss.clear();
            ss.str ( std::string() );
            std::getline ( file, line );
            ss << line;
            ss >> x[n] >> y[n];
        }

        CorrelationFunctions::pairCorrelation ( x.data(), y.data(), lx, ly, lengthscale_gr, n_particles, gr, bins_gr );
        for ( size_t n = 0; n < gr.size(); ++n )
        {
            gr_avg[n] += gr[n];
        }


        if ( do_sq )
        {
            CorrelationFunctions::structureFactor ( x.data(), y.data(), lx, ly, lengthscale_sq, n_shell_avg_sq, n_particles, sq, bins_sq );
            for ( size_t n = 0; n < sq.size(); ++n )
            {
                sq_avg[n] += sq[n];
            }
        }

        n_systems++;
    }

    std::ofstream outfile_gr;
    outfile_gr.open ( path_str + "/gr.dat" );
    outfile_gr.precision ( 10 );
    
    for ( size_t n = 0; n < gr.size(); ++n )
    {
        outfile_gr << bins_gr[n] << ' ' << gr_avg[n] /n_systems << std::endl;
    }

    if ( do_sq )
    {
        std::ofstream outfile_sq;
        outfile_sq.open ( path_str + "/sq.dat" );
        outfile_sq.precision ( 10 );

        for ( size_t n = 0; n < sq.size(); ++n )
        {
            outfile_sq << bins_sq[n] << ' ' << sq_avg[n] /n_systems << std::endl;
        }
    }
    std::cout << "done" << std::endl;
}
