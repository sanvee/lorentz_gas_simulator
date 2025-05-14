#pragma once

#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <vector>

#include "container.h"
#include "globals.h"

namespace CorrelationFunctions
{

void minimalImage ( double &dx, double &dy, double lx, double ly )
{
    if ( dx <= -0.5 * lx )
    {
        dx += lx;
    }
    if ( dx >   0.5 * lx )
    {
        dx -= lx;
    }

    if ( dy <= -0.5 * ly )
    {
        dy += ly;
    }
    if ( dy >   0.5 * ly )
    {
        dy -= ly;
    }
}

void pairCorrelation ( double* x, double* y, double lx, double ly, double lengthscale, int npart, std::vector<double> &gr, std::vector<double> &bins )
{

    //std::cout << "computing pair correlation function ..." << std::endl;

    //set the number of bins!!
    int nbins = gr.size();
    
    //choose the smallest dimension on the container
    
    double length;
    lx > ly ? length = lx : length = ly; 
    //std::cout << std::setprecision(10) << ly << ' ' << ly << std::endl;
    
    length *= 0.5; // half of container lenght will the considered max distance
    if(length > lengthscale)
    {
        length = lengthscale;
    }
    
    double bin_length = length/ ( nbins );

    bins.resize ( nbins );
    int current_bin;
    double dx, dy;
    double dr;

    for ( size_t i = 0; i < gr.size(); i++ )
    {
        gr[i] = 0.0;
    }

    for ( int i = 0; i < npart; ++i )
    {
        for ( int j = 0; j < npart; ++j )
        {

            if ( i==j )
            {
                continue;
            }

            dx = x[i]-x[j];
            dy = y[i]-y[j];

            minimalImage ( dx, dy, lx, ly );

            dr = sqrt ( dx * dx + dy * dy );
            
            //if(dr < 1)
            //{
            //    std::cout << "problem" << std::endl;
            //}

            current_bin = floor ( dr/bin_length );
            if ( current_bin < nbins )
            {
                gr[current_bin] = gr[current_bin] + 1.0;
            }
        }
    }

    //nomrmalization and output
    double shell_volume;
    for ( int i = 0; i < nbins; ++i )
    {
        shell_volume =  M_PI * ( bin_length * bin_length ) * ( 2 * i + 1 );
        gr[i] /= shell_volume * npart * ( npart/ ( lx*ly ) );
        bins[i] = i * bin_length;
    }
    //std::cout << "done" << std::endl;
    return;
}

void structureFactor ( double* x, double* y, double lx, double ly, double lenghtscale, int n_kavg, int npart, std::vector<double> &sq, std::vector<double> &bins )
{

    //std::cout << "computing structure factor ..." << std::endl;

    // Important Parameters ++++++++++++++++++++++++++++++++++++++
    
    int datapoints = sq.size(); // Number of kbins.
    // n_kavg Number of k to averages per shell (if enought)
    
    double l;
    lx > ly ? l = ly : l = lx; 
    
    double kunit = 2 * M_PI / l;
    
    //max lenghtscale to compute 
    if(l > lenghtscale)
    {
        l = lenghtscale;
    }
    
    double lmax = 1.0 * l; //maximal lenghscale to cover

    int nmax = int ( lmax/kunit ); // Maximum index of ni
    
    //if(nmax > 200)
    //{
    //    nmax = 200;
    //}
    
    //int nmax = 100;
    double kmax = kunit * nmax;
    int lastshell = datapoints;
    double kbinwidth = kmax/datapoints;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    int k_shell_max = 300000;

    int * kx = new int[k_shell_max];
    int * ky = new int[k_shell_max];

    int * kxs = new int[lastshell * n_kavg];
    int * kys = new int[lastshell * n_kavg];

    bins.resize ( datapoints );

    for ( int i = 0; i < lastshell * n_kavg; i++ )
    {
        kxs[i] = -1;
        kys[i] = -1;
    }

    double * results = new double[lastshell]();

    int count = 0;

    //Loop over all the k-shells

    double knormmin, knormmax;

    for ( int shell = 0; shell < lastshell; ++shell ) // Shell loop
    {
        knormmax = ( shell + 1 ) * kbinwidth;
        knormmin = shell * kbinwidth;
        // pick all the kvectors in the shell:
        double knormsq;

        //int n_start = knormmin / kunit;
        int n_stop  = std::ceil(knormmax / kunit);
        
        long k_index = 0;
        long k_index_max = 0;

        for ( int i = 1; i<n_stop; i++ )
        {
            for ( int j = 0; j<n_stop; j++ )
            {
                knormsq = ( i*i+j*j ) * kunit * kunit;
                if ( knormsq <= knormmax*knormmax and knormsq > knormmin * knormmin )
                {
                    kx[k_index] = i;
                    ky[k_index] = j;
                    k_index++;
                    if ( k_index >= k_shell_max )
                    {
                        std::cout << "k_index gerater than k_dim_max, reduce bin size ??" << std::endl;

                        break;
                    }
                }
            }
        }
        k_index_max = k_index;

        // Now in our shell we choose at most nk random kvectors

        int ik;

        std::uniform_int_distribution<int>random_k ( 0,k_index_max-1 );

        for ( int i = 0; i < n_kavg; i++ )
        {
            if ( k_index_max <= n_kavg ) // choose every k vector for small k numbers.
            {
                ik = i;
                if ( i == k_index_max )
                {
                    break;
                }
            }
            else
            {
                ik = random_k (Globalrand::e1);
            }

            kxs[shell*n_kavg+i]=kx[ik];
            kys[shell*n_kavg+i]=ky[ik];
        }
    }
    
    double sqmax = 0.0; // q with the max s(q)
    //int shellmax;
    for ( int shell = 0; shell < lastshell; ++shell ) // Shell loop
    {

        int index;
        count = 0;
        for ( int k = 0 ; k < n_kavg; ++k )
        {

            index = shell * n_kavg + k;
            double tempresult = 0;
            //std::complex<double>tempresult = std::complex<double>(0.0,0.0);

            if ( kxs[index] < 0 )
            {
                break;
            }

            for ( int p = 0; p < npart; p ++ )
            {
                //tempresult += std::complex<double>(cos(kr1),sin(kr1));
                double kr1 = ( x[p] * kxs[index] + y[p] * kys[index] ) * kunit;
                double kr2 = ( x[p] * kys[index] - y[p] * kxs[index] ) * kunit;
                tempresult +=  cos ( kr1 ) + cos ( kr2 );
            }
            //results[shell] += std::norm(tempresult);
            results[shell] += tempresult * tempresult;
            count = k;
        }

        if ( results[shell] > 0.0 ) // avoid datapoints at zero
        {
            results[shell] /= ( count + 1 );

            if ( results[shell] > sqmax )
            {
                sqmax = results[shell];
                //shellmax = shell;
            }

            bins[shell] = shell * kbinwidth;
            sq[shell] = results[shell]/npart;
        }
    }

    delete[] kx;
    delete[] ky;

    delete[] kxs;
    delete[] kys;

    delete[] results;

    //std::cout << "done" << std::endl;
}

}





