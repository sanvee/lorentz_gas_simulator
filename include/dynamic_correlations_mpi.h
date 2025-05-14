#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>

#include "globals.h"
#include "timescales.h"
#include "propagators.h"
#include "neighbourlist.h"

namespace CorrelationFunctions_Mpi
{
struct vacmsd_results
{
    double *msd_;
    double *vaf_xx_;
    double *vaf_yy_;
    double *vaf_xy_;
    double *vaf_yx_;
    
    double *m4_;
    double *m6_;
    double *m8_;
    double *m10_;
//    double *m12_;
//    double *m14_;
//    double *m16_;

    double *dr_sq_all_;
    int64_t *van_hove_self_count_;
    long double *van_hove_self_;
    long double r_base_;
    long double ln_b_;
    long double rmax_;
    long double rmin_;
    long double ln_rmin_;

    double *fsqt_;

    int *kx_;
    int *ky_;

    std::vector<double> q_;
    std::vector<int> q_index_;

    double q_unit_;
    double qmin_;
    double qmax_;

    int64_t *normalisator_;

    //rbinwdth
    //double *rbinwidth_;

    //--
    int task_;
    int nbint_;
    //number of dr
    int64_t nbindr_;

    //number of bins of the van vove function (time * r_nins)
    int nbinvh_;

    //nbin of the space part of the function
    int nbingr_;

    //max r
    double vhrmax_;

public:

    vacmsd_results ( int nbin )
    {
        task_ = 0;
        msd_ = new double[nbin]();

        vaf_xx_ = new double[nbin]();
        vaf_yy_ = new double[nbin]();
        vaf_xy_ = new double[nbin]();
        vaf_yx_ = new double[nbin]();

        m4_ = new double[nbin]();
        m6_ = NULL;
        m8_ = NULL;
        m10_ = NULL;
        //m12_ = NULL;
        //m14_ = NULL;
        //m16_ = NULL;

        normalisator_ = new int64_t[nbin]();
        nbint_ = nbin;

        dr_sq_all_ = NULL;
        van_hove_self_count_ = NULL;
        //rbinwidth_=NULL;

        fsqt_= NULL;
        kx_ = NULL;
        ky_ = NULL;

        double mem_size = sizeof ( double ) *  nbin * 6.0 + sizeof ( int64_t ) * nbin;
        //GlobalCounters::memsize += mem_size;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Results---" << std::endl;
            printf ( "total heapsize of the results is %.2f kb \n", mem_size / 1000.0 );
        }
    }

    vacmsd_results ( int nbin, TimeScale::LogScale &scale, double rmax, int nbingr )
    {
        task_ = 1;
        //number of timestamps_
        nbindr_ = scale.n_stamps_;
        
        //number of bin per van hove function
        nbingr_ = nbingr;

        //totals bins for functions different times
        nbinvh_ = scale.time_.size() * nbingr;

        nbint_ = nbin;

        msd_ = new double[nbin]();

        vaf_xx_ = new double[nbin]();
        vaf_yy_ = new double[nbin]();
        vaf_xy_ = new double[nbin]();
        vaf_yx_ = new double[nbin]();
        m4_ = new double[nbin]();
        normalisator_ = new int64_t[nbin]();

        dr_sq_all_ = new double[nbindr_]();
        van_hove_self_ = new long double[nbinvh_]();
        van_hove_self_count_ = new int64_t[nbinvh_]();
        
        m6_ = new double[nbin]();
        m8_ = new double[nbin]();
        m10_ = new double[nbin]();
        
        //m12_ = NULL;
        //m14_ = NULL;
        //m16_ = NULL;
        
        fsqt_= NULL;
        kx_ = NULL;
        ky_ = NULL;

        // l propto k * v * tbegin ... v = 1
        //double dr_min = 10 * scale.t_begin_ / nbingr;
        rmax_ = rmax;
        // setting this to 0.1
        rmin_ = 0.1; //scale.t_begin_;
        ln_rmin_ = std::log(rmin_);
        
        //double dr_max = rmax / nbingr;
        //double rbin_base = pow ( dr_max / dr_min, 1.0 / ( nbin - 1 ) ); // intepolate bin width during time

        //r_base_ = std::pow(nbingr_, 1.0/nbingr_); // bin width for log scale histogram at time t
        r_base_ = std::pow ( static_cast<long double>(rmax_/rmin_), static_cast<long double>(1.0/nbingr_)); // bin width for log scale histogram at time t
        ln_b_ = std::log ( r_base_ );

        //for ( int i = 0; i < nbin; i++ )
        //{
            //rbinwidth_[i] = dr_max;
            //rbinwidth_[i] = dr_min * pow ( rbin_base, pow ( i, 2.0 - pow ( ( double ) i / ( ( double ) nbint_ - 1 ), 1. / 3.8 ) ) );
        //}

        double mem_size = sizeof ( double ) *  nbin * 6.0 + sizeof ( int64_t ) * nbin ;

        //Van Hove part:
        double mem_size_vh_dr = sizeof ( double ) * nbindr_;
        double mem_size_vh = sizeof ( double ) *  nbingr_ * nbin * 2;

        //GlobalCounters::memsize += mem_size;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Results with Van-Hove---" << std::endl;
            std::cout << "R max " << rmax_ << std::endl;
            std::cout << "R_min " << rmin_ << std::endl;
            
            //std::cout << "Space resolution small " << rbinwidth_[0] << std::endl;
            //std::cout << "Space resolution big " << rbinwidth_[nbint_ - 1] << std::endl;
            
            std::cout << "Space bins " << nbingr_ << std::endl;
            std::cout << "Time bins " << scale.time_.size() << std::endl;
            std::cout << "Van Hove bins " << nbinvh_ << std::endl;
            std::cout << "Buffer bins " << nbindr_ << std::endl;
            printf ( "heapsize vac and msd: %.2f MB \n", mem_size / 1000000.0 );
            printf ( "heapsize vh_dr_all  %.2f MB \n", mem_size_vh_dr / 1000000.0 );
            printf ( "heapsize vh %.2f MB \n", mem_size_vh / 1000000.0 );
            printf ( "TOTAL %.2f MB \n", ( mem_size + mem_size_vh + mem_size_vh_dr ) / 1000000.0 );
        }
    }

    vacmsd_results ( int nbin, TimeScale::LogScale &scale, std::vector<double> &q, double l, bool modes )
    {
        task_ = 2;
        q_ = q;
        nbindr_ = scale.n_stamps_;

        nbint_ = nbin;

        msd_ = new double[nbin]();
        fsqt_ = new double[nbin * q.size()]();

        vaf_xx_ = new double[nbin]();
        vaf_yy_ = new double[nbin]();
        vaf_xy_ = new double[nbin]();
        vaf_yx_ = new double[nbin]();
        m4_ = new double[nbin]();
        normalisator_ = new int64_t[nbin]();
        
        m6_ = NULL;
        m8_ = NULL;
        m10_ = NULL;
        //m12_ = NULL;
        //m14_ = NULL;
        //m16_ = NULL;

        //rbinwidth_ = NULL;
        dr_sq_all_ = NULL;
        van_hove_self_count_ = NULL;

        double mem_size = sizeof ( double ) *  nbin * ( 6.0 + q.size() ) + sizeof ( int64_t ) * nbin ;

        //Fs(q,t) part:
        // Important Parameters ++++++++++++++++++++++++++++++++++++++
        q_unit_ = 2.0 * M_PI / l;
        
        if ( modes == false )
        {
            for ( uint i = 0; i < q_.size(); ++i )
            {
                q_[i] = std::round(l/q_[i]);
                if (q_[i] <= 0)
                {
                    std::cout << "Warning in Fs(q,t): Lowest mode bigger than container length setting to one." << std::endl;
                    q_[i] = 1;
                }
            }
        }
        //we need to set the modes corresponding to the q:
        
        qmax_ = *std::max_element ( q_.begin(), q_.end() ) * q_unit_;
        qmin_ = *std::min_element ( q_.begin(), q_.end() ) * q_unit_;
        
        //std::cout << "q_min " << qmin_ << std::endl;
        //std::cout << "q_max " << qmax_ << std::endl;
        
        if ( qmin_ < q_unit_ * 0.1 ) //just to avoid carelessness ...
        {
            std::cout << "lowest q to small aborting" << std::endl;
            std::abort();
        }

        int n_qavg = 50; // Number of k to averages per shell (if enought)

        double kbinwidth = q_unit_ / 5.0;

        int nmax = std::ceil ( ( qmax_ + kbinwidth ) / q_unit_); // Maximum index of k_i

        int k_shell_max = 300000;

        int *kx = new int[k_shell_max * static_cast<int> ( q_.size() )];
        int *ky = new int[k_shell_max * static_cast<int> ( q_.size() )];

        std::vector<int>q_index ( q_.size() + 1, 0 );

        //select only first quadrant
        int count = 0;
        //loop over the qs
        //std::ofstream kfile;
        //kfile.open("kfile.txt");
        
        for ( int n = 0; n < static_cast<int> ( q.size() ); ++n )
        {
            double    knormmax = q_[n] * q_unit_ + kbinwidth;
            double    knormmin = q_[n] * q_unit_ - kbinwidth;

            for ( int i = 0; i <= nmax; i++ )
            {
                for ( int j = 1; j <= nmax; j++ )
                {
                    double knormsq = ( i * i + j * j ) * q_unit_ * q_unit_;
                    
                    if ( knormsq <= knormmax *knormmax and knormsq > knormmin * knormmin )  // stupid and slow binning ...
                    {
                        if ( q_index[n + 1] >= k_shell_max )
                        {
                            std::cout << "k_index gerater than k_dim_max, reduce bin size ??" << std::endl;
                            std::abort();
                        }
                        else
                        {
                            //std::cout << i << ' ' << j << std::endl;
                            kx[count] = i;
                            ky[count] = j;
                            //kfile << kx[count] << ' ' << ky[count] << std::endl;
                            ++count;
                            ++q_index[n + 1];
                        }
                    }
                }
            }
        }

        // Now in our shell we choose at most nk random kvectors
        q_index_ = q_index;

        int sum_avg = 0;
        for ( int i = 0; i < static_cast<int> ( q_index_.size() ); ++i )
        {
            if ( q_index_[i] < n_qavg )
            {
                sum_avg += q_index_[i];
            }
            else
            {
                sum_avg += n_qavg;
                q_index_[i] = n_qavg;
            }
        }

        for ( int i = 0; i < static_cast<int> ( q_index_.size() - 1 ); i++ )
        {
            q_index_[i + 1] += q_index_[i];
            q_index[i + 1] += q_index[i];
        }

        std::random_device rd;
        std::mt19937_64 engine ( rd() );
        std::uniform_real_distribution<double> uni_dist ( 0, 1 );

        kx_ = new int[sum_avg];
        ky_ = new int[sum_avg];

        int ik;
        count = 0;
        for ( int n = 0; n < static_cast<int> ( q_.size() ) ; n++ )
        {
            for ( int i = q_index_[n]; i < q_index_[n + 1]; i++ )
            {
                if ( q_index[n + 1] - q_index[n] < n_qavg ) // choose every k vector for small k numbers.
                {
                    ik = i;
                }
                else
                {
                    ik = uni_dist ( engine ) * ( q_index[n + 1] - q_index[n] ) + q_index[n];
                }
                kx_[count] = kx[ik];
                ky_[count] = ky[ik];
                ++count;
            }
        }

        delete[] kx;
        delete[] ky;


        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Results with fsqt---" << std::endl;
            std::cout << "Time bins " << scale.time_.size() << std::endl;
            std::cout << "fsqt: " << std::endl;
            std::cout << "modes: " ;
            for (std::size_t i = 0; i < q_.size(); i++)
            {
                std::cout << q_[i] << ' '; 
            }
            std::cout << std::endl;
            
            std::cout << "wavelenght: " ;
            for (std::size_t i = 0; i < q_.size(); i++)
            {
                std::cout << l/q_[i] << ' '; 
            }
            std::cout << std::endl;
            
            std::cout << "kunit     " << q_unit_ << std::endl;
            std::cout << "kbinwidth " << kbinwidth << std::endl;
            std::cout << "q_max " << qmax_ << std::endl;
            std::cout << "n_max " << nmax << std::endl;
            printf ( "heapsize vac and msd: %.2f MB \n", mem_size / 1000000.0 );
            printf ( "TOTAL %.2f MB \n", ( mem_size / 1000000.0 ) );
        }
    }
    
    vacmsd_results ( int nbin, int task )//yes I know ... we should use named constructors here
    {
        task_ = 3;
        msd_ = new double[nbin]();

        vaf_xx_ = new double[nbin]();
        vaf_yy_ = new double[nbin]();
        vaf_xy_ = new double[nbin]();
        vaf_yx_ = new double[nbin]();

        m4_ = new double[nbin]();
        m6_ = new double[nbin]();
        m8_ = new double[nbin]();
        m10_ = new double[nbin]();
        //m12_ = new double[nbin]();
        //m14_ = new double[nbin]();
        //m16_ = new double[nbin]();

        normalisator_ = new int64_t[nbin]();
        nbint_ = nbin;

        dr_sq_all_ = NULL;
        van_hove_self_count_ = NULL;
        //rbinwidth_=NULL;

        fsqt_= NULL;
        kx_ = NULL;
        ky_ = NULL;

        double mem_size = sizeof ( double ) *  nbin * 6.0 + sizeof ( int64_t ) * nbin;
        //GlobalCounters::memsize += mem_size;

        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Results---" << std::endl;
            printf ( "total heapsize of the results is %.2f kb \n", mem_size / 1000.0 );
        }
    }
    

    ~vacmsd_results()
    {
        delete[] msd_;
        delete[] vaf_xx_;
        delete[] vaf_yy_;
        delete[] vaf_xy_;
        delete[] vaf_yx_;
        delete[] m4_;
        delete[] normalisator_;

        if ( task_ == 1 )
        {
            //delete[] rbinwidth_;
            delete[] dr_sq_all_;
            delete[] van_hove_self_count_;
            delete[] m6_;
            delete[] m8_;
            delete[] m10_;
        }
        if ( task_ == 2 )
        {
            delete[] fsqt_;
            delete[] kx_;
            delete[] ky_;
        }
        if ( task_ == 3 )
        {
            delete[] m6_;
            delete[] m8_;
            delete[] m10_;
            //delete[] m12_;
            //delete[] m14_;
            //delete[] m16_;
        }
    }

    void normalize()
    {
        for ( auto i = 0; i < nbint_ ; i++ )
        {
            long double norm = 1.0 / static_cast<long double>(normalisator_[i]);

            msd_[i] *= norm;
            vaf_xx_[i] *= norm;
            vaf_yy_[i] *= norm;
            vaf_xy_[i] *= norm;
            vaf_yx_[i] *= norm;
            m4_[i] *= norm;

            if ( task_ == 1 )
            {
                m6_[i] *= norm;
                m8_[i] *= norm;
                m10_[i] *= norm;

                van_hove_self_[ i * nbingr_] = static_cast<long double>(van_hove_self_count_[ i * nbingr_]) * ( norm / ( M_PI * rmin_ * rmin_));
                for ( auto j = 0 ; j < nbingr_-1; j++ )
                {
                    //logscale histogram
                    van_hove_self_[i*nbingr_+j+1] = static_cast<long double>(van_hove_self_count_[i*nbingr_+j+1]) * (norm/(M_PI*rmin_*rmin_ * std::pow ( r_base_, 2.0*j ) * ( r_base_ * r_base_ - 1 ))) ;
                }
            }

            if ( task_ == 2 )
            {
                for ( auto j = 0; j < static_cast<int> ( q_.size() ); ++j )
                {
                    double normq = 1.0 / ( q_index_[j+1] - q_index_[j] );
                    fsqt_[nbint_*j + i] *=norm * normq;
                }
            }
            if (task_ == 3)
            {
                m6_[i] *= norm;
                m8_[i] *= norm;
                m10_[i] *= norm;
                //m12_[i] *= norm;
                //m14_[i] *= norm;
                //m16_[i] *= norm;
            }
        }
    }
};

//////////////
/////B!=0/////
//////////////

template <class T_Simbox, class T_Propagator>
void VacMsd_Mpi ( T_Simbox &simbox, Trajectory::Circular &trajectory, T_Propagator &propagator, vacmsd_results &results, TimeScale::LogScale &scale )
{
    //this routine adds the results to the input result no normalization is done.
    int n_averages = scale.n_averages_;
    int_fast64_t n_stamps = scale.n_stamps_;

    //setup calculation buffers
    double *cx_start = new double[n_averages]();
    double *cy_start = new double[n_averages]();
    double *vx_start = new double[n_averages]();
    double *vy_start = new double[n_averages]();

    switch ( results.task_ )
    {

    case 0:
    {
        for ( auto i = 0; i < trajectory.n_trajectories_; ++i )
        {
            trajectory.engage ( i );

            for ( auto stamp = 0; stamp < n_stamps; stamp++ )
            {
                double cx, cy, vx, vy;

                int relative_step = scale.timestamps_[stamp].relative_step_;
                int offset =  scale.timestamps_[stamp].step_offset_;
                double time = scale.timestamps_[stamp].absolute_time_;

                propagator.propagate ( time, *simbox.obstacles_ );

                trajectory.get_current_position_and_velocity ( time, cx, cy, vx, vy );

                if ( relative_step == 0 )
                {
                    trajectory.get_current_position_and_velocity ( time, cx_start[offset], cy_start[offset], vx_start[offset], vy_start[offset] );
                }

                double dx = cx - cx_start[offset];
                double dy = cy - cy_start[offset];

                double dx_2 = dx * dx;
                double dy_2 = dy * dy;
                
                double d_sum = dx_2 + dy_2;

                results.msd_[relative_step] += d_sum;
                
                results.vaf_xx_[relative_step] += dx * vx;
                results.vaf_yy_[relative_step] += dy * vy;
                
                results.vaf_xy_[relative_step] += dx * vy;
                results.vaf_yx_[relative_step] += dy * vx;

                results.m4_[relative_step] += d_sum * d_sum;

                results.normalisator_[relative_step]++;
            }
        }
        break;
    }
    case 1:
    {
        for ( auto i = 0; i < trajectory.n_trajectories_; ++i )
        {
            trajectory.engage ( i );

            for ( auto stamp = 0; stamp < n_stamps; ++stamp )
            {
                double cx, cy, vx, vy;

                int relative_step = scale.timestamps_[stamp].relative_step_;
                int offset =  scale.timestamps_[stamp].step_offset_;
                double time = scale.timestamps_[stamp].absolute_time_;

                propagator.propagate ( time, *simbox.obstacles_ );

                trajectory.get_current_position_and_velocity ( time, cx, cy, vx, vy );
                if ( relative_step == 0 )
                {
                    trajectory.get_current_position_and_velocity ( time, cx_start[offset], cy_start[offset], vx_start[offset], vy_start[offset] );
                }

                double dx = cx - cx_start[offset];
                double dy = cy - cy_start[offset];

                double dx_2 = dx * dx;
                double dy_2 = dy * dy;
                
                double d_sum = dx_2 + dy_2;
                double d_sum_2 = d_sum * d_sum;
                double d_sum_3 = d_sum * d_sum_2;

                results.msd_[relative_step] += d_sum;
                results.vaf_xx_[relative_step] += dx * vx;
                results.vaf_yy_[relative_step] += dy * vy;
                results.vaf_xy_[relative_step] += dx * vy;
                results.vaf_yx_[relative_step] += dy * vx;

                results.m4_[relative_step] += d_sum_2;
                results.m6_[relative_step] += d_sum_3;
                results.m8_[relative_step] += d_sum_2 * d_sum_2;
                results.m10_[relative_step] += d_sum_3 * d_sum_2;

                results.dr_sq_all_[stamp] = d_sum;

                results.normalisator_[relative_step]++;
            }

            //process van-hove
            for ( auto stamp = 0; stamp < n_stamps; stamp++ )
            {
                int relative_step = scale.timestamps_[stamp].relative_step_;
                //linear histogram
                //int bin_number = sqrt ( results.dr_sq_all_[stamp] ) / ( results.rbinwidth_[relative_step] );

                //log scale histogram-----
                int bin_number = 0;
                
                //taking care of the case dr_sq_all_ = 0 (undefined logs)
                if(results.dr_sq_all_[stamp] > 0.0)
                {
                    bin_number = (0.5 * std::log (results.dr_sq_all_[stamp]) - results.ln_rmin_) / results.ln_b_;
                }
                //smallest bin has a hole towards minus infinity :-(
                if (bin_number < 0)
                {
                    bin_number = 0;
                }
                
                if ( bin_number >= results.nbingr_ )
                {
                    bin_number = results.nbingr_ - 1;
                }
                results.van_hove_self_count_[relative_step * results.nbingr_ + bin_number ] += 1;
                //------------------------
            }
        }
        break;
    }
    case 2:
    {
        for ( auto i = 0; i < trajectory.n_trajectories_; ++i )
        {
            trajectory.engage ( i );
            double q_unit = results.q_unit_;
            for ( auto stamp = 0; stamp < n_stamps; ++stamp )
            {
                double cx, cy, vx, vy;

                int relative_step = scale.timestamps_[stamp].relative_step_;
                int offset =  scale.timestamps_[stamp].step_offset_;
                double time = scale.timestamps_[stamp].absolute_time_;

                propagator.propagate ( time, *simbox.obstacles_ );

                trajectory.get_current_position_and_velocity ( time, cx, cy, vx, vy );

                if ( relative_step == 0 )
                {
                    trajectory.get_current_position_and_velocity ( time, cx_start[offset], cy_start[offset], vx_start[offset], vy_start[offset] );
                }

                double dx = cx - cx_start[offset];
                double dy = cy - cy_start[offset];

                double dx_2 = dx * dx;
                double dy_2 = dy * dy;
                
                double d_sum = dx_2 + dy_2;

                results.msd_[relative_step] += d_sum;
                results.vaf_xx_[relative_step] += dx * vx;
                results.vaf_yy_[relative_step] += dy * vy;
                results.vaf_xy_[relative_step] += dx * vy;
                results.vaf_yx_[relative_step] += dy * vx;

                //ngp
                results.m4_[relative_step] += d_sum * d_sum;

                // fsqt
                for ( int qi = 0; qi < static_cast<int> ( results.q_.size() ); ++qi )
                {
                    double expqr_r = 0.0; // adding the also the directions (kx -ky)(-kx -ky) and (-kx-ky)
                    
                    for ( int index = results.q_index_[qi]; index < results.q_index_[qi+1]; ++index )
                    {
                        
                        double kdr =  (dx * results.kx_[index] + dy * results.ky_[index] ) * q_unit;                        
                        double kdr_=  (dx * results.ky_[index] - dy * results.kx_[index] ) * q_unit;

                        expqr_r += 0.5 * (std::cos(kdr) + std::cos(kdr_)); //factor 1/4 average over 4 dierctions
                        
                    }
                    
                    results.fsqt_[qi * results.nbint_ +relative_step] += expqr_r;
                }
                results.normalisator_[relative_step]++;
            }
        }
        break;
    }
    case 3:
    {
        for ( auto i = 0; i < trajectory.n_trajectories_; ++i )
        {
            trajectory.engage ( i );

            for ( auto stamp = 0; stamp < n_stamps; stamp++ )
            {
                double cx, cy;

                int relative_step = scale.timestamps_[stamp].relative_step_;
                int offset =  scale.timestamps_[stamp].step_offset_;
                double time = scale.timestamps_[stamp].absolute_time_;

                propagator.propagate ( time, *simbox.obstacles_ );

                trajectory.get_current_position( time, cx, cy);

                if ( relative_step == 0 )
                {
                    trajectory.get_current_position ( time, cx_start[offset], cy_start[offset]);
                }

                double dx = cx - cx_start[offset];
                double dy = cy - cy_start[offset];

                double dx_2 = dx * dx;
                double dy_2 = dy * dy;
                
                double d_sum = dx_2 + dy_2;
                double d_sum_2 = d_sum * d_sum;
                double d_sum_3 = d_sum * d_sum_2;
                double d_sum_4 = d_sum_2 * d_sum_2;

                results.msd_[relative_step] += d_sum;
                results.m4_[relative_step] += d_sum_2;
                results.m6_[relative_step] += d_sum_3;
                results.m8_[relative_step] += d_sum_4;
                results.m10_[relative_step] += d_sum_4 * d_sum;
                //results.m12_[relative_step] += d_sum_4 * d_sum_2;
                //results.m14_[relative_step] += d_sum_4 * d_sum_3;
                //results.m16_[relative_step] += d_sum_4 * d_sum_4;

                results.normalisator_[relative_step]++;
            }
        }
        break;
    }
    }
    delete[] cx_start;
    delete[] cy_start;
    delete[] vx_start;
    delete[] vy_start;
};


////////////
////B=0/////
////////////

template <class T_Simbox, class T_Propagator>
void VacMsd_Mpi ( T_Simbox &simbox, Trajectory::Linear &trajectory, T_Propagator &propagator, vacmsd_results &results, TimeScale::LogScale &scale )
{
    //this routine adds the results to the input result no normalization is done.
    int n_averages = scale.n_averages_;
    int64_t n_stamps = scale.n_stamps_;

    //setup calculation buffers
    double *cx_start = new double[n_averages]();
    double *cy_start = new double[n_averages]();
    switch ( results.task_ )
    {
    case 0:
    {
        for ( auto i = 0; i < trajectory.n_trajectories_; ++i )
        {
            propagator.engage ( i );
            
            double cx = trajectory.cx_;
            double cy = trajectory.cy_;

            for ( auto stamp = 0; stamp < n_stamps; ++stamp )
            {
                int relative_step = scale.timestamps_[stamp].relative_step_;
                int offset =  scale.timestamps_[stamp].step_offset_;
                double time = scale.timestamps_[stamp].absolute_time_;

                propagator.propagate ( time, simbox.obstacles_ );
                
                //
                //trajectory.get_current_position_inbox ( time, cx, cy );
                //gl_tr_file << cx << ' ' << cy << std::endl;
                //
                
                trajectory.get_current_position ( time, cx, cy );
                
                if ( relative_step == 0 )
                {
                    trajectory.get_current_position ( time, cx_start[offset], cy_start[offset] );
                }
                
                

                double dx = cx - cx_start[offset];
                double dy = cy - cy_start[offset];
                
                double dx_2 = dx*dx;
                double dy_2 = dy*dy;
                
                double d_sum = dx_2 + dy_2;

                results.msd_[relative_step] += d_sum;
                results.m4_[relative_step] += d_sum * d_sum;

                results.normalisator_[relative_step]++;
            }
        }
        break;
    }
    case 1:
    {
        for ( auto i = 0; i < trajectory.n_trajectories_; ++i )
        {
            propagator.engage ( i );

            double cx = trajectory.cx_;
            double cy = trajectory.cy_;

            for ( auto stamp = 0; stamp < n_stamps; ++stamp )
            {
                int relative_step = scale.timestamps_[stamp].relative_step_;
                int offset =  scale.timestamps_[stamp].step_offset_;
                double time = scale.timestamps_[stamp].absolute_time_;

                propagator.propagate ( time, simbox.obstacles_ );

                //
                //trajectory.get_current_position_inbox ( time, cx, cy );
                //gl_tr_file << cx << ' ' << cy << std::endl;
                //
                
                trajectory.get_current_position ( time, cx, cy );
                
                if ( relative_step == 0 )
                {
                    trajectory.get_current_position ( time, cx_start[offset], cy_start[offset] );
                }

                double dx = cx - cx_start[offset];
                double dy = cy - cy_start[offset];
                
                double dx_2 = dx*dx;
                double dy_2 = dy*dy;
                
                double d_sum = dx_2 + dy_2;
                double d_sum_2 = d_sum * d_sum;
                double d_sum_3 = d_sum * d_sum_2;

                //msd
                results.msd_[relative_step] += d_sum;

                //ngp
                results.m4_[relative_step] += d_sum_2;
                results.m6_[relative_step] += d_sum_3;
                results.m8_[relative_step] += d_sum_2 * d_sum_2;
                results.m10_[relative_step] += d_sum_3 * d_sum_2;
                

                //van_hove
                results.dr_sq_all_[stamp] = d_sum;

                results.normalisator_[relative_step]++;
            }

            //process van-hove
            for ( auto stamp = 0; stamp < n_stamps; ++stamp )
            {
                int relative_step = scale.timestamps_[stamp].relative_step_;

                //log scale histogram-----
                int bin_number = 0;
                
                //taking care of the case dr_sq_all_ = 0 (undefined logs)
                if(results.dr_sq_all_[stamp] > 0.0)
                {
                    bin_number = (0.5 * std::log (results.dr_sq_all_[stamp]) - results.ln_rmin_) / results.ln_b_;
                }
                //smallest bin has a hole towards minus infinity
                if (bin_number < 0)
                {
                    bin_number = 0;
                }
                
                if ( bin_number >= results.nbingr_ )
                {
                    bin_number = results.nbingr_ - 1;
                }
                results.van_hove_self_count_[relative_step * results.nbingr_ + bin_number ] += 1;
                //------------------------
            }
        }
        break;
    }
    case 2:
    {
        for ( auto i = 0; i < trajectory.n_trajectories_; ++i )
        {
            propagator.engage ( i );

            double cx = trajectory.cx_;
            double cy = trajectory.cy_;
            double q_unit = results.q_unit_;

            for ( auto stamp = 0; stamp < n_stamps; ++stamp )
            {
                int relative_step = scale.timestamps_[stamp].relative_step_;
                int offset =  scale.timestamps_[stamp].step_offset_;
                double time = scale.timestamps_[stamp].absolute_time_;

                propagator.propagate ( time, simbox.obstacles_ );

                trajectory.get_current_position ( time, cx, cy );
                if ( relative_step == 0 )
                {
                    trajectory.get_current_position ( time, cx_start[offset], cy_start[offset] );
                }

                double dx = cx - cx_start[offset];
                double dy = cy - cy_start[offset];

                double dx_2 = dx * dx;
                double dy_2 = dy * dy;
                
                double d_sum = dx_2 + dy_2;

                //msd
                results.msd_[relative_step] += d_sum;

                //ngp
                results.m4_[relative_step] += d_sum * d_sum;

                results.normalisator_[relative_step]++;

                // fsqt
                for ( int qi = 0; qi < static_cast<int> ( results.q_.size() ); ++qi )
                {
                    double expqr_r = 0.0;
                    //double expqr_i = 0.0; is zero
                    
                    for ( int index = results.q_index_[qi]; index < results.q_index_[qi+1]; ++index )
                    {
                        
                        double kdr =  (dx * results.kx_[index] + dy * results.ky_[index] ) * q_unit;                        
                        double kdr_=  (dx * results.ky_[index] - dy * results.kx_[index] ) * q_unit;

                        expqr_r += 0.5 * (std::cos(kdr) + std::cos(kdr_)); //factor 1/4 average over 4 dierctions
                    }

                    results.fsqt_[qi * results.nbint_ + relative_step] += expqr_r;
                }
            }
        }
        break;
    }
    case 3:
    {
        for ( auto i = 0; i < trajectory.n_trajectories_; ++i )
        {
            propagator.engage ( i );

            for ( auto stamp = 0; stamp < n_stamps; stamp++ )
            {
                double cx, cy;

                int relative_step = scale.timestamps_[stamp].relative_step_;
                int offset =  scale.timestamps_[stamp].step_offset_;
                double time = scale.timestamps_[stamp].absolute_time_;

                propagator.propagate ( time, simbox.obstacles_ );

                trajectory.get_current_position( time, cx, cy);

                if ( relative_step == 0 )
                {
                    trajectory.get_current_position ( time, cx_start[offset], cy_start[offset]);
                }

                double dx = cx - cx_start[offset];
                double dy = cy - cy_start[offset];

                double dx_2 = dx * dx;
                double dy_2 = dy * dy;
                
                double d_sum = dx_2 + dy_2;
                double d_sum_2 = d_sum * d_sum;
                double d_sum_3 = d_sum * d_sum_2;
                double d_sum_4 = d_sum_2 * d_sum_2;

                results.msd_[relative_step] += d_sum;
                results.m4_[relative_step] += d_sum_2;
                results.m6_[relative_step] += d_sum_3;
                results.m8_[relative_step] += d_sum_2 * d_sum_2;
                results.m10_[relative_step] += d_sum_4 * d_sum;
                //results.m12_[relative_step] += d_sum_4 * d_sum_2;
                //results.m14_[relative_step] += d_sum_4 * d_sum_3;
                //results.m16_[relative_step] += d_sum_4 * d_sum_4;

                results.normalisator_[relative_step]++;
            }
        }
        break;
    }
    }
    delete[] cx_start;
    delete[] cy_start;
};
};


