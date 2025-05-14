#pragma once
#pragma once

#include <memory>
#include <stdint.h>
#include <iostream>

class CellList // cell_nr = ox + oy * nx_;!!
{
public:

    double l_;
    double l_inv_;

    int n_obs_tot_;
    int n_cells_;
    int nx_, ny_;

    int *head_cells_;
    int *list_cells_;
    
    bool sorted_ = false;
    
    //std::unique_ptr<int[]> head_cells_;
    //std::unique_ptr<int[]> list_cells_;

private:

    double cell_length_; // should not be used since the real cell length lx_ and ly_  will differ.

public:

    //constructors and destructors
    CellList()
    {
        l_ = 0.0;
        l_inv_ = 0.0;

        n_obs_tot_ = 0;
        n_cells_ = 0;
        nx_ = 0;
        ny_ = 0;

        head_cells_ = NULL;
        list_cells_ = NULL;
    }

    CellList ( int n_obs_tot, double container_lx, double container_ly, double cell_length )
    {
        // as the simbox contains obstacles outside the box, we add them to the edge cells
        n_obs_tot_ = n_obs_tot;

        // the cells have square shape and the x length of the container is the reference.
        nx_ = container_lx / cell_length;
        ny_ = container_ly / cell_length;
        
        l_ = container_lx / nx_;
        
        if (nx_ <= 0 or ny_ <= 0)
        {
            std::cout << "Pathological cell and container lengh: l_cell >= l_container ... ABORTING" << std::endl;
            std::abort();
        }
        
        n_cells_ = nx_ * ny_;

        l_inv_ = 1.0 / l_;
        
        head_cells_ = new int[n_cells_+1];
        list_cells_ = new int[n_obs_tot_];

        double mem_size = sizeof ( int ) * n_cells_ + sizeof ( int ) * n_obs_tot_;

        //GlobalCounters::memsize += mem_size;
        if ( GlobalFlags::verbose )
        {
            std::cout << std::endl << "---Cell List---" << std::endl;
            printf ( "the celllist hast  %i (%i x %i) cells of dimension %.6f x %.6f \n", n_cells_, nx_, ny_, l_, l_ );
            printf ( "total heapsize of the celllist is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    ~CellList()
    {
        delete[] head_cells_;
        
        if(!sorted_)
        {
            delete[] list_cells_;
        }
    }
   /*
    template<class T_Obstacles>
    T_Obstacles *initialyze ( T_Obstacles &obstacles )
    {
        return fill_sorted( obstacles );
    }
*/
    inline double cellDistanceSquared ( int i, int j )
    {
        //compute the closest distance between central bin (0,0) and bin (i,j)

        double dx, dy;
        //-------------------------------x
        if ( i < 0 )       dx = l_ * ( i + 1 );
        else if ( i == 0 ) dx = 0.0;
        else               dx = l_ * ( i - 1 );

        //-------------------------------y
        if ( j < 0 )       dy = l_ * ( j + 1 );
        else if ( j == 0 ) dy = 0.0;
        else               dy = l_ * ( j - 1 );

        return dx * dx + dy * dy;
    }

    template<typename T_int>
    inline void remapCells ( T_int &i, T_int &j )
    {
        if ( i < 0 )    i += nx_;
        if ( i >= nx_ ) i -= nx_;
        if ( j < 0 )    j += ny_;
        if ( j >= ny_ ) j -= ny_;
    };

    template<typename T_int>
    inline void getCell ( double x, double y, T_int &cellx, T_int &celly )
    {
        cellx = static_cast<T_int> ( x * l_inv_ );
        celly = static_cast<T_int> ( y * l_inv_ );
        
        //cellx = std::floor ( x * l_inv_ );
        //celly = std::floor ( y * l_inv_ );
    }

    template<typename T_int>
    inline void getCellx ( double x, T_int &cellx )
    {
        cellx = static_cast<T_int> ( x * l_inv_ );
    }

    template<typename T_int>
    inline void getCelly ( double y, T_int &celly )
    {
        celly = static_cast<T_int> ( y * l_inv_ );
    }

//private:
    template<class T_Obstacles>
    void fill ( T_Obstacles &obstacles )
    {
        // initializing the linked cellist
        // if -1 is reached as next obstacle then the last element in the cell has been reached.

        for ( auto i = 0; i < n_cells_; i++ )
        {
            head_cells_[i] = -1;
        }

        //putting the obstacles in the cells
        int ox, oy;
        int bin_nr = 0;
        double x, y;

        for ( auto i = 0; i < n_obs_tot_; i++ )
        {
            x = obstacles.x_[i];
            y = obstacles.y_[i];

            // as the simbox contains obstacles outside the box, we add them to the edbge cells
            ox = static_cast<int> ( x * l_inv_ );
            if ( ox < 0 ) ox = 0;
            if ( ox >= nx_ ) ox = nx_ - 1;

            oy = static_cast<int> ( y * l_inv_ );
            if ( oy < 0 ) oy = 0;
            if ( oy >= ny_ ) oy = ny_ - 1;

            bin_nr = ox + oy * nx_;
            list_cells_[i] = head_cells_[bin_nr];
            head_cells_[bin_nr] = i;
        }
    };
    
    /*
    template<class T_Obstacles>
    T_Obstacles* fill_sorted ( T_Obstacles &obstacles )
    {
        T_Obstacles * sorted_obstacles = new T_Obstacles(obstacles);
        // initializing the linked cellist
        // if -1 is reached as next obstacle then the last element in the cell has been reached.
        for ( auto i = 0; i < n_cells_; i++ )
        {
            head_cells_[i] = -1;
        }

        //putting the obstacles in the cells
        int ox, oy;
        int bin_nr = 0;
        double x, y;

        for ( auto i = 0; i < n_obs_tot_; i++ )
        {
            x = obstacles.x_[i];
            y = obstacles.y_[i];

            // as the simbox contains obstacles outside the box, we add them to the edbge cells
            ox = static_cast<int> ( x * l_inv_ );
            if ( ox < 0 ) ox = 0;
            if ( ox >= nx_ ) ox = nx_ - 1;

            oy = static_cast<int> ( y * l_inv_ );
            if ( oy < 0 ) oy = 0;
            if ( oy >= ny_ ) oy = ny_ - 1;

            bin_nr = ox + oy * nx_;
            list_cells_[i] = head_cells_[bin_nr];
            head_cells_[bin_nr] = i;
        }
        
        //now we sort the obstacles according to the cells 
        int j = 0;
        for ( auto i = 0; i < n_cells_; ++i )
        {
            int obs_num = head_cells_[i];
            head_cells_[i] = j;
            
            while ( obs_num >= 0 )
            {
                obstacles.move_to(obs_num, j, *sorted_obstacles);
                ++j;
                obs_num = list_cells_[obs_num];
            }
        }
        head_cells_[n_cells_] = j;
        
        delete[] list_cells_;
        sorted_ = true;
        return sorted_obstacles;
    };
    */
    
};
