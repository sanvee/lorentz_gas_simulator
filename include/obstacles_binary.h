#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

#include "globals.h"

namespace Obstacles
{

template<class T_o1, class T_o2>
class Binary
{
public:

    int n_tot1_;
    int n_tot_;
    int n_obstacles_;
    int sides1_;
    int sides2_;
    int sides_;

    double max_radius_;

    int type1_;
    int type2_;

    std::vector<double> x_;    // center position
    std::vector<double> y_;

    //int *type_;

    T_o1 *obstacles1_;
    T_o2 *obstacles2_;

    Binary()
    {

    }

    Binary ( T_o1 *obstacles1, T_o2 *obstacles2 )
    {
        obstacles1_ = obstacles1;
        obstacles2_ = obstacles2;

        obstacles1->max_radius_ > obstacles2->max_radius_ ? max_radius_ = obstacles1->max_radius_ : max_radius_ = obstacles2->max_radius_;

        sides1_ = obstacles1_->sides_;
        sides2_ = obstacles2_->sides_;

        sides1_ > sides2_ ? sides_ = sides1_: sides_ = sides2_;

        n_tot_ = obstacles1_->n_tot_ + obstacles2_->n_tot_;
        n_tot1_ = obstacles1->n_tot_;
        n_obstacles_ = obstacles1_->n_obstacles_ + obstacles2_->n_obstacles_;

        x_.resize ( n_tot_ );
        y_.resize ( n_tot_ );
        
        type1_ = obstacles1->type_;
        type2_ = obstacles2->type_;

        for ( int i = 0; i < obstacles1_->n_tot_; ++i )
        {
            x_[i]=obstacles1->x_[i];
            y_[i]=obstacles1->y_[i];
        }

        for ( int i = 0; i < obstacles2_->n_tot_; ++i )
        {
            x_[i+n_tot1_]=obstacles2->x_[i];
            y_[i+n_tot1_]=obstacles2->y_[i];
        }

        if ( GlobalFlags::verbose )
        {
            double mem_size = sizeof ( double ) * n_tot_ * 2.0 + sizeof ( int ) * n_tot_;
            std::cout << std::endl << "---Obstacles (Binary)---" << std::endl;
            std::cout << "n Obstacles " << n_obstacles_ << std::endl;

            printf ( "total heapsize of the Obstacles is %.2f kB \n", mem_size / 1000.0 );
        }
    }

    inline int getType ( int i )
    {
        if ( i < n_tot1_ )
        {
            return type1_;
        }
        else
        {
            return type2_;
        }
    }
    
    inline int getRealNumber( int i)
    {
        if ( i < n_tot1_ )
        {
            return i;
        }
        else
        {
            return i - n_tot1_;
        }
    }

    int getpart ( int i )
    {
        if ( i < n_tot1_ )
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }

    bool isInside ( int *obstacles_list, int last, double x, double y )
    {
        for ( auto i = 0; i < last; i++ )
        {
            int n = obstacles_list[i];
            if ( n < n_tot1_ )
            {
                if ( obstacles1_->isInside ( n, x, y ) )
                {
                    return true;
                }
            }
            else
            {
                if ( obstacles2_->isInside ( n - n_tot1_, x, y ) )
                {
                    return true;
                }
            }
        }
        return false;
    }
    
    void writeObstacles ( std::string filename)
    {
        obstacles1_->writeObstacles(filename, false);
        obstacles2_->writeObstacles(filename, true);
    }
    
};



//template<class T_o1, class T_o2>
//Binary(T_o1 o1, T_o2 o2) -> Binary<typename T_o1::value_type,typename T_o2::value_type>;


};


