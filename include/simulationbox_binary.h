#pragma once

#include "container.h"
#include "config.h"
#include "obstacles.h"
#include "globals.h"
#include "cell_list.h"
#include "obstacles.h"
#include "mc_obstacles.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

namespace SimulationBox // sets up a Periodic box with ghost obstacles
{

template <class T_box1, class T_box2, class T_obs1, class T_obs2>
class BinaryBox
{
public:
    int nr_ghosts_approx_;
    int nr_total_approx_;

    Container *container_;
    CellList *celllist_;

    Obstacles::Binary<T_obs1, T_obs2> *obstacles_;

public:
    ~BinaryBox()
    {
        delete container_;
        delete celllist_;
        delete obstacles_;
    };

    BinaryBox()
    {
        nr_ghosts_approx_ = 0;
        nr_total_approx_ = 0;

        container_ = NULL;
        celllist_ = NULL;
        obstacles_ = NULL;
    };

public:

    BinaryBox ( T_box1 *simbox1, T_box2 *simbox2, T_obs1 *obstacles1, T_obs2 *obstacles2, double target_cell_length )
    {
        if ( simbox1->container_->lx_!=simbox2->container_->lx_ or simbox1->container_->ly_!=simbox2->container_->ly_ )
        {
            std::cout << "Error Container sizes do not match in Binary obstacles" << std::endl;
            std::abort();
        }

        obstacles_ = new Obstacles::Binary<T_obs1, T_obs2> ( simbox1->obstacles_, simbox2->obstacles_ );
        container_ = new Container ( simbox1->container_->lx_, simbox1->container_->ly_ );
        celllist_ = new CellList ( obstacles_->n_tot_, container_->lx_, container_->ly_, target_cell_length );
        celllist_->fill ( *obstacles_ );
    }

    void writeObstacles_config ( std::string folder, std::string filename )
    {
        mkdir ( folder.c_str(), 0755 );
        std::string sub_folder = folder + "/configurations";
        mkdir ( sub_folder.c_str(), 0755 );
        obstacles_->writeObstacles ( sub_folder + '/' + filename );
    }

//template <class T_box1, class T_box2, class T_obs1, class T_obs2>
//BinaryBox( T_box1 *simbox1, T_box2 *simbox2) -> BinaryBox <decltype(simbox1), decltype(simbox2), decltype(simbox1->obstacles_),decltype(simbox2->obstacles_)>;
};
}
