#pragma once

#include <fstream>
#include <string>

#include "obstacles.h"

namespace IO
{

//wirtes a header with informations on the simulation
//contains 10 lines per default
void writeHeader ( std::ofstream &outfile )
{
    outfile << "" << std::endl;
    outfile << "" << std::endl;
    outfile << "" << std::endl;
    outfile << "" << std::endl;
    outfile << "" << std::endl;
    outfile << "" << std::endl;
    outfile << "" << std::endl;
    outfile << "" << std::endl;
    outfile << "" << std::endl;
    outfile << "" << std::endl;
}

void makeFilenameUnique ( std::string &name, std::string ending )
{
    std::string prefix = name;
    std::string filename = prefix + ending;
    int flag = 1;

    while ( std::ifstream ( filename ) )
    {
        prefix = name + "_(" + std::to_string ( flag ) + ")";
        filename =  prefix + ending;
        flag ++;
    }
    name = prefix;
}
};
