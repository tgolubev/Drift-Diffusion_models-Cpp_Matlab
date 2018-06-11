#ifndef PHOTOGENERATION_H
#define PHOTOGENERATION_H

#include<vector>
#include "constants.h"
#include "parameters.h"
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>

class Photogeneration
{
public:
    Photogeneration(Parameters &params, double Photogen_scaling, std::string GenRateFileName);   //constructor will get the generation rate from file
    std::vector<double> getPhotogenRate() {return PhotogenRate;}
private:
    std::vector<double> PhotogenRate;
    double PhotogenRate_max;
};



#endif // PHOTOGENERATION_H
