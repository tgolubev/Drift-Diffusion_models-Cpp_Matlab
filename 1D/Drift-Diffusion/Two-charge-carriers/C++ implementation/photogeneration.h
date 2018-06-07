#ifndef PHOTOGENERATION_H
#define PHOTOGENERATION_H

#include<vector>
#include "constants.h"
#include "parameters.h"

class Photogeneration
{
public:
    Photogeneration(Parameters &params, double Photogen_scaling);   //constructor will get the generation rate from file
    std::vector<double> getPhotogenRate() {return PhotogenRate;}
private:
    std::vector<double> PhotogenRate;
    double PhotogenRate_max;
};



#endif // PHOTOGENERATION_H
