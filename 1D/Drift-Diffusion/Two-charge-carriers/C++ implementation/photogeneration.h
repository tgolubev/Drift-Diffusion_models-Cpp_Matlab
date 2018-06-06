#ifndef PHOTOGENERATION_H
#define PHOTOGENERATION_H

#include<vector>
#include "parameters.h"
#include "simulation.h"

class Photogeneration
{
public:
    Photogeneration(Simulation &simul, double Photogen_scaling);   //constructor will get the generation rate from file
    std::vector<double> getPhotogenRate() {return PhotogenRate;}
private:
    std::vector<double> PhotogenRate;
    double PhotogenRate_max;
};



#endif // PHOTOGENERATION_H
