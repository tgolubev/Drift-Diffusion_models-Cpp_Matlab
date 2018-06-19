#ifndef PHOTOGENERATION_H
#define PHOTOGENERATION_H

#include<vector>
#include "constants.h"
#include "parameters.h"
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>

#include <Eigen/Dense>

class Photogeneration
{
public:

    //!Constructor will get the generation rate from file.
    //!Generation rate file should contain num_cell -2 number of entries in a single column, corresponding to
    //!the the generation rate at each mesh point (except the endpoints).
    //! \param photogen_scaling is the scaling factor obtained from fit to get the correct short-circuit current.
    Photogeneration(const Parameters &params, double photogen_scaling, const std::string gen_rate_file_name);

    Eigen::MatrixXd getPhotogenRate() {return PhotogenRate;}

private:
    Eigen::MatrixXd PhotogenRate;
    double PhotogenRate_max;
};



#endif // PHOTOGENERATION_H
