#ifndef RUN_DD_H
#define RUN_DD_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>   //allows to use fill and min
#include <fstream>
#include <chrono>
#include <string>
#include <time.h>

#include "constants.h"        //these contain physics constants only
#include "parameters.h"
#include "poisson.h"
#include "continuity_p.h"
#include "continuity_n.h"
#include "recombination.h"
#include "photogeneration.h"
#include "thomas_tridiag_solve.h"
#include "Utilities.h"

std::vector<double> run_DD(Parameters &params);




#endif // RUN_DD_H
