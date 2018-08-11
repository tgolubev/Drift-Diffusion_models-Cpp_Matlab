#include "Utilities.h"
#include "parameters.h"

Utilities::Utilities()
{
}

std::vector<double> Utilities::linear_mix(const Parameters &params, const std::vector<double> &new_values,const std::vector<double> &old_values)
{
    static std::vector<double> result(params.num_cell+1);  //static so only allocate at 1st fnc call

    for (int i = 1; i < params.num_cell; i++) {  //note: all the changing values in vectors start from 1 (0th index is a BC)
        result[i] = new_values[i]*params.w + old_values[i]*(1.0 - params.w);
    }

    return result;
}


void Utilities::write_details(const Parameters &params, double Va, const std::vector<double> &V, const std::vector<double> &p,  const std::vector<double> &n, const std::vector<double> &J_total, const std::vector<double>  &Un, const std::vector<double> &PhotogenRate, const std::vector<double> &R_Langevin)
{
        static std::ofstream VaData;
        //Write charge densities, recombination rates, etc
        std::string filename = std::to_string(Va);
        filename += ".txt";  //add .txt extension
        VaData.open(filename); //this will need to have a string as file name
        for (int i = 1; i <= params.num_cell; i++) {
            VaData << std::setw(15) << std::setprecision(8) << params.dx*i;
            VaData << std::setw(15) << std::setprecision(8) << Vt*V[i];
            VaData << std::setw(15) << std::setprecision(8) << params.N*p[i];         //setprecision(8) sets that use 8 sigfigs
            VaData << std::setw(15) << std::setprecision(8) << params.N*n[i];
            VaData << std::setw(15) << std::setprecision(8) << J_total[i];
            VaData << std::setw(15) << std::setprecision(8) << Un[i];
            VaData << std::setw(15) << std::setprecision(8) << PhotogenRate[i];
            VaData << std::setw(15) << std::setprecision(8) << R_Langevin[i];
            VaData << std::setw(15) << std::setprecision(8) << params.w;
            VaData << std::setw(15) << std::setprecision(8) << params.tolerance <<std::endl;
        }
        VaData.close();
}

void Utilities::write_JV(const Parameters &params, std::ofstream &JV, double iter, double Va, const std::vector<double> &J_total)
{
    if (JV.is_open())
        JV << Va << " " << J_total[static_cast<int>(floor(params.num_cell/2))] << " " << iter << "\n";
}
