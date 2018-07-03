#include "Utilities.h"
#include "parameters.h"

Utilities::Utilities()
{
}

std::vector<double> Utilities::linear_mix(const Parameters &params, const std::vector<double> &new_values,const std::vector<double> &old_values)
{
    static std::vector<double> result(params.num_elements+1);  //static so only allocate at 1st fnc call

    for (int i = 1; i <= params.num_elements; i++) {  //note: all the changing values in vectors start from 1 (0th index is a BC)
        result[i] = new_values[i]*params.w + old_values[i]*(1.0 - params.w);
    }

    return result;
}


void Utilities::write_details(const Parameters &params, double Va, const Eigen::Tensor<double, 3> &V_matrix, const std::vector<double> &p, const std::vector<double> &n, const Eigen::Tensor<double, 3> &J_total_Z, const std::vector<double>  &Un)
{

    //CAN ADD AN n and p to matrix conversion inside here, b/c that's the only place that it's used.


        static std::ofstream VaData;
        //Write charge densities, recombination rates, etc
        std::string filename = std::to_string(Va);
        filename += ".txt";  //add .txt extension
        VaData.open(filename); //this will need to have a string as file name
        //for now, write out only information for a line profile along the z direction, and use middle of the device in x direction.
        for (int k = 1; k < params.num_cell; k++) {
            int i =  static_cast<int>(floor(params.num_cell/2));  //just use middle index for now
            int j =  static_cast<int>(floor(params.num_cell/2));
            VaData << std::setw(15) << std::setprecision(8) << params.dx*i;
            VaData << std::setw(15) << std::setprecision(8) << params.dx*j;
            VaData << std::setw(15) << std::setprecision(8) << params.dx*k;
            VaData << std::setw(15) << std::setprecision(8) << Vt*V_matrix(i,j,k);
            //VaData << std::setw(15) << std::setprecision(8) << params.N_dos*p_matrix(i,j);
            //VaData << std::setw(15) << std::setprecision(8) << params.N_dos*n_matrix(i,j);
            VaData << std::setw(15) << std::setprecision(8) << J_total_Z(i,j,k);
            //VaData << std::setw(15) << std::setprecision(8) << Un_matrix(i,j);
            VaData << std::setw(15) << std::setprecision(8) << params.w;
            VaData << std::setw(15) << std::setprecision(8) << params.tolerance <<std::endl;
        }
        VaData.close();
}

void Utilities::write_JV(const Parameters &params, std::ofstream &JV, double iter, double Va, const Eigen::Tensor<double, 3> &J_total_Z)
{
    if (JV.is_open()) {
        int i =  static_cast<int>(floor(params.num_cell/2));
        JV << Va << " " << J_total_Z(i,i,i) << " " << iter << "\n";
    }
}
