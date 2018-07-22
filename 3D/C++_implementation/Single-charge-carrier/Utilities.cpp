#include "Utilities.h"
#include "parameters.h"

Utilities::Utilities()
{
}

Eigen::Tensor<double, 3> Utilities::linear_mix(const Parameters &params, const Eigen::Tensor<double, 3> &new_values, const Eigen::Tensor<double, 3> &old_values)
{
    static Eigen::Tensor<double, 3> result(params.num_elements, 1, 1);  //static so only allocate at 1st fnc call

    for (int i = 0; i < params.num_elements; i++) {
        result(i,0,0) = new_values(i,0,0)*params.w + old_values(i,0,0)*(1.0 - params.w);  //since now dealing with Tensor's, neeed to use (), otherwise will get "YOU MADE A PROGRAMMING MISTAKE" ERROR, about brackets....
    }

    return result;
}


void Utilities::write_details(const Parameters &params, double Va, const Eigen::Tensor<double, 3> &V_matrix, Eigen::Tensor<double, 3> &p_matrix, const Eigen::Tensor<double, 3> &J_total_Z, const std::vector<double>  &Up)
{

    //CAN ADD AN n and p to matrix conversion inside here, b/c that's the only place that it's used.


        static std::ofstream VaData;
        //Write charge densities, recombination rates, etc
        std::string filename = std::to_string(Va);
        filename += ".txt";  //add .txt extension
        VaData.open(filename); //this will need to have a string as file name
        //for now, write out only information for a line profile along the z direction, and use middle of the device in x direction.
        for (int k = 1; k < params.num_cell_z; k++) {
            int i =  static_cast<int>(floor(params.num_cell_x/2));  //just use middle index for now
            int j =  static_cast<int>(floor(params.num_cell_y/2));
            VaData << std::setw(15) << std::setprecision(8) << params.dx*i;
            VaData << std::setw(15) << std::setprecision(8) << params.dy*j;
            VaData << std::setw(15) << std::setprecision(8) << params.dz*k;
            VaData << std::setw(15) << std::setprecision(8) << Vt*V_matrix(i,j,k);
            VaData << std::setw(15) << std::setprecision(8) << params.N_dos*p_matrix(i,j,k);
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
        int i =  static_cast<int>(floor(params.num_cell_x/2));
        int j =  static_cast<int>(floor(params.num_cell_y/2));
        int k =  static_cast<int>(floor(params.num_cell_z/2));
        JV << Va << " " << J_total_Z(i,i,i) << " " << iter << "\n";
    }
}
