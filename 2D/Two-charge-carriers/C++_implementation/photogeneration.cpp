#include <vector>
#include <iostream>

#include "photogeneration.h"

//constructor definition
Photogeneration::Photogeneration(const Parameters &params, double photogen_scaling, const std::string gen_rate_file_name){

    PhotogenRate.resize(params.num_cell, params.num_cell);  //need 0 through N indices, and N = num_cell-1
    PhotogenRate_max = photogen_scaling;

    std::vector<double> Photogen_vector; //temporary vector, to input gen rate from file

    std::ifstream GenRateFile;

     GenRateFile.open(gen_rate_file_name);
     //check if file was opened
     if (!GenRateFile) {
         std::cerr << "Unable to open file gen_rate.txt";
         exit(1);   // call system to stop
     }

     for (int i = 1; i <= params.num_cell-1; i++) {
         GenRateFile >> Photogen_vector[i];
     }

     //----normalize the photogen rate--------------------------

     double maxOfGPhotogenRate = *std::max_element(Photogen_vector.begin(),Photogen_vector.end());

     for (int i= 1; i <= params.num_cell-1; i++) {
         Photogen_vector[i] = PhotogenRate_max*Photogen_vector[i]/maxOfGPhotogenRate;
         //std::cout << "G(i) " << G[i] <<std::endl;
     }

     //-------------------------------------
     //fill PhotogenRate matrix
     for (int i = 1; i <= params.num_cell - 1; i++)
         for (int j = 1; j <= params.num_cell - 1; j++)
             PhotogenRate(i,j) = Photogen_vector[i];

     GenRateFile.close();

     //Using constant generation rate
     //std::fill(G.begin(), G.end(), G_max);

}
