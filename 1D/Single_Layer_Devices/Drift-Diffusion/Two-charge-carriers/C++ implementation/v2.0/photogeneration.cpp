#include <vector>
#include <iostream>

#include "photogeneration.h"

//constructor definition
Photogeneration::Photogeneration(const Parameters &params, double photogen_scaling, const std::string gen_rate_file_name){

    PhotogenRate.resize(params.num_cell);
    PhotogenRate_max = photogen_scaling;


    std::ifstream GenRateFile;

     GenRateFile.open(gen_rate_file_name);
     //check if file was opened
     if (!GenRateFile) {
         std::cerr << "Unable to open file gen_rate.txt";
         exit(1);   // call system to stop
     }

     for (int i = 1; i <= params.num_cell-1; i++) {
         GenRateFile >> PhotogenRate[i];
         //std::cout << "G(i) " << G[i] <<std::endl;
     }
     double maxOfGPhotogenRate = *std::max_element(PhotogenRate.begin(),PhotogenRate.end());

     for (int i= 1; i <= params.num_cell-1; i++) {
         PhotogenRate[i] = PhotogenRate_max*PhotogenRate[i]/maxOfGPhotogenRate;
         //std::cout << "G(i) " << G[i] <<std::endl;

         GenRateFile.close();
     }



     //Using constant generation rate
     //std::fill(PhotogenRate.begin(), PhotogenRate.end(), G_max);

}
