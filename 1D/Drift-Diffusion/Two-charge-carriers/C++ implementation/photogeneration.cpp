#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "photogeneration.h"


//Generation rate file should contain num_cell -2 number of entries in a single column, corresponding to
//the the generation rate at each mesh point (except the endpoints).

//constructor definition
Photogeneration::Photogeneration(Parameters &params, double Photogen_scaling){  //requires input of the scaling factor

    PhotogenRate.resize(params.get_num_cell());  //if vector already initialized (in .h), then can't use {} initialization again! Use resize, to give it a size.
    PhotogenRate_max = Photogen_scaling;

    std::ifstream GenRateFile;

     GenRateFile.open("gen_rate.txt");
     //check if file was opened
     if (!GenRateFile) {
         std::cerr << "Unable to open file gen_rate.txt";
         exit(1);   // call system to stop
     }

     for(int i=1;i<=params.get_num_cell()-1;i++){
         GenRateFile >> PhotogenRate[i];
         //std::cout << "G(i) " << G[i] <<std::endl;
     }
     double maxOfGPhotogenRate = *std::max_element(PhotogenRate.begin(),PhotogenRate.end());

     for(int i= 1;i<=params.get_num_cell()-1;i++){
         PhotogenRate[i] = PhotogenRate_max*PhotogenRate[i]/maxOfGPhotogenRate;
         //std::cout << "G(i) " << G[i] <<std::endl;

         GenRateFile.close();
     }


     //Using constant generation rate
     //std::fill(G.begin(), G.end(), G_max);

}
