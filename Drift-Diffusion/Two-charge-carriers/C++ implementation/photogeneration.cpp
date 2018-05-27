#include <vector>
#include <iostream>
#include <iomanip>
#include "parameters.h"
#include <fstream>
#include <algorithm>

std::vector<double> PhotogenerationRate(){

    std::vector<double> G(num_cell);

    //Uncomment below code for reading in generation rate from file
    /*
     std::ifstream GenRateFile;   

     GenRateFile.open("gen_rate.txt");
     //check if file was opened
     if (!GenRateFile) {
         std::cerr << "Unable to open file gen_rate.txt";
         exit(1);   // call system to stop
     }
     for(int i=1;i<=G.size()-1;i++){
         GenRateFile >> G[i];
         //std::cout << "G(i) " << G[i] <<std::endl;
     }
     double maxOfG = *std::max_element(G.begin(),G.end());

     for(int i=l_HTL_int + 1;i<=G.size()-1;i++){
         G[i] = G_max*G[i]/maxOfG;
         //std::cout << "G(i) " << G[i] <<std::endl;

         GenRateFile.close();
     }
     */

     //Using constant generation rate
     std::fill(G.begin(), G.end(), G_max);


     return G;
}
