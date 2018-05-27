#include <vector>
#include <iostream>
#include <iomanip>
#include "parameters.h"
#include <fstream>
#include <algorithm>

std::vector<double> PhotogenerationRate(){
     std::ifstream GenRateFile;
     std::vector<double> G(num_cell);

     GenRateFile.open("gen_rate.txt");
     //check if file was opened
     if (!GenRateFile) {
         std::cerr << "Unable to open file gen_rate.txt";
         exit(1);   // call system to stop
     }
     for(int i=l_HTL_int + 1;i<=G.size()-1;i++){
         GenRateFile >> G[i];
         //std::cout << "G(i) " << G[i] <<std::endl;
     }
     double maxOfG = *std::max_element(G.begin(),G.end());
     for(int i=l_HTL_int + 1;i<=G.size()-1;i++){
         G[i] = G_max*G[i]/maxOfG;
         //std::cout << "G(i) " << G[i] <<std::endl;
     }
     GenRateFile.close();
     //std::fill(G.begin(), G.end(), G_max);  //for now use a constant generation rate
     return G;
}
