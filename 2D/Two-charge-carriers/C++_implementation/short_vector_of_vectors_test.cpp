// Example program
#include <iostream>
#include <string>
#include <vector>


int main()
{

std::vector<std::vector<double> > matrix;

matrix.resize(10, std::vector<double>(10));

std::cout << matrix[9].size() << std::endl;

 for (int i = 0; i < 10; i++) {
        std::fill(matrix[i].begin(), matrix[i].end(),5);
    }
    
    std::cout << matrix[8][4] << std::endl;
}
