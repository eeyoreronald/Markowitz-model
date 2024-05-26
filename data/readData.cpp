#include "readData.h"
#include "string_to_double.h"
#include <fstream>
#include <iostream>
#include "csv.h" 

void readData(std::vector<std::vector<double>>& data, const std::string& fileName)
{
    std::ifstream file(fileName);
    Csv csv(file);
    std::string line;
    if (file.is_open())
    {
        int i = 0;
        while (csv.getline(line) != 0) {
            for (int j = 0; j < csv.getnfield(); j++)
            {
                double temp = string_to_double(csv.getfield(j));
                //std::cout << "Asset " << j << ", Return " << i << "=" << temp << "\n";
                data[j][i] = temp;
            }
            i++;
        }
        file.close();

    }
    else {
        std::cout << fileName << " missing\n";
        exit(0);
    }
    return;
}