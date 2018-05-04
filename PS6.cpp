/**
 * Code for 22.211 Pset 6.
 * Copyright 2018 Micah Gale.
 */

#include <iostream>
#include "diffusion.h"

std::vector<std::vector<double>> getMaterials(int); 

int main () {
    std::vector<double> test{1,2,3,4,5};
    std::cout<<sumArray(test);
    return 0;

}

/**
 *Returns an array of material properties based on the problem specifed mat type
 *
 * @param requested. The index of the material per the PS6 problem statement
 * @return the material data formatted in a way that the library expects
 */
std::vector<std::vector<double>> getMaterials(int requested) {
    std::vector<std::vector<double>> ret(2,std::vector<double>(5,0));

    switch(requested)  {
        case 1:
            ret={{1.43,0.0079,0,0.0195, 1,0.0034},
                {0.37, 0.0605,0 , 0, 0, 0.0711}};
            break;
        case 2: 
            ret={{1.43,0.0084,0,0.0185,1,0.0054},
                {0.37,0.0741, 0, 0,    0,0.1}};
            break;
        case 3:
            ret={{1.43,0.0089,0 ,0.0178, 1, 0.0054},
                {0.37,0.0862, 0, 0,      0, 0.1}};
            break;
        case 4:
            ret={{1.43,0.0088, 0, 0.0188, 1, 0.0062},
                {0.37, 0.0852, 0, 0,      0, 0.1249}};

            break;
        case 5:
            ret={{1.26, 0.0025, 0, 0.0294, 0, 0},
                {0.27, 0.02,    0, 0,      0, 0}};
            break;
        case 6:
            ret={{1.0, 0.0054, 0, 0.0009, 0, 0},
                {0.34, 0.13,   0, 0,      0, 0}};

            break;
        case 7:
            ret={{1.55, 0.001, 0, 0.045, 0, 0},
                {0.27,  0.0286, 0, 0,    0, 0}};
            break;

    }

}
