/*
 * The library for a multi-group 1D diffusion solver
 */
#include <vector>

/**
 * Sums all of the elements of a 1D array
 *
 * @param summee the array to be summed.
 * @return the sum of the entire array
 */
double sumArray(const std::vector<double> &summee) {
    double sum=0;
    for(double buff: summee) {
        sum+=buff;
    }
    return sum;
}
int sumArray(const std::vector<int> &summee) {
    int sum=0;
    for(int buff:summee) {
        sum+=buff;
    }
    return sum;
}

/**
 *Builds 2D vector of material properties for every mesh point.
 *
 *This vector can then be fed into the matrix builder functions.
 *
 *@param cellSpec, the material properties in each physical cell level 1-cell
                    level 2- groups
                    level 3 properties
                        0-D
                        1-Sigma_a
                        2-Sigma_upscat
                        3-Sigma_downscat
                        4-chi
                        5-nu-Sigma_f
 *@param cellThick the thickness in unit mesh cells for each spacial region above
 *@param delta     the thickness in cm of each mesh grid
 *
 *@return the top level of the vector is each material the lower level mimicks cellSpec
 *
 */

std::vector<std::vector<std::vector<double>>> buildPropArray(
        const std::vector<std::vector<std::vector<double>>> &cellSpec, 
        const std::vector<int> &cellThick) {
    int total,groups,outI;
     
    total=sumArray(cellThick);
    groups=cellSpec[0].size(); //get the number of groups being used
    std::vector<std::vector<std::vector<double>>> ret
        (total,std::vector<std::vector<double>>(groups,
            std::vector<double>(6))); //init the return vector

    outI=0;
    for(int i =0; i<cellSpec.size(); i++) { //iterate over all cells
        for(int j=0;j<cellThick[i];j++) { //count over all of the meshpoints in a cell
            ret[outI]=cellSpec[i];  //copies in the material properties for that cell
            outI++; //move onto the next mesh point
        }
    }
    return ret;
}

/**
 *
 *
 *
 *
 */

std::vector<std::vector<double>> generateH (std::vector<
        std::vector<std::vector<double>>> propArray, double delta) {

}
