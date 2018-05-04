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
 *@param cellSpec, the material properties in each physical cell
 *@param cellThick the thickness in unit mesh cells for each spacial region above
 *@param delta     the thickness in cm of each mesh grid
 *
 *@return the top level of the vector is each material the lower level mimicks cellSpec
 *
 */

std::vector<std::vector<double>> buildPropArray(
        const std::vector<std::vector<double>> &cellSPec, 
        const std::vector<int> &cellThick, 
        double delta) {
    int total;
     
    total=sumArray(cellThick);
    std::vector<std::vector<double>> ret(total,std::vector<double>(5)); //init the return vector

} 
