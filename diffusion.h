/*
 * The library for a multi-group 1D diffusion solver
 */
#include <vector>
#include <cmath>
#include <Eigen/Dense>

class eigenSolut {
    public:
        double K;
        int loopCount;
        Eigen::MatrixXd flux;
        Eigen::MatrixXd fission;

        eigenSolut(double K, int loop, Eigen::MatrixXd flux,Eigen::MatrixXd fission) {
            this->K=K;
            this->loopCount=loop;
            this->flux=flux;
            this->fission=fission;
        }
};

void printArray(const Eigen::MatrixXd&);
bool checkConverg(const Eigen::MatrixXd&,  const Eigen::MatrixXd&,
        const Eigen::MatrixXd&, const Eigen::MatrixXd&, double,int);
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

std::vector<std::vector<std::vector<double>>> 
   buildPropArray(const std::vector<std::vector<std::vector<double>>> &cellSpec, 
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
 * Creates the H matrix for solving H\phi=F\phi+S.
 *
 * Assumes that fluxes are in order by group the space. So 1g1, 1g2, 2g1, ....
 *
 * @param leftBound, the left boundary (0 side) condition. True= vacuum,
 * false=reflective
 * @param righBouund same but for the right side.
 *@return the H matrix assumes first dimension is row. 
 */

Eigen::MatrixXd generateH (std::vector<std::vector<std::vector<double>>> propArray, 
        double delta, bool leftBound, bool rightBound) {
    
    int cells, groups,cellTarget;
    double DtildaMinus, DtildaPlus,sigmaRemoval;

    cells=propArray.size(); //gets the number mesh cells
    groups=propArray[0].size(); //gets the number of groups in model

    Eigen::MatrixXd ret(groups*cells, groups*cells); //creates square return matrix
    for(int i=0; i<cells; i++) { //iterate over all mesh cells
        for(int j=0;j<groups;j++) { //iterate over all groups
            if(i>0) { //if D-tilda-minus exists find it
                DtildaMinus=2*propArray[i][j][0]*propArray[i-1][j][0]/
                    delta*(propArray[i][j][0]+propArray[i-1][j][0]);
            } else if(leftBound) { //if at vacuum boundary update the stufff
                DtildaMinus=2*propArray[i][j][0]/(delta*(1+4*propArray[i][j][0]/delta));
            } else {  //if at reflective boundary set it to 0
                DtildaMinus=0;
            }
            if(i<cells-i) {//same for dTitldaPlus
                DtildaPlus=2*propArray[i][j][0]*propArray[i+1][j][0]/
                    delta*(propArray[i][j][0]+propArray[i+1][j][0]);
            } else if (rightBound) { //if at right vacuum boundary
                DtildaPlus=2*propArray[i][j][0]/(delta*(1+4*propArray[i][j][0]/delta));
            } else { //if at reflective boundary
                DtildaPlus=0;
            }

            cellTarget=groups*i+j; //finds where this cell and group lives in flux array
           //
            /******************Now to start filling the array*************/ 
            //
            if(i>0) { //fill in the left difusion
                //same energy group 1 cell to the left though
                ret(cellTarget,(i-1)*groups+j)= -DtildaMinus;
            } 

            if(i<cells -1 )  {// fill in the right difusion. 
                ret(cellTarget,(i+1)*groups+j)= -DtildaPlus;
            }
            sigmaRemoval=propArray[i][j][1]+propArray[i][j][2]+propArray[i][j][3];

            //fill in flux removal term
            ret(cellTarget,cellTarget)=sigmaRemoval*delta+DtildaMinus+DtildaPlus;
        }
    }
    return ret;
}

/**
 *Generates the fission, F, matrix.
 *
 *Nearly identical params and return values as generateH().
 *
 */

Eigen::MatrixXd generateF(std::vector<std::vector<std::vector<double>>> &propArray, 
              double delta) {
    
    int cells, groups,cellTarget;
    double DtildaMinus, DtildaPlus,chi;
    
    cells=propArray.size(); //get number of cells
    groups=propArray[0].size(); //get number of groups
    
    Eigen::MatrixXd ret(cells*groups,cells*groups);

    for(int i=0; i<cells;i++) { //iterate over all mesh
        for(int j=0;j<groups; j++) { //iterate over all groups
            cellTarget=i*groups+j; //get proper row for this flux
            
            chi=propArray[i][j][4]; //get the chi emission for this group
            //
             /***************************Fill values****************************/
            //
            for(int k=0;k<groups;k++) { //iterate over all groups again. Find all fission sources
                ret(cellTarget,i*groups+k)=chi*propArray[i][k][5]; //chi*nu*Sigma_f,g
            }
            
            //Handle inscatter now from groups above and below
            if(j>0) { //if not group one look for downscatter
                ret(cellTarget,cellTarget-1)+=propArray[i][j-1][3];
            }
            if(j<groups-1) { //if not the lowest thermal group upscatter
                ret(cellTarget,cellTarget+1)+=propArray[i][j+1][2];

            }

        }
    }
    return ret;
}
/**
 *Solves the problem for the primary eigen solution. 
 *
 * Returns the K value and the normalized flux. Flux is "normalized"
 * by dividing it by the average.
 *
 * @param H the H matrix
 * @param F the fision matrix
 * @return an object with K and the flux
 *
 */
eigenSolut solveEigenProblem(const Eigen::MatrixXd &H, const Eigen::MatrixXd &F, int groups) {
    Eigen::MatrixXd flux,histFlux, source,histSource;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp;
    double K,histK;
    bool converged;
    int loopCount;
    K=1.0;

    //initial guesss. Probably really wrong
    flux=Eigen::MatrixXd(H.rows(),1);
    flux.setOnes(); //filll it with 1s
    converged=false; 
    decomp=Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(H); //decompose the problem
    source=flux;
    source.setZero();
    loopCount=0;
    while(!converged&&loopCount<=1000) {
        histSource=source;   //store the old fission source
        source=(1/K)*F*flux; //update the problem source term
        histFlux=flux;   //update old flux

        flux=decomp.solve(source); //solve H\phi=1/K(F\phi)
        
        K*=flux.norm()/histFlux.norm();

        flux=flux/flux.mean(); //renormalize

        std::cout<<"Loop:"<<loopCount<<" K: "<<K<<std::endl;
        converged=checkConverg(flux, histFlux, source,histSource,1e-7,groups);
        loopCount++;
    }
    if(loopCount>=1000) {
        std::cerr<<"Max loop executions reached"<<std::endl;
    }

    return eigenSolut(K,loopCount, flux,source);
}
/**
 *Checks convergence of the fission source.
 *
 * Checks that every term is converged within the tolerance
 *
 * @param flux the current flux
 * @param fluxHist the historic flux
 * @param source the fission source matrix
 * @param histSource the previous iteration's source matrix
 * @param tolerance the fraction of variation that is allowed
 * @return true iff all cells are converged by RMS of flux and fission source
 */

bool checkConverg(const Eigen::MatrixXd &flux, const Eigen::MatrixXd &histFlux,
        const Eigen::MatrixXd &source, const Eigen::MatrixXd &histSource, 
        double tolerance,int groups) {
    double fluxRMS, sourceRMS,fission,histFission;

    for(int i =0; i<source.rows();i++) { //iterate over all rows
        for(int j=0;j<source.cols();j++) { //over all cells
            if(flux(i,j)!=0) { //avoid division by zero
                fluxRMS+=std::pow((flux(i,j)-histFlux(i,j))/flux(i,j),2); //add to RMS
            }
            if(i%groups==0) { //if first energy group sum the fission source
               //it is poluted with cross-group scattering, but meh close
               //enough
               fission=source(i,j); 
               histFission=histSource(i,j);
               if(fission!=0) //avoid divide by 0
                    sourceRMS+=std::pow((fission-histFission)/fission,2);
            }
        }
    }
        //test the convergence condition. Because RMS is soooo fancy.
    if(std::sqrt(fluxRMS/source.rows())<100*tolerance&&
            std::sqrt(sourceRMS/(source.rows()/groups))<tolerance) {
        return true;
    } else {
        return false;
    }
}
/**
 * Prints out the structure of the 2D array to stdout
 */
void printArray(const Eigen::MatrixXd &inspectee) {
    for(int i=0; i<inspectee.rows();i++) {
        std::cout<<'[';
        for(int j=0; j<inspectee.cols();j++) {
            if(inspectee(i,j)==0)
                std::cout<<"-";
            else
                std::cout<<"1";
        }
        std::cout<<']'<<std::endl;
    }
}
