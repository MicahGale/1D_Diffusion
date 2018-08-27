/**
 * Code for 22.211 Pset 6.
 * Copyright 2018 Micah Gale.
 */

#include <iostream>
#include "diffusion.h"
#include <fstream>

using namespace Eigen;
std::vector<std::vector<double>> getMaterials(int); 
eigenSolut Q1();
eigenSolut Q2();
eigenSolut Q3();
eigenSolut Q4();
eigenSolut Q5(int);
void parseResults(const eigenSolut&,std::ofstream&,double);

//runs all problems and creates the table
void partA() {
    std::ofstream csv;

    csv.open("PartA.csv");
    csv<<"Problem,$K_{eff}$,Peak Fission, Peak fission Location, Number of Iterations"<<std::endl;
    csv<<"1,";
    parseResults(Q1(),csv,5.0);
    parseResults(Q2(),csv,5.0);
    parseResults(Q3(),csv,1.0);
    parseResults(Q4(),csv,1.0);
    parseResults(Q5(2),csv,1.0);
    csv.close();

}

void parseResults(const eigenSolut& answer, std::ofstream& csv,double delta) {
    MatrixXd fission,source;
    MatrixXd::Index maxRow,maxCol;
    double maxVal,distance;
    int i;

    fission=MatrixXd((answer.flux.rows()/2),1);
    source=answer.fission;
    //dump the K_eff term right away
    csv<<answer.K<<",";
    
    for(int i=0;i<source.rows();i+=2) { //iterate over all of the fast groups
        fission(i/2,0)=source(i,0);//pull out only the fission sources
    }
    fission=fission/fission.mean(); //normalize by the average
    maxVal=fission.maxCoeff(&maxRow, &maxCol);
    csv<<maxVal<<","; //prints out the maximum fission peak
    distance=delta*(maxRow-(double)(fission.rows()/2.0)); //calculate distance from center
    
    csv<<distance<<","; //print that stuff

    csv<<answer.loopCount<<std::endl;

}
int main () {
    partA();
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
    return ret;

}

eigenSolut Q1() {
    std::vector<std::vector<std::vector<double>>> cellProp, meshProp;
    MatrixXd H, F;
    std::vector<int> pinWidth;
    double delta;

    cellProp={getMaterials(1)}; //build the cell properties
    pinWidth={60};
    delta=5.0;
    meshProp=buildPropArray(cellProp,pinWidth);
    H=generateH(meshProp,delta,true,true);
    F=generateF(meshProp,delta); 
    return solveEigenProblem(H,F,2);

}
eigenSolut Q2() {
    std::vector<std::vector<std::vector<double>>> cellProp, meshProp;
    MatrixXd H, F;
    std::vector<int> pinWidth;
    double delta;

    cellProp={getMaterials(5),getMaterials(2), getMaterials(5)};
    pinWidth={5,50,5};
    delta=5.0;
    meshProp=buildPropArray(cellProp,pinWidth);
    H=generateH(meshProp,delta,true,true);
    F=generateF(meshProp,delta);
    return solveEigenProblem(H,F,2);
}

eigenSolut Q3() {
    std::vector<std::vector<std::vector<double>>> cellProp, meshProp;
    MatrixXd H, F;
    std::vector<int> pinWidth;
    double delta;

    cellProp={getMaterials(5),getMaterials(2), getMaterials(5)};
    pinWidth={25,250,25};
    delta=1.0;
    meshProp=buildPropArray(cellProp,pinWidth);
    H=generateH(meshProp,delta,true,true);
    F=generateF(meshProp,delta);
    return solveEigenProblem(H,F,2);
}

eigenSolut Q4() {
    std::vector<std::vector<std::vector<double>>> cellProp, meshProp;
    MatrixXd H, F;
    std::vector<int> pinWidth;
    double delta;

    cellProp={getMaterials(5),getMaterials(4), getMaterials(3), getMaterials(4),
		getMaterials(5)};
    pinWidth={25,15,220,15,25};
    delta=1.0;
    meshProp=buildPropArray(cellProp,pinWidth);
    H=generateH(meshProp,delta,true,true);
    F=generateF(meshProp,delta);
    return solveEigenProblem(H,F,2);
}

eigenSolut Q5(int baffleThick) {
    std::vector<std::vector<std::vector<double>>> cellProp, meshProp;
    MatrixXd H, F;
    std::vector<int> pinWidth;
    double delta;

    cellProp={getMaterials(7),getMaterials(6), getMaterials(4),getMaterials(3),
    getMaterials(4),getMaterials(6),getMaterials(7)};
    pinWidth={23,baffleThick,15,220,15,baffleThick,23};
    delta=1.0;
    meshProp=buildPropArray(cellProp,pinWidth);
    H=generateH(meshProp,delta,true,true);
    F=generateF(meshProp,delta);
    return solveEigenProblem(H,F,2);
}

