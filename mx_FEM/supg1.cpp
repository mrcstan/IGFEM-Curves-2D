/*
Created by Marcus Tan on 1/11/2014
Copyright 2014 University of Illinois 
Purpose: this function calculates the body source at a give position
INPUT:
    B: [dN1/dx,dN2/dx,...dNn/dx;
        dN1/dy,dN2/dy,...dNn/dy];
        n is the number of shape functions
*/
#include "assemble.h"
#include "armadillo"
#include <math.h>
#include <iostream>

namespace igfem
{
void streamwise_elem_length(double& he, 
                            arma::vec&Bsw1, 
                            const arma::vec& channelVec,
                            const arma::mat& B)
{
    Bsw1 = B.t()*channelVec;
    double denom = arma::sum(arma::abs(Bsw1));
    if (fabs(denom) < DENOMTOL)
        std::cerr << "streamwise_elem_length: denom close to zero" << std::endl;
    he = 2.0/denom;
    
    //he = 2.0/arma::sum(arma::abs(Bsw1));
    return;
}

}
