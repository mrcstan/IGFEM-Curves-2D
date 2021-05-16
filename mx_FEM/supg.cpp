/*
Created by Marcus Tan on 1/11/2014
Updated on 10/7/2015
Copyright 2014 University of Illinois 
Purpose: this function calculates the body source at a give position
INPUT:
    B: [dN1/dx,dN2/dx,...dNn/dx;
        dN1/dy,dN2/dy,...dNn/dy];
        n is the number of shape functions
    he: streamwise length = HALF that presented in papers
    Bsw1: unit tangent vector of channel dot with B matrix
*/
#include "assemble.h"
#include "armadillo"
#include <math.h>
#include <mex.h>
//#define HUGHES_SUPG

namespace igfem
{
void streamwise_elem_length(double& he, 
                            arma::vec& Bsw1, 
                            const arma::vec& channelVec,
                            const arma::mat& B)
{
    Bsw1 = B.t()*channelVec;
    double denom = arma::sum(arma::abs(Bsw1));
    if (fabs(denom) < DENOMTOL) 
    {
        //std::cerr << "streamwise_elem_length: denom close to zero" << std::endl;
        mexErrMsgIdAndTxt("streamwise_elem_length:denom_zero",
                           "streamwise_elem_length denom close to zero");
        
        he = 0.0;
        return; 
    }
    he = 1.0/denom; //Important: this is half the streamwise length presented in papers 
    
    //he = 2.0/arma::sum(arma::abs(Bsw1));
    return;
}

void supg_weighing_function(arma::vec& W,
                            const arma::mat& channelUnitVecs,
                            const arma::vec& N,
                            const arma::mat& B,
                            const arma::vec& mcf)
{
    arma::vec Bsw;
    W = N;
    double he;

    #ifdef HUGHES_SUPG 
        //Brooks-Hughes SUPG
        double paramsupg = 0.0;
        double summdotsq = 0.0; 
        for (std::size_t i = 0; i < channelUnitVecs.n_cols;i++)
        {
            if (fabs(mcf(i)) < MCFTOL)
                continue;
            streamwise_elem_length(he,Bsw,sgn(mcf(i))*channelUnitVecs.col(i),B);
            paramsupg += he*mcf(i);
            summdotsq += mcf(i)*mcf(i);
        }
        paramsupg /= summdotsq;
        for (std::size_t i = 0; i < channelUnitVecs.n_cols;i++)
        {
            if (fabs(mcf(i)) < MCFTOL)
                continue;
            streamwise_elem_length(he,Bsw,sgn(mcf(i))*channelUnitVecs.col(i),B);
            W += paramsupg*mcf(i)*Bsw;
        }
        // End of Brooks-Hughes SUPG
    #else
        // My SUPG
        for (std::size_t i = 0; i < channelUnitVecs.n_cols;i++)
        {
            if (fabs(mcf(i)) < MCFTOL)
                continue;
            streamwise_elem_length(he,Bsw,sgn(mcf(i))*channelUnitVecs.col(i),B);
            W += he*Bsw;
        }
    #endif

    return;  
}


}
