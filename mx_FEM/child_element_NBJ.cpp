/*
Created by Marcus Tan on 8/17/2014
Modified on 10/26/2014
Copyright 2014 University of Illinois 
Purpose: this function calculates shape functions, B matrix of the child 
        and parent elements as well as the Jacobian of the child element
        mapping to the global space
INPUT:
OUTPUT: 
    N: shape functions
    B: derivative of shape functions wrt global coordinates 
		(note: size of matrix is ndim x nNodesPerElem. 
				in comparison, the size of DN is nNodes x nNodesPerElem)
			
    detJ: determinant of Jacobian
    X: global coordinates corresponding to the local coordinates
INPUT_OUTPUT: 
    calcBJch: flag indicating whether the B matrix of the child element should be calculated. if true, it calculates the
    matrix and returns a false
    calcBel: flag indicating whether the B matrix of the parent element should be calcul             ated. if true, it
    calculates the matrix and returns a false
*/
#include "assemble.h"
#include "armadillo"
#include <cstddef>

namespace igfem
{
void child_element_NBJ(arma::vec& N,
                       arma::mat& B,
                       double& detJ,
                       arma::vec& X,
                       bool& calcBJch,
                       arma::mat& Bch,
                       bool& calcBel,
                       arma::mat& Bel,
                       std::size_t nNodes,
                       const arma::uvec& paLocOrigNodes, // parent local number of original nodes
                       const arma::uvec& chLocPaEnrichNodes, // child local number of enrichnment nodes wrt parent
                       const arma::uvec& chLocEnrichNodes, // child local number of enrichment nodes wrt itself
                       const arma::mat& Xel,
                       const arma::mat& Xch,
                       const arma::vec& locCoord,
                       shapes shape)
{
   
    N = arma::zeros<arma::vec>(nNodes);
    B = arma::zeros<arma::mat>(Xel.n_rows,nNodes);

    // calculate shape functions, B matrix and Jacobian of child element
    arma::vec Nch; // child shape function
    arma::mat DNch; 
    shape_function_2D(Nch,DNch,locCoord,shape);
    if (calcBJch)
    {
        //std::cout.precision(16);
        //std::cout.setf(std::ios::scientific);
        //Xch.print("Xch = ");
        arma::mat Jch = Xch*DNch;
        //Jch.raw_print(std::cout,"Jch = ");
        detJ = arma::det(Jch);
        //std::cout << "detJ = " << detJ << std::endl;
        if (detJ > JACTOL)
            Bch = arma::solve(Jch.t(), DNch.t());
        else
            Bch = arma::zeros<arma::mat>(DNch.n_cols,DNch.n_rows);
       if (shape == TRIANGLE)   
            calcBJch = false;
    }
    N(chLocPaEnrichNodes) = Nch(chLocEnrichNodes);
    //B(rowInds,chLocPaEnrichNodes) = Bch(rowInds,chLocEnrichNodes);
    B.cols(chLocPaEnrichNodes) = Bch.cols(chLocEnrichNodes);
    
    // calculate shape functions, B matrix of nodes belonging to parent element
    X = Xch*Nch; // global coordinates
    arma::vec locCoordPa= local_coord_2D(X,Xel);
  
    arma::vec Nel;
    arma::mat DNel;
    shape_function_2D(Nel,DNel,locCoordPa);
    
    if (calcBel)
    {
        arma::mat Jel = Xel*DNel;
        Bel = arma::solve(Jel.t(),DNel.t());
        calcBel = false;
    }
   
    N(paLocOrigNodes) = Nel;
    //B(rowInds,paLocOrigNodes) = Bel;
    B.cols(paLocOrigNodes) = Bel;
    
}
}
