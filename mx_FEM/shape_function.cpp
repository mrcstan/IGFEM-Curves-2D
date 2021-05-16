/*
Created by Marcus Tan on 8/11/2014
Copyright 2014 University of Illinois 
Purpose: this function calculates the shape functions of a triangular 
         or tetrahedral element
INPUT: 
    locCoords: a length number of nodes column vector of local coordinate
OUTPUT:
    N: a length number of nodes column vector of shape functions
        N = [1-sum(locCoords), locCoords(0), ...locCoords(end))]
    DN: a number of nodes x number of dimensions column vector
        DN = [-1,-1,-1,...;1,0,0...;0,1,0,...;0,0,1,...]
*/
#include "assemble.h"
#include "armadillo"
#include <iostream>

namespace igfem
{
// linear shape functions for tetrahedral elements
void shape_function_3D(arma::vec& N, arma::mat& DN, const arma::vec& locCoord)
{
    if (locCoord.n_elem != 3)
        std::cerr << "shape_function_3D: dimension of locCoord must be 3" << std::endl;
    N.set_size(4);
    DN.set_size(4, 3);
        
    N(0) = 1.0-arma::sum(locCoord);
    N(arma::span(1,3)) = locCoord;
    
    DN.row(0) = -arma::ones<arma::rowvec>(3);
    DN(arma::span(1,3),arma::span(0,2)) 
        = arma::eye<arma::mat>(3,3);
}

// linear shape functions for triangular elements
void shape_function_2D(arma::vec& N, arma::mat& DN, const arma::vec& locCoord)
{
    if (locCoord.n_elem != 2)
        std::cerr << "shape_function_2D: dimension of locCoord must be 2" << std::endl;
    N.set_size(3);
    DN.set_size(3, 2);
        
    N(0) = 1.0-arma::sum(locCoord);
    N(arma::span(1,2)) = locCoord;
    
    DN.row(0) = -arma::ones<arma::rowvec>(2);
    DN(arma::span(1,2),arma::span(0,1)) 
        = arma::eye<arma::mat>(2,2);
}

void shape_function_2D(arma::vec& N, arma::mat& DN, const arma::vec& locCoord, shapes shape)
{
    if (locCoord.n_elem != 2)
        std::cerr << "shape_function_2D: dimension of locCoord must be 2" << std::endl;
    if (shape == TRIANGLE)
    {
        N.set_size(3);
        DN.set_size(3, 2);
        
        N(0) = 1.0-arma::sum(locCoord);
        N(arma::span(1,2)) = locCoord;
    
        DN.row(0) = -arma::ones<arma::rowvec>(2);
        DN(arma::span(1,2),arma::span(0,1)) 
            = arma::eye<arma::mat>(2,2);
    }
    else if (shape == QUADRILATERAL)
    {
        N.set_size(4);
        DN.set_size(4,2);
        N(0) = 0.25 * (1.0 - locCoord(0)) * (1.0 - locCoord(1));
        N(1) = 0.25 * (1.0 + locCoord(0)) * (1.0 - locCoord(1));
        N(2) = 0.25 * (1.0 + locCoord(0)) * (1.0 + locCoord(1));
        N(3) = 0.25 * (1.0 - locCoord(0)) * (1.0 + locCoord(1));
        DN(0,0) = -0.25 * (1.0 - locCoord(1));
        DN(1,0) =  0.25 * (1.0 - locCoord(1));        
        DN(2,0) =  0.25 * (1.0 + locCoord(1));
        DN(3,0) = -0.25 * (1.0 + locCoord(1));
        DN(0,1) = -0.25 * (1.0 - locCoord(0));
        DN(1,1) = -0.25 * (1.0 + locCoord(0));
        DN(2,1) =  0.25 * (1.0 + locCoord(0));
        DN(3,1) =  0.25 * (1.0 - locCoord(0));

    }
    else
        std::cerr << "shape_function_2D: unknown shape" << std::endl;

}

// linear 1D shape functions
void shape_function_1D(arma::vec& N, arma::mat& DN, double locCoord)
{
    
    N.set_size(2);
    DN.set_size(2);
    
    N(0) = 0.5 * (1.0 - locCoord);
    N(1) = 0.5 * (1.0 + locCoord);
    DN(0) = -0.5;
    DN(1) = 0.5;
}
}
