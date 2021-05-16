/*
Created by Marcus Tan on 8/7/2014
Update on 10/31/2014
Copyright 2014 University of Illinois 
Purpose: add Neumann BC to load vector
*/

#include "assemble.h"
#include <cstddef>
#include "armadillo"
#include <iostream>

namespace igfem
{
void impose_Neumann_BC(arma::vec& Pglo,
                       const Neumann& neumann,
                       const arma::umat& elemNodes,
                       const arma::mat& nodeCoords)
{
    if (neumann.elem.n_elem == 0)
        return;

    if (arma::max(neumann.elem) >= elemNodes.n_cols)
        std::cerr << "error: impose_Neumann_BC: a Neumann element number exceeds size of elemNodes" << std::endl;
    arma::uvec nodes(3);
    arma::uvec edgeNodes(2);
    arma::vec edgeVec(2);
    for (std::size_t i = 0; i < neumann.elem.n_elem; i++)
    {
        nodes = elemNodes.col(neumann.elem(i));
        edgeNodes = nodes(neumann.surfLocNodes(arma::span::all,neumann.surf(i)));
           
        edgeVec = nodeCoords(arma::span::all,edgeNodes(1)) 
                 -nodeCoords(arma::span::all,edgeNodes(0)),
        Pglo(edgeNodes) += 0.5*arma::norm(edgeVec,2)*neumann.vals(i);
    }

}

}

