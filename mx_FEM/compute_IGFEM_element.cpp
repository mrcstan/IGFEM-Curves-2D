/*
Created by Marcus Tan on 8/17/2014
Modified on 10/7/2015
Copyright 2014 University of Illinois 
Purpose: this function calculates the IGFEM element stiffness matrix and load vector
INPUT:
    nodeCoords:
    paLocOrigNodes: local number of original nodes of the parent element wrt to itself
    elemHeatSource:
    gauss1:
    parent1:
        parent1.children(i).locPaNodes: local number of child nodes wrt parent element
        parent1.children(i).locPaEnNodes: local number of enrichment node in child wrt parent
        parent1.children(i).locEnNodes: local number of enrichment node in child wrt itself
        parent1.children(i).conductivity
    channels:
OUTPUT:
    Kel:
    Pel:
    errFlag = 0: no negative or zero Jacobian
            < 0: -errFlag is the largest child element number with neg Jacobian
            > 0: errFlag is the largest child element number with vanishing Jacobian
           
*/
#include "assemble.h"
#include "armadillo"
#include <cstddef>
//#include <iostream>
#include <mex.h>

//#define SOURCE_FUNCTION

namespace igfem
{
void compute_IGFEM_element(arma::mat& Kel,
                           arma::vec& Pel,
                           int& errFlag,
                           const arma::mat& nodeCoords,
                           const arma::uvec& paLocOrigNodes,
                           double elemHeatSource,
                           const Convection& convect,
                           const gauss& gauss1,
                           const parent& parent1,
                           const chanNetwork& channels,
                           bool supg)
{
    Kel.zeros(parent1.nodes.n_elem, parent1.nodes.n_elem);
    Pel.zeros(parent1.nodes.n_elem);
    errFlag = 0;
    
    double factor;
    arma::vec N;
    arma::mat B;
    double detJ;
    bool calcBJch, calcBel;
    arma::mat Xpa = nodeCoords.cols(parent1.nodes);
    arma::mat Xel = Xpa.cols(paLocOrigNodes);
    arma::mat Xch;
    arma::vec Xglo(nodeCoords.n_rows);
    arma::vec locCoord(nodeCoords.n_rows);
    
    // for applying constraints
    arma::uvec cstrApplied = arma::zeros<arma::uvec>(parent1.cstrLocNodes.n_elem);
    arma::mat cstrEqns(parent1.nodes.n_elem,parent1.cstrLocNodes.n_elem);

    arma::mat Cmat = parent1.conductivity*arma::eye<arma::mat>(
                                          nodeCoords.n_rows,
                                          nodeCoords.n_rows);

    // parent element is a linear element.
    // calcBel ensures that the B matrix is only calculated once
    calcBel = true;
    arma::mat Bel;

    // for supg
    // double he;
    arma::vec Bsw;
    arma::vec Wfunc;
    
    elemHeatSource += convect.coef*convect.Tref;

    for (std::size_t ch = 0; ch < parent1.children.n_elem; ch++)
    {
        arma::mat channelVecs;
        if (parent1.children(ch).channelNum.n_elem > 0)
            channelVecs.set_size(nodeCoords.n_rows,parent1.children(ch).channelNum.n_elem);

        for (std::size_t i = 0; i < parent1.children(ch).channelNum.n_elem; i++)
        {
                channelVecs.col(i)
                        = nodeCoords.col(parent1.children(ch).channelNodes(1,i))
                         -nodeCoords.col(parent1.children(ch).channelNodes(0,i));
        }
        arma::mat unitChannelVecs;
        if (supg && parent1.children(ch).channelNum.n_elem > 0)
        {
            unitChannelVecs.set_size(channelVecs.n_rows,channelVecs.n_cols);
            double denom;
            for (std::size_t i = 0; i < channelVecs.n_cols; i++)
            {
                denom = arma::norm(channelVecs.col(i));
                if (fabs(denom) < DENOMTOL)
                {
                    //std::cerr << "compute IGFEM element: channel vec mag close to 0 " << std::endl;
                    mexWarnMsgIdAndTxt("compute_IGFEM_element:chanVec","compute_IGFEM_element:vanishing chan vec mag");
                }
                unitChannelVecs.col(i) = channelVecs.col(i)
                                        /denom;
            }
        }    
        // parent1.children(ch).channelNum.print("channel num = "); 
        // unitChannelVecs.print("unit channel vecs =  ");
        // since the child element is a linear element,
        // its B matrix and Jacobian are constants
        // calcBJch ensures that these quantities are calcuated
        // once only per child element
        calcBJch = true;
        arma::mat Bch;

        Xch = Xpa.cols(parent1.children(ch).locPaNodes);
        //Cmat = parent1.children(ch).conductivity
        //      *arma::eye<arma::mat>(nodeCoords.n_rows,
        //                            nodeCoords.n_rows);


        if (supg && parent1.children(ch).channelNum.n_elem > 0)
            if (parent1.children(ch).shape == TRIANGLE)
            {
                for (std::size_t i = 0; i < gauss1.elem.n_cols; i++)
                {
                    child_element_NBJ(N,
                                      B,
                                      detJ,
                                      Xglo, // global coordinates of gauss point
                                      calcBJch,
                                      Bch,
                                      calcBel,
                                      Bel,
                                      parent1.nodes.n_elem,
                                      paLocOrigNodes,
                                      parent1.children(ch).locPaEnNodes,
                                      parent1.children(ch).locEnNodes,
                                      Xel,
                                      Xch,
                                      gauss1.elem(arma::span(0,gauss1.elem.n_rows-2),i),
                                      parent1.children(ch).shape);
                    factor = detJ*gauss1.elem(gauss1.elem.n_rows-1,i);
                    
                    supg_weighing_function(Wfunc,
                                           unitChannelVecs,
                                           N,
                                           B,
                                           channels.mcf(parent1.children(ch).channelNum));
 
                    
                    Kel += (B.t()*Cmat*B + convect.coef*Wfunc*N.t())*factor;
                                      
                    Pel += Wfunc*elemHeatSource*factor;
                    #ifdef SOURCE_FUNCTION
                        Pel += Wfunc*body_source(Xglo)*factor;
                    #endif

                }
                if (fabs(detJ) <= JACTOL)
                    errFlag = ch;
                else if (detJ <= -JACTOL)
                    errFlag = -ch;

            }
            else if (parent1.children(ch).shape == QUADRILATERAL)
            {
                for (std::size_t i = 0; i < gauss1.quadElem.n_cols; i++)
                {
                    child_element_NBJ(N,
                                      B,
                                      detJ,
                                      Xglo, // global coordinates of gauss point
                                      calcBJch,
                                      Bch,
                                      calcBel,
                                      Bel,
                                      parent1.nodes.n_elem,
                                      paLocOrigNodes,
                                      parent1.children(ch).locPaEnNodes,
                                      parent1.children(ch).locEnNodes,
                                      Xel,
                                      Xch,
                                      gauss1.quadElem(arma::span(0,gauss1.quadElem.n_rows-2),i),
                                      parent1.children(ch).shape);
                    factor = detJ*gauss1.quadElem(gauss1.quadElem.n_rows-1,i);
                    supg_weighing_function(Wfunc,
                                           unitChannelVecs,
                                           N,
                                           B,
                                           channels.mcf(parent1.children(ch).channelNum));                    
                    Kel += (B.t()*Cmat*B + convect.coef*Wfunc*N.t())*factor;
                    Pel += Wfunc*elemHeatSource*factor;
                    #ifdef SOURCE_FUNCTION
                        Pel += Wfunc*body_source(Xglo)*factor;
                    #endif
                    if (fabs(detJ) <= JACTOL)
                        errFlag = ch;
                    else if (detJ <= -JACTOL)
                        errFlag = -ch;

                }

            }
            else
            {
                //std::cerr << "compute_IGFEM_element: unknown child element shape" << std::endl;
                mexErrMsgIdAndTxt("compute_IGFEM_element:child_shape","compute_IGFEM_element:unknown child shape");
            }
        else
        {
            if (parent1.children(ch).shape == TRIANGLE)
            {
                for (std::size_t i = 0; i < gauss1.elem.n_cols; i++)
                {
                    child_element_NBJ(N,
                                      B,
                                      detJ,
                                      Xglo, // global coordinates of gauss point
                                      calcBJch,
                                      Bch,
                                      calcBel,
                                      Bel,
                                      parent1.nodes.n_elem,
                                      paLocOrigNodes,
                                      parent1.children(ch).locPaEnNodes,
                                      parent1.children(ch).locEnNodes,
                                      Xel,
                                      Xch,
                                      gauss1.elem(arma::span(0,gauss1.elem.n_rows-2),i),
                                      parent1.children(ch).shape);
                    factor = detJ*gauss1.elem(gauss1.elem.n_rows-1,i);
                    Kel += (B.t()*Cmat*B + convect.coef*N*N.t())*factor;
                    Pel += N*elemHeatSource*factor;
                    
                    #ifdef SOURCE_FUNCTION
                        Pel += N*body_source(Xglo)*factor;
                    #endif
                }
                
                if (fabs(detJ) <= JACTOL)
                    errFlag = ch;
                else if (detJ <= -JACTOL)
                    errFlag = -ch;

            }
            else if (parent1.children(ch).shape == QUADRILATERAL)
            {
                for (std::size_t i = 0; i < gauss1.quadElem.n_cols; i++)
                {
                    child_element_NBJ(N,
                                      B,
                                      detJ,
                                      Xglo, // global coordinates of gauss point
                                      calcBJch,
                                      Bch,
                                      calcBel,
                                      Bel,
                                      parent1.nodes.n_elem,
                                      paLocOrigNodes,
                                      parent1.children(ch).locPaEnNodes,
                                      parent1.children(ch).locEnNodes,
                                      Xel,
                                      Xch,
                                      gauss1.quadElem(arma::span(0,gauss1.quadElem.n_rows-2),i),
                                      parent1.children(ch).shape);
                    factor = detJ*gauss1.quadElem(gauss1.quadElem.n_rows-1,i);
                    Kel += (B.t()*Cmat*B + convect.coef*N*N.t())*factor;
                    Pel += N*elemHeatSource*factor;
                    
                    #ifdef SOURCE_FUNCTION
                        Pel += N*body_source(Xglo)*factor;
                    #endif
                    
                    if (fabs(detJ) <= JACTOL)
                        errFlag = ch;
                    else if (detJ <= -JACTOL)
                        errFlag = -ch;
                }

            }
            else
            {
                //std::cerr << "compute_IGFEM_element: unknown child element shape" << std::endl;
                mexErrMsgIdAndTxt("compute_IGFEM_element:child_shape","compute_IGFEM_element:unknown child shape");
            }
        }

        
        
        if (channels.model == MEAN_TEMP)
        {
            int edge;
            for (std::size_t i = 0; i < parent1.children(ch).channelNum.n_elem; i++)
            {
                //channelVec = nodeCoords.col(parent1.children(ch).channelNodes(1,i))
                //            -nodeCoords.col(parent1.children(ch).channelNodes(0,i));
                edge = edge_number(parent1.children(ch).channelLocNodes.col(i),
                                   parent1.children(ch).shape);  
                if (supg)
                {
                    for (std::size_t j = 0; j < gauss1.line.n_cols; j++)
                    {
                        
                        Xglo = nodeCoords.col(parent1.children(ch).channelNodes(0,i)) 
                            + gauss1.line(0,j)*channelVecs.col(i);
                        locCoord = local_coord_2D_along_edge(Xglo, 
                                                             Xch, 
                                                             parent1.children(ch).shape,
                                                             edge);
                        child_element_NBJ(N,
                                          B,
                                          detJ,
                                          Xglo,  
                                          calcBJch,
                                          Bch,
                                          calcBel,
                                          Bel,
                                          parent1.nodes.n_elem,
                                          paLocOrigNodes,
                                          parent1.children(ch).locPaEnNodes,
                                          parent1.children(ch).locEnNodes,
                                          Xel,
                                          Xch,
                                          locCoord,
                                          parent1.children(ch).shape);
                        /*
                        Wfunc = N;
                        for (std::size_t k = 0; k < parent1.children(ch).channelNum.n_elem; k++)
                        {
                            streamwise_elem_length(he,Bsw,unitChannelVecs.col(k),B);
                            Wfunc = Wfunc + he*Bsw;
                        }
                        */

                        supg_weighing_function(Wfunc,
                                               unitChannelVecs,
                                               N,
                                               B,
                                               channels.mcf(parent1.children(ch).channelNum));
                        Kel += 0.5*channels.mcf(parent1.children(ch).channelNum(i))
                               *Wfunc*arma::trans(channelVecs.col(i))
                               *B*gauss1.line(1,j);                        
                    }

                }
                else
                {
                    for (std::size_t j = 0; j < gauss1.line.n_cols; j++)
                    {
                        Xglo = nodeCoords.col(parent1.children(ch).channelNodes(0,i)) 
                            + gauss1.line(0,j)*channelVecs.col(i);
                        locCoord = local_coord_2D_along_edge(Xglo,
                                                             Xch, 
                                                             parent1.children(ch).shape,
                                                             edge);
                        child_element_NBJ(N,
                                          B,
                                          detJ,
                                          Xglo,  
                                          calcBJch,
                                          Bch,
                                          calcBel,
                                          Bel,
                                          parent1.nodes.n_elem,
                                          paLocOrigNodes,
                                          parent1.children(ch).locPaEnNodes,
                                          parent1.children(ch).locEnNodes,
                                          Xel,
                                          Xch,
                                          locCoord,
                                          parent1.children(ch).shape);
                       // if (!arma::is_finite(B))
                       //    errFlag = ch;
                        Kel += 0.5*channels.mcf(parent1.children(ch).channelNum(i))
                               *N*arma::trans(channelVecs.col(i))
                               *B*gauss1.line(1,j);                        
                    }
                }
            }
        }
        else if (channels.model == CONST_HEAT)
        {

            int edge;  
            arma::vec N1D(2), DN1D(2);
            double channelLength;
            for (std::size_t i = 0; i < parent1.children(ch).channelNum.n_elem; i++)
            {    

                edge = edge_number(parent1.children(ch).channelLocNodes.col(i),
                                   parent1.children(ch).shape);  
                channelLength = arma::norm(channelVecs.col(i));
                for (std::size_t j = 0; j < gauss1.line.n_cols; j++)
                {
                    Xglo = nodeCoords.col(parent1.children(ch).channelNodes(0,i)) 
                        + gauss1.line(0,j)*channelVecs.col(i);
                    locCoord = local_coord_2D_along_edge(Xglo, 
                                                         Xch, 
                                                         parent1.children(ch).shape,
                                                         edge);
                    child_element_NBJ(N,
                                      B,
                                      detJ,
                                      Xglo,
                                      calcBJch,
                                      Bch,
                                      calcBel,
                                      Bel,
                                      parent1.nodes.n_elem,
                                      paLocOrigNodes,
                                      parent1.children(ch).locPaEnNodes,
                                      parent1.children(ch).locEnNodes,
                                      Xel,
                                      Xch,
                                      locCoord,
                                      parent1.children(ch).shape);
                    shape_function_1D(N1D,DN1D,gauss1.line(0,j));
                    factor = arma::dot(parent1.children(ch).channelNurbsParam.col(i),N1D)
                            * channels.lengths(parent1.children(ch).channelNum(i));
                    factor = conductance(factor,
                                         channels.kapf(parent1.children(ch).channelNum(i)),
                                         channels.mcf(parent1.children(ch).channelNum(i)),
                                         channels.eigvalsq,
                                         channels.CR1s);
                    factor *= 0.5*channelLength*gauss1.line(1,j);
                    Kel += (N*arma::trans(N))*factor;
                    Pel += N*channels.Tin*factor;
                }

            }
        } // if channels.model == MEAN_TEMP
        
        for (std::size_t i = 0; i < parent1.cstrLocNodes.n_elem; i++)
        {
            if (!cstrApplied(i) && any(parent1.children(ch).locPaNodes == parent1.cstrLocNodes(i)))
            {
                Xglo = nodeCoords.col(parent1.nodes(parent1.cstrLocNodes(i)));
                locCoord = local_coord_2D(Xglo, Xch);
                child_element_NBJ(N,
                                  B,
                                  detJ,
                                  Xglo,
                                  calcBJch,
                                  Bch,
                                  calcBel,
                                  Bel,
                                  parent1.nodes.n_elem,
                                  paLocOrigNodes,
                                  parent1.children(ch).locPaEnNodes,
                                  parent1.children(ch).locEnNodes,
                                  Xel,
                                  Xch,
                                  locCoord,
                                  parent1.children(ch).shape);
                cstrEqns.col(i) = N;
                cstrApplied(i) = 1;
            }
        }
       
  
    } // loop of child elements
   
        


   // apply Dirichlet node at channel inlet by replacing row corresponding to enrichment node by constraint equation   
   if (parent1.cstrLocNodes.n_elem)
   {
        arma::mat DN;
        arma::vec cstrEqnsPad = arma::zeros<arma::vec>(parent1.cstrLocNodes.n_elem);
        std::size_t KelNewSize = parent1.nodes.n_elem + parent1.cstrLocNodes.n_elem;
        Kel.resize(KelNewSize,KelNewSize);
        Pel.resize(KelNewSize);
        for (std::size_t i = 0; i < parent1.cstrLocNodes.n_elem; i++)
        {   
            // Don't need this part since the contribution of the original
            // shape functions is already considered when looping through
            // the child elements
            /*
            Xglo = nodeCoords.col(parent1.nodes(parent1.cstrLocNodes(i)));
            locCoord = local_coord_2D(Xglo, Xel);
            shape_function_2D(N,DN,locCoord);
            cstrEqns(arma::span(0,N.n_elem-1),i) = N;
            */
            Kel.col(parent1.nodes.n_elem+i) = arma::join_cols(cstrEqns.col(i),
                                                              cstrEqnsPad);
            Pel(parent1.nodes.n_elem+i) = parent1.cstrVals(i);             
            Kel.row(parent1.nodes.n_elem+i) = arma::trans(Kel.col(parent1.nodes.n_elem+i));
        }
    }
    
}   


}
