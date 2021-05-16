/*
Created by Marcus Tan on 8/7/2014
Update on 1/21/2015
Copyright 2014 University of Illinois 
Purpose: this function assembles the local stiffness matrix of each element 
         and output the results KFF, KFP, KPF and KPP in sparse matrix format
INPUT: 
    elemNodes:
OUTPUT:
    ii:
    jj:
    Kval:
    Pglo: global load vector
    estSize: estimated size of ii, jj and Kval
    nodeCoords: a number of dimension x number of nodes matrix of nodal coodinates
    elemNodes: a number of nodes per element x number of elements matrix of the element nodal connectivity
    conductivity: a vector of length number of elements of the element conductivities
    elemHeatSource: a vector of the element heat sources
    eqNum: a vector of length number of nodes of the equation number corresponding to each node
    gauss1: see header file "assemble.h" for def
    parents: see header file "assemble.h" for def
    channels: see header file "assemble.h" for def

*/

#include "assemble.h"
#include <omp.h>
#include <cstddef> // NULL, std::size_t
//#include <iostream>
#include <mex.h>
#include "armadillo"

namespace igfem
{
void assemble(arma::ivec& ii,
              arma::ivec& jj, 
              arma::vec& Kval,
              arma::vec& Pglo,
              std::size_t estSize, 
              const arma::mat& nodeCoords,
              const arma::umat& elemNodes,
              const arma::vec& elemHeatSource,
              const Convection& convect,
              const arma::ivec& eqNum,
              const gauss& gauss1,
              const arma::field<parent>& parents,
              const chanNetwork& channels,
              const Neumann& neumann,
              bool supg)
{
    // number of equations
    //std::size_t nEqs = eqNum.n_elem;
    
    // initialize Pglo to 0
    Pglo.zeros(eqNum.n_elem);
   
    
    // set sizes of ii, jj and Kval
    // std::size_t estSize = 2*elemNodes.n_cols*elemNodes.n_rows*elemNodes.n_rows;
    ii.set_size(estSize);
    jj.set_size(estSize);
    Kval.set_size(estSize);
    
    int nThreads = omp_get_max_threads();
    //std::cout << "number of threads = " << nThreads << std::endl;
    // locVecSize[tid] is the size of the vector locii of thread tid
    arma::ivec locVecSize = arma::zeros<arma::ivec>(nThreads);
    // locCumSUm is the cummulative sum of {0, locVecSize}
    arma::ivec locCumSum = arma::zeros<arma::ivec>(nThreads+1); 
    

    // assume that the parent element nodes are arranged such that the first nodes
    // are the original nodes
    arma::uvec paLocOrigNodes(elemNodes.n_rows);
    for (std::size_t i = 0; i < elemNodes.n_rows; i++)
        paLocOrigNodes(i) = i;
   
   
    omp_lock_t writelock;
    omp_init_lock(&writelock);
    omp_set_num_threads(1); // REMOVE THIS LINE TO USE MULTIPLE THREADS
    #pragma omp parallel  
    {

        int threadNum = omp_get_thread_num();
        int numThreads = omp_get_num_threads();
        std::size_t newVecSize = 0, oldVecSize = 0;
        int errFlag;

        arma::ivec locii = arma::zeros<arma::ivec>(estSize);
        arma::ivec locjj = arma::zeros<arma::ivec>(estSize);
        arma::vec locKval = arma::zeros<arma::vec>(estSize);
        //arma::vec locPglo = arma::zeros<arma::vec>(nEqs);
        arma::mat Kel(elemNodes.n_rows,elemNodes.n_rows);
        arma::vec Pel(elemNodes.n_rows);
        arma::uvec elii(elemNodes.n_rows*elemNodes.n_rows);
        arma::uvec eljj(elemNodes.n_rows*elemNodes.n_rows);
        arma::uvec elemNodesCol(elemNodes.n_rows);
        arma::uvec parentNodesnCstrs; 
        // calculation of element stiffness matrix and load vector is distributed to multiple threads
        int istart = threadNum*elemNodes.n_cols/numThreads;
        int iend;
        if (threadNum == numThreads-1)
            iend = (int) elemNodes.n_cols;
        else
            iend = istart + elemNodes.n_cols/numThreads;
       
        /* 
        omp_set_lock(&writelock);
        std::cout << "thread = " << threadNum 
                  << ", istart = " << istart 
                  << ", iend = " << iend << std::endl;
        std::cout << "running in parallel = " << omp_in_parallel() << std::endl;
        omp_unset_lock(&writelock);
        */

        for(int i = istart; i < iend; ++i)   
        {
            if (parents(i).type == REGULAR)
            {
                //std::cout << "thread " << threadNum << ", regular element " << i << std::endl;
                compute_regular_element(Kel,  
                                        Pel,
                                        nodeCoords, 
                                        elemNodes.col(i), 
                                        elemHeatSource(i),
                                        convect,
                                        gauss1,
                                        parents(i),
                                        channels,
                                        supg);
                
                // update thread copies of ii, jj and Kval
                elemNodesCol = elemNodes.col(i);           
                combination_two_vectors(elii, eljj,elemNodesCol ,elemNodesCol);
                #pragma omp critical
                {
                    Pglo(elemNodesCol) += Pel;
                }
                //locPglo(elemNodesCol) += Pel;
            }
            else if (parents(i).type == IGFEM)
            {
                //std::cout << "thread " << threadNum << ", IGFEM element " << i << std::endl;
                compute_IGFEM_element(Kel,
                                      Pel,
                                      errFlag,
                                      nodeCoords,
                                      paLocOrigNodes,
                                      elemHeatSource(i),
                                      convect,
                                      gauss1,
                                      parents(i),
                                      channels,
                                      supg);
                /*   
                omp_set_lock(&writelock);
                std::cout << "thread " << threadNum << ", IGFEM element " << i << std::endl;
                Kel.print("Kel = ");
                omp_unset_lock(&writelock);
                */
                if (errFlag > 0)
                {
                    omp_set_lock(&writelock);
                    {
                        mexWarnMsgIdAndTxt("assemble:IGFEM_zero_Jacobian",
                                        "compute_IGFEM_element:vanishing Jacobian in parent %i, child %i\n",i,errFlag);
                    }
                    omp_unset_lock(&writelock);
                } else if (errFlag < 0)
                {
                    omp_set_lock(&writelock);
                    {
                        mexWarnMsgIdAndTxt("assemble:IGFEM_negative_Jacobian",
                                        "compute_IGFEM_element:negative Jacobian in parent %i, child %i\n",i,-errFlag);
                    }
                    omp_unset_lock(&writelock);
                }

                /*
                else if (errFlag > 0)
                {
                    std::cerr << "warning: compute_IGFEM_element: inf or nan in B matrix of parent " 
                    << i << ", child "
                    << errFlag << std::endl;
                }
                */
                parentNodesnCstrs = arma::join_cols(arma::conv_to<arma::umat>::from(parents(i).nodes),
                                               arma::conv_to<arma::umat>::from(parents(i).cstrRows));
                //combination_two_vectors(elii, eljj, parents(i).nodes, parents(i).nodes);
                combination_two_vectors(elii, eljj, parentNodesnCstrs, parentNodesnCstrs);

                #pragma omp critical
                {
                    Pglo(parentNodesnCstrs) += Pel;
                }
                //locPglo(parentNodesnCstrs) += Pel;
              
            }
            newVecSize = oldVecSize + elii.n_elem;
            if (newVecSize > locii.n_elem)
            {
                mexWarnMsgIdAndTxt("assemble:growing_locii",
                                   "assemble:growing locii",i,-errFlag);
                locii.resize(2*newVecSize);
                locjj.resize(2*newVecSize);
                locKval.resize(2*newVecSize);
            }
            locii(arma::span(oldVecSize,newVecSize-1)) =  eqNum(elii);
            locjj(arma::span(oldVecSize,newVecSize-1)) =  eqNum(eljj);
            locKval(arma::span(oldVecSize,newVecSize-1)) = arma::vectorise(Kel);
            oldVecSize = newVecSize;
        }    
        
        locVecSize(threadNum) = newVecSize;                    
        
        // debug !!!
        /*
        locii.resize(newVecSize);
        locKval.resize(newVecSize);
        omp_set_lock(&writelock);
        locii.print("locii = ");
        locKval.print("locKval = ");
        omp_unset_lock(&writelock);
        */
        #pragma omp barrier // wait until all threads have update locVecSize, lociiPtr, locjjPtr and locKvalPtr 
        
        //#pragma omp flush (locVecSize,locCumSum,lociiPtr,locjjPtr,locKvalPtr)
        
        #pragma omp master
        {
            locVecSize.print("number of computations of ii entries by each thread = ");
            locCumSum(arma::span(1,locCumSum.n_elem-1)) = locVecSize;
	        locCumSum = arma::cumsum(locCumSum);
            if (locCumSum(numThreads) != (int)ii.n_elem)
            {
                ii.set_size(locCumSum(numThreads));
                jj.set_size(locCumSum(numThreads));
                Kval.set_size(locCumSum(numThreads));
            }

        }
                
        //#pragma omp flush (locCumSum)

        #pragma omp barrier
        /*
        omp_set_lock(&writelock);
        std::cout << "thread = " << threadNum << std::endl;
        locCumSum.print("locCumSum = ");
        locVecSize.print("locVecSize = ");
        std::cout << "locii_size = " << locii.n_elem 
                  << ", locjj_size = " << locjj.n_elem
                  << ", locKval_size = " << locKval.n_elem << std::endl;
        omp_unset_lock(&writelock);
        */
        // aggregate the private copies of vector
        #pragma omp critical (triplets)
        for(int j = 0; j < locVecSize(threadNum); j++)         
        {
            ii(locCumSum(threadNum)+j) =  locii(j);
            jj(locCumSum(threadNum)+j) =  locjj(j);
            Kval(locCumSum(threadNum)+j) = locKval(j);
        
        }
        /*
        #pragma omp critical (globalForce)
        for (int j = 0; j < nEqs; j++)
            Pglo(j) += locPglo(j);
       
        */
        
        
       
    }

    impose_Neumann_BC(Pglo, neumann, elemNodes, nodeCoords);
}
}
