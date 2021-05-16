/*
Created by Marcus Tan on 8/12/2014
Copyright 2014 University of Illinois
Purpose: Partition the global stiffness matrix and load vector into 
         free and prescribed blocks according to ii, jj and eqNum.
         Kff, Kfp, Kpf and Kpp are given in sparse matrix format
Input:
    ii:
    jj:
    Kval:
    Pglo:
    eqNum:
Output:
    Kff:
    Kfp:
    Kpf:
    Kpp:
    Pf:
    Pp:
         
*/
#include "assemble.h"
#include <cstddef> // NULL, Std::size_t
#include "armadillo"
#include <iostream>
#include "mex.h"

namespace igfem
{
void partition_stiffness_mat_n_load_vec(arma::umat& ffInd,
                                        arma::vec& ffKval,
                                        arma::umat& fpInd,
                                        arma::vec& fpKval,
                                        arma::umat& pfInd,
                                        arma::vec& pfKval,
                                        arma::umat& ppInd,
                                        arma::vec& ppKval,
                                        arma::vec& Pf,
                                        arma::vec& Pp,
                                        const arma::ivec& ii,
                                        const arma::ivec& jj,
                                        const arma::vec& Kval,
                                        const arma::vec& Pglo,
                                        const arma::ivec& eqNum,
                                        std::size_t estNff,
                                        std::size_t estNfp,
                                        std::size_t estNpf,
                                        std::size_t estNpp)
{
    /* buggy
    arma::uvec fInd = arma::find(eqNum > 0);
    arma::uvec pInd = arma::find(eqNum < 0);
    Pf = Pglo(fInd);   
    Pp = Pglo(pInd);
    fInd.clear();
    pInd.clear();
    */
    
    Pf.set_size(Pglo.n_elem);
    Pp.set_size(Pglo.n_elem);
    std::size_t nf = 0, np = 0;
    for (std::size_t i = 0; i < eqNum.n_elem; i++)
    {
        if (eqNum(i) > 0)
        {
            Pf(eqNum(i)-1) = Pglo(i);
            nf++;
        }
        else if (eqNum(i) < 0)
        {
            Pp(-eqNum(i)-1) = Pglo(i);
            np++;
        }
        else
        {
            //std::cerr << "error: partition_stiffness_mat_n_load_vec: eqNum cannot be zero" << std::endl;
            mexErrMsgIdAndTxt("partition_stiffness_mat_n_load_vec:eqNum_zero",
                              "partition_stiffness_mat_n_load_vec eqNum cannot be zero");
        }
    }
    Pf.resize(nf);
    Pp.resize(np);

    ffInd.set_size(2,estNff);
    ffKval.set_size(estNff);
    fpInd.set_size(2,estNfp);
    fpKval.set_size(estNfp);
    pfInd.set_size(2,estNpf);
    pfKval.set_size(estNpf);
    ppInd.set_size(2,estNpp);
    ppKval.set_size(estNpp);
    
    std::size_t nffInd = 0, nfpInd = 0, npfInd = 0, nppInd = 0;
    #pragma omp parallel for   
    for (int i = 0; i < (int) ii.n_elem; i++)
    {
        if (ii(i) > 0 && jj(i) > 0)
        {
            #pragma omp critical (FF)
            {
                if (nffInd >= estNff)
                {
                    estNff *= 2;
                    ffInd.resize(2,estNff);
                    ffKval.resize(estNff);
                    //std::cerr << "partition: warning: ffInd resized" << std::endl;
                    mexWarnMsgIdAndTxt("partition_stiffness_mat_n_load_vec:ffInd_resized",
                                      "partition_stiffness_mat_n_load_vec ffInd resized");

                }
                ffInd(0,nffInd) = ii(i)-1;
                ffInd(1,nffInd) = jj(i)-1;
                ffKval(nffInd++) = Kval(i);
            }
        } 
        else if (ii(i) > 0 && jj(i) < 0)
        {
            #pragma omp critical (FP)
            {
                if (nfpInd >= estNfp)
                {
                    estNfp *= 2;
                    fpInd.resize(2,estNfp);
                    fpKval.resize(estNfp);
                    //std::cerr << "partition: warning: fpInd resized" << std::endl;
                    mexWarnMsgIdAndTxt("partition_stiffness_mat_n_load_vec:fpInd_resized",
                                       "partition_stiffness_mat_n_load_vec fpInd resized");
                }
                fpInd(0,nfpInd) = ii(i)-1;
                fpInd(1,nfpInd) = -jj(i)-1;
                fpKval(nfpInd++) = Kval(i);
            }

        } 
        else if (ii(i) < 0 && jj(i) > 0)
        {
            #pragma omp critical (PF)
            {
                
                if (npfInd >= estNpf)
                {
                    estNpf *= 2;                   
                    pfInd.resize(2,estNpf);
                    pfKval.resize(estNpf);
                    //std::cerr << "partition: warning: pfInd resized" << std::endl;
                    mexWarnMsgIdAndTxt("partition_stiffness_mat_n_load_vec:pfInd_resized",
                                      "partition_stiffness_mat_n_load_vec pfInd resized");
                }
                pfInd(0,npfInd) = -ii(i)-1;
                pfInd(1,npfInd) = jj(i)-1;
                pfKval(npfInd++) = Kval(i);
            }
        }
        else if (ii(i) < 0 && jj(i) < 0)
        {
            #pragma omp critical (PP)
            {
                if (nppInd >= estNpp)
                {
                    estNpp *= 2;
                    ppInd.resize(2,estNpp);
                    ppKval.resize(estNpp);
                    //std::cerr << "partition: warning: ppInd resized" << std::endl;
                    mexWarnMsgIdAndTxt("partition_stiffness_mat_n_load_vec:ppInd_resized",
                                      "partition_stiffness_mat_n_load_vec ppInd resized");
                }

                ppInd(0,nppInd) = -ii(i)-1;
                ppInd(1,nppInd) = -jj(i)-1;
                ppKval(nppInd++) = Kval(i);
            }
        }
    }
    ffInd.resize(2,nffInd);
    ffKval.resize(nffInd);
    fpInd.resize(2,nfpInd);
    fpKval.resize(nfpInd);
    pfInd.resize(2,npfInd);
    pfKval.resize(npfInd);
    ppInd.resize(2,nppInd);
    ppKval.resize(nppInd);
   
   
}
}

