/*
Created by Marcus Tan on 8/7/2014
Copyright 2014 University of Illinois 
Purpose: this function outputs all the possible combination of two vectors.
         If uin = {u1,u2,..,um} and vin = {v1,v2,...,vn}, then
         uout = {u1,u2,...,um, u1,u2,...,um,......, u1, u2,...,um}
         vout = {v1,v1,...,v1, v2,v2,...,v2,......, vn, vn,...,vn} 
Remarks: only works for vector types defined in Armadillo library
*/
// #include "assemble.h" DO NOT INCLUDE THIS AS THIS HAS ALREADY BEEN INCLUDED IN ASSEMBLE.H
#include <cstddef> // NULL, std::size_t
#include <stdexcept>
#include "armadillo"

template<class T>                             
void combination_two_vectors(T& uout, 
                             T& vout, 
                             const T& uin,
                             const T& vin)                             
{
    std::size_t mn = uin.n_elem*vin.n_elem;
    if (uout.n_elem != mn)
        uout.set_size(mn);
    if (vout.n_elem != mn)
        vout.set_size(mn);
    
    std::size_t k;
    for (std::size_t j = 0; j < vin.n_elem; j++)
    {
        k = j*uin.n_elem;
        for (std::size_t i = 0; i < uin.n_elem; i++)
        {
            uout(k+i) = uin(i);
            vout(k+i) = vin(j); 
        }
    }
       
}

