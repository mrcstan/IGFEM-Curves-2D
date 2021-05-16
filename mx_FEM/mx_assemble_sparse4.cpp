/*
Created by Marcus Tan on 8/7/2014
Modified on 12/23/2014
Copyright 2014 University of Illinois 
Purpose: this function assembles the local stiffness matrix of each element 
         and output the results KFF, KFP, KPF and KPP in sparse matrix format
Compile:
    setenv OMP_NUM_THREADS numThreads // must set this before compiling so that 
                                         mex function uses OMP_NUM_THREADS of threads.
                                         otherwise, it uses the max number of threads
    mex -v -larmadillo -lgfortran mx_assemble_sparse.cpp assemble.h assemble.cpp ...
    compute_regular_element.cpp shape_function.cpp combination_two_vectors.tcc etc etc...
    CXXFLAGS="\$CXXFLAGS -fopenmp  " LDFLAGS="\$LDFLAGS -fopenmp " % 

Input:
    NOTE: all numeric inputs are of type double
    (1) nodeCoords: a number of dimension x number of nodes matrix of nodal coodinates
    (2) elemNodesIn: a number of nodes per element x number of elements matrix of the element nodal connectivity
    (3) elemHeatSource: a vector of the element heat sources
    (4) eqNumIn: a vector of length number of nodes of the equation number corresponding to each node
    (5) gauss: a structural variable containing the gauss points and weights for 
               for the integration.
               The variable has the following fields:
        gauss.elem: a (number of dimension+1) x number of gauss points matrix,
                    the last row being the weights associated with each gauss point and 
                    the remainding rows are the gauss point local coordinates. 
                    This field is for integration over an element
        gauss.line: a 2 x number of gauss points matrix,
                    the first row is the gauss point local coordinates and
                    the last row is the weights
     (6) parent: 
         parent.channelNodes
         parent.channelNum
         parent.nShared
     (7) channel:
         channel.mcf
     (8) Dirichlet:
         Dirichlet.temp_node:
         Dirichlet.temp_val
     (9) Neumann:
          Neumann.heatFlux_elem:
          Neuamnn.heatFlux_surface:
          Neumann.heatFlux_value:
     (10) supg: apply supg, true or false
Output:
    (1): Kff
    (2): Kfp
    (3): Kpf
    (4): Kpp
    (5): Pf
    (6): Pp
     
*/
#include <iostream>
#include <cstddef> // for size_t
#include <math.h>
#include <algorithm>
#include "assemble.h"
#include "armaMex.hpp"

using namespace igfem;

enum chFields {ISTRIANGLE,LOCPANODES,LOCPAENNODES,LOCENNODES,CHANNELNUM,
               CHANNELNODES, CHANNEL_LOC_NODES,CHANNELNURBSPARAM,CHLAST};
const int nChFields = CHLAST - LOCPANODES;

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
    /* check for proper number of arguments */
    if(nrhs < 10) 
        mexErrMsgIdAndTxt("mx_assemble_sparse:nrhs","10 inputs required.");
    
    if(nlhs > 6) 
        mexErrMsgIdAndTxt("mx_assemble_sparse:nlhs","Too many outputs requested.");
    
    if (!mxIsDouble(prhs[0]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input1","Input 1, nodeCoords must be double");
    
    if (!mxIsDouble(prhs[1]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input2","Input 2, elemNodes  must be double");
    
    
    if (!mxIsDouble(prhs[2]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input3","Input 3, elemHeatSource must be double");
    
    if (!mxIsDouble(prhs[3]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input4","Input 4, eqNum must be double");
    
    if (!mxIsStruct(prhs[4]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input5","Input 5, gauss must be structural");
    
    if (!mxIsStruct(prhs[5]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input6","Input 6, parent must be structural");
 
    if (!mxIsStruct(prhs[6]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input7","Input 7, channel must be structural");
    
    if (!mxIsStruct(prhs[7]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input8","Input 8, Dirichlet must be structural");
    
    if (!mxIsStruct(prhs[8]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input9","Input 9, Neumann must be structural");
    
    if (!mxIsLogicalScalar(prhs[9]))
        mexErrMsgIdAndTxt("mx_assemble_sparse:input10","Input 10, supg must be scalar boolean");
    
    bool supg = mxIsLogicalScalarTrue(prhs[9]); 
    
    if (mxGetM(prhs[0]) > 2)
        mexErrMsgIdAndTxt("mx_assemble_sparse:input0","Input 0, each col of nodeCoords correspond to a node");
    arma::mat nodeCoords = armaGetPr(prhs[0]);

    // indices start with 0 in C
    if (mxGetM(prhs[1]) > 3)
        mexErrMsgIdAndTxt("mx_assemble_sparse:input1","Input 1, each col of elemNodes correspond to an element");
    arma::umat elemNodes = arma::conv_to<arma::umat>::from(armaGetPr(prhs[1])-1);     

    arma::vec elemHeatSource;
    if (mxGetM(prhs[2]) == 1)
        elemHeatSource = arma::trans(armaGetPr(prhs[2]));
    else if(mxGetN(prhs[2]) == 1)
        elemHeatSource = armaGetPr(prhs[2]);
    else      
        mexErrMsgIdAndTxt("mx_assemble_sparse:input3","Input 3, elemHeatSource must be row or column vector");

    arma::ivec eqNum;
    if (mxGetM(prhs[3]) == 1)
        eqNum = arma::conv_to<arma::ivec>::from(arma::trans(armaGetPr(prhs[3])));  
    else if (mxGetN(prhs[3]) == 1)
        eqNum = arma::conv_to<arma::ivec>::from(armaGetPr(prhs[3]));
    else 
        mexErrMsgIdAndTxt("mx_assemble_sparse:input4","Input 4, eqNum must be row or column vector");

    // the number of elements of eqNum must be the same as the number of nodes
    if (nodeCoords.n_cols > eqNum.n_elem)
         mexErrMsgIdAndTxt("mx_assemble_sparse:eqNum_n_elem_no_of_nodes",
                           "eqNum.n_elem must be greater than or equal the number of nodes");
    
    mxArray* fieldPtr;
    int fieldNum;
     
    // validating and getting gauss integration points
    gauss gauss1;
    mexPrintf("validating and getting gauss integration points \n");
    // element integration
    fieldNum = mxGetFieldNumber(prhs[4],"elem");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_gauss_elem",
                           "gauss points for element integration not found");
    fieldPtr = mxGetFieldByNumber(prhs[4],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:gauss_elem_empty",
                            "gauss points for element integration must be > 0");
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:gauss_elem_not_double",
                            "gauss points for element integration must be of type double");
    gauss1.elem = armaGetPr(fieldPtr);  
    if(gauss1.elem.n_rows != nodeCoords.n_rows+1)
          mexErrMsgIdAndTxt("mx_assemble_sparse:gauss_elem_wrong_dim",
                            "wrong dimension for element integration");
    
    // quadrilateral element integration
    fieldNum = mxGetFieldNumber(prhs[4],"quadElem");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_gauss_quad_elem",
                           "gauss points for quad element integration not found");
    fieldPtr = mxGetFieldByNumber(prhs[4],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:gauss_quad_elem_empty",
                            "gauss points for quad element integration must be > 0");
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:gauss_quad_elem_not_double",
                            "gauss points for quad element integration must be of type double");
    gauss1.quadElem = armaGetPr(fieldPtr);  
    if(gauss1.quadElem.n_rows != nodeCoords.n_rows+1)
          mexErrMsgIdAndTxt("mx_assemble_sparse:gauss_quad_elem_wrong_dim",
                            "wrong dimension for quad element integration");
    
    // line integration
    fieldNum = mxGetFieldNumber(prhs[4],"line");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_gauss_line",
                           "gauss points for line integration not found");
    fieldPtr = mxGetFieldByNumber(prhs[4],0,fieldNum);
    if (!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:gauss_line_empty",
                            "gauss points for line integration must be > 0");
    
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:gauss_line_not_double",
                            "gauss points for line integration must be of type double");
   
    gauss1.line = armaGetPr(fieldPtr);
    if (gauss1.line.n_rows != 2)
          mexErrMsgIdAndTxt("mx_assemble_sparse:gauss_line_wrong_dim",
                            "wrong dimension for line integration");
    
    
    
    // validating and getting parent data
    mexPrintf("validating and getting parent data \n");
    if (mxGetNumberOfElements(prhs[5]) != elemNodes.n_cols)
        mexErrMsgIdAndTxt("mx_assemble_sparse:wrong_parent_length",
                          "parent structural array must be of same length as number of elements");
    arma::field<parent> parents(elemNodes.n_cols);
    

    // parent type
    fieldNum = mxGetFieldNumber(prhs[5],"type");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_type",
                           "parent type field not found");
    // fieldPtr = mxGetFieldByNumber(prhs[5],0,fieldNum);  
    // if (!mxIsDouble(fieldPtr))
    //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_type_not_double",
    //                        "parent type must be of type double");
    
    int lastIGFEM = -1;
    for (std::size_t j = 0; j < elemNodes.n_cols; j++)
    {
        fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum);
        if (fieldPtr && !mxIsEmpty(fieldPtr) && mxGetScalar(fieldPtr) == 2)
        {
            parents(j).type = IGFEM;
            lastIGFEM = j;
        }
    }
;
    // parent node
    std::size_t estSize = 0; // estimated size of ii,jj and Kval
     // number of entries of the stiffness matrix of a regular element
    std::size_t KelRegSize = elemNodes.n_rows*elemNodes.n_rows;
    fieldNum = mxGetFieldNumber(prhs[5],"nodes");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_node",
                           "parent node field not found");
    // fieldPtr = mxGetFieldByNumber(prhs[5],0,fieldNum);  
    // if (!mxIsDouble(fieldPtr))
    //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_node_not_double",
    //                        "parent node must be of type double");
    
    for (std::size_t j = 0; j < elemNodes.n_cols; j++)
    {
        fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum);
        
        if (fieldPtr && !mxIsEmpty(fieldPtr) && parents(j).type == IGFEM )
        {
            parents(j).nodes = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
            estSize += parents(j).nodes.n_elem*parents(j).nodes.n_elem;
        }
        else
        {            
            estSize += KelRegSize;
        }
    }
    
    // channelNum
    fieldNum = mxGetFieldNumber(prhs[5],"channelNum");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_channelNum",
                           "parent channelNum field not found");
    //fieldPtr = mxGetFieldByNumber(prhs[5],0,fieldNum);  
    //if (!mxIsDouble(fieldPtr))
    //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_channelNum_not_double",
    //                        "parent channelNum must be of type double");
    for (std::size_t j = 0; j < elemNodes.n_cols; j++)
    {
        fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum);
        if (fieldPtr && !mxIsEmpty(fieldPtr))
            parents(j).channelNum = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
    }
    
    // channelNodes
    fieldNum = mxGetFieldNumber(prhs[5],"channelNodes");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_channelNodes",
                           "parent(i).channelNodes field not found");
    //fieldPtr = mxGetFieldByNumber(prhs[5],0,fieldNum);  
    //if (!mxIsDouble(fieldPtr))
    //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_channelNodes_not_double",
    //                        "parent(i).channelNodes must be of type double");
    
    for (std::size_t j = 0; j < elemNodes.n_cols; j++)
    {
        fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum);
        if (fieldPtr && !mxIsEmpty(fieldPtr))
        {
            parents(j).channelNodes = arma::conv_to<arma::umat>::from(armaGetPr(fieldPtr)-1);
            if(parents(j).channelNodes.n_rows != 2)
                mexErrMsgIdAndTxt("mx_assemble_sparse:parent_channelNodes_no_2_rows",
                            "parent(i).channelNodes must be a 2 x number of channelNodes matrix");
        }
    }
    
    
    // nSharedElems
    //fieldNum = mxGetFieldNumber(prhs[5],"nSharedElems");
    //if (fieldNum < 0)
    //     mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_nSharedElems",
    //                       "parent(i).nSharedElems field not found");
    //fieldPtr = mxGetFieldByNumber(prhs[5],0,fieldNum);  
    //if (!mxIsDouble(fieldPtr))
    //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_nSharedElems_not_double",
    //                        "parent(i).nSharedElems must be of type double");

    //for (std::size_t j = 0; j < elemNodes.n_cols; j++)
    //{
    //    fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum);
    //    if (fieldPtr && !mxIsEmpty(fieldPtr))
    //    {
    //        parents(j).nSharedElems = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr));
    //    }
    //}
    
    // cstrLocNodes
    fieldNum = mxGetFieldNumber(prhs[5],"cstrLocNodes");
    if (fieldNum < 0)
         mexWarnMsgIdAndTxt("mx_assemble_sparse:no_parent_cstrLocNodes",
                           "parent(i).cstrLocNodes field not found");
    else
    {
        int fieldNum2 = mxGetFieldNumber(prhs[5],"cstrVals");
        if (fieldNum2 < 0)
            mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_cstrVals",
                              "parent(i).cstrVals field not found");   

        int fieldNum3 = mxGetFieldNumber(prhs[5],"cstrRows");
        if (fieldNum3 < 0)
            mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_cstrRows",
                              "parent(i).cstrVals field not found");   


        for (std::size_t j = 0; j < elemNodes.n_cols; j++)
        {
            fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum);
            if (fieldPtr && !mxIsEmpty(fieldPtr))
            {
                parents(j).cstrLocNodes = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
                estSize += 2*parents(j).cstrLocNodes.n_elem*parents(j).nodes.n_elem + 1;
            }
            fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum2);
            if (fieldPtr && !mxIsEmpty(fieldPtr))
                parents(j).cstrVals = armaGetPr(fieldPtr);
            
            fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum3);            
            if (fieldPtr && !mxIsEmpty(fieldPtr))   
                parents(j).cstrRows = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
             
        }
    }
   
    // conductivity
    fieldNum = mxGetFieldNumber(prhs[5],"conductivity");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_conductivity",
                           "parent(i).conductivity field not found");
    
    
    for (std::size_t j = 0; j < elemNodes.n_cols; j++)
    {
        fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum);
        if (fieldPtr && !mxIsEmpty(fieldPtr))
            parents(j).conductivity = mxGetScalar(fieldPtr);
       
    }
    
    // parent children
    mxArray* chFieldPtr;
    int chFieldNum[nChFields];
    mexPrintf("validating and getting children \n");
    fieldNum = mxGetFieldNumber(prhs[5],"child");
    // if child info is not available, assume all elements are regular
    if (fieldNum < 0 || lastIGFEM < 0)
    {
       
         mexWarnMsgIdAndTxt("mx_assemble_sparse:no_parent_child",
                            "parent.child field not found, all elements are assumed regular");
        if (lastIGFEM >= 0)
            for (std::size_t j = 0; j < elemNodes.n_cols; j++)
                parents(j).type = REGULAR;       
    }
    // get child info if field is available
    else
    {
        fieldPtr = mxGetFieldByNumber(prhs[5],lastIGFEM,fieldNum);  
        if (!fieldPtr || !mxIsStruct(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:parent_child_not_struct",
                                "parent.child is not initialized or not of type struct");
        // validate child fields
        // isTriangle
        chFieldNum[ISTRIANGLE] = mxGetFieldNumber(fieldPtr,"isTriangle");
        if (chFieldNum[ISTRIANGLE] < 0)
             mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_isTriangle",
                                "parent.child.isTriangle field not found");
        // locPaNodes
        chFieldNum[LOCPANODES] = mxGetFieldNumber(fieldPtr,"locPaNodes");
        if (chFieldNum[LOCPANODES] < 0)
            mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_locPaNodes",
                               "parent.child.locPaNodes field not found");
        //chFieldPtr = mxGetFieldByNumber(fieldPtr,0,chFieldNum[LOCPANODES]);   
        //if (!mxIsDouble(chFieldPtr))
        //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_child_locPaNodes_not_double",
        //                        "parent.child.locPaNodes must be of type double");

        // locPaEnNodes
        chFieldNum[LOCPAENNODES] = mxGetFieldNumber(fieldPtr,"locPaEnNodes");
        if (chFieldNum[LOCPAENNODES] < 0)
            mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_locPaEnNodes",
                               "parent.child.locPaEnNodes field not found"); 
        //chFieldPtr = mxGetFieldByNumber(fieldPtr,0,chFieldNum[LOCPAENNODES]);   
        //if (!mxIsDouble(chFieldPtr))
        //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_child_locPaEnNodes_not_double",
        //                        "parent.child.locPaEnNodes must be of type double");

        // locEnNodes
        chFieldNum[LOCENNODES] = mxGetFieldNumber(fieldPtr,"locEnNodes");
        if (chFieldNum[LOCENNODES] < 0)
            mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_locEnNodes",
                               "parent.child.locEnNodes field not found");
        //chFieldPtr = mxGetFieldByNumber(fieldPtr,0,chFieldNum[LOCENNODES]);   
        //if (!mxIsDouble(chFieldPtr))
        //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_child_locEnNodes_not_double",
        //                        "parent.child.locEnNodes must be of type double");

        // channelNum
        chFieldNum[CHANNELNUM] = mxGetFieldNumber(fieldPtr,"channelNum");
        if (chFieldNum[CHANNELNUM] < 0)
            mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_channelNum",
                               "parent.child.channelNum field not found");
        //chFieldPtr = mxGetFieldByNumber(fieldPtr,0,chFieldNum[CHANNELNUM]);   
        //if (!mxIsDouble(chFieldPtr))
        //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_child_channelNum_not_double",
        //                        "parent.child.channelNum must be of type double");

        // channelNodes
        chFieldNum[CHANNELNODES] = mxGetFieldNumber(fieldPtr,"channelNodes");
        if (chFieldNum[CHANNELNODES] < 0)
            mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_channelNodes",
                               "parent.child.channelNodes field not found");
        //chFieldPtr = mxGetFieldByNumber(fieldPtr,0,chFieldNum[CHANNELNODES]);   
        //if (!mxIsDouble(chFieldPtr))
        //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_child_channelNodes_not_double",
        //                        "parent.child.channelNodes must be of type double");
        
        chFieldNum[CHANNEL_LOC_NODES] = mxGetFieldNumber(fieldPtr,"channelLocNodes");
        if (chFieldNum[CHANNEL_LOC_NODES] < 0)
            mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_channelLocNodes",
                               "parent.child.channelLocNodes field not found");

        // channelNurbsParam
        chFieldNum[CHANNELNURBSPARAM] = mxGetFieldNumber(fieldPtr,"channelNurbsParam");
        if (chFieldNum[CHANNELNURBSPARAM] < 0)
            mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_channelNurbsParam",
                               "parent.child.channelNurbsParam field not found");
        //chFieldPtr = mxGetFieldByNumber(fieldPtr,0,chFieldNum[CHANNELNURBSPARAM]);   
        //if (!mxIsDouble(chFieldPtr))
        //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_child_channelNurbsParam_not_double",
        //                        "parent.child.channelNurbsParam must be of type double");
        //conductivity
        //chFieldNum[CONDUCTIVITY] = mxGetFieldNumber(fieldPtr,"conductivity");
        //if (chFieldNum[CONDUCTIVITY] < 0)
        //    mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_conductivity",
        //                       "parent.child.conductivity field not found");
        //chFieldPtr = mxGetFieldByNumber(fieldPtr,0,chFieldNum[CONDUCTIVITY]);   
        //if (!mxIsDouble(chFieldPtr))
        //      mexErrMsgIdAndTxt("mx_assemble_sparse:parent_child_conductivity_not_double",
        //                       "parent.child.conductivity must be of type double");

       
        //nSharedChildren
        //chFieldNum[NSHAREDCHILDREN] = mxGetFieldNumber(fieldPtr,"nSharedChildren");
        //if (chFieldNum[NSHAREDCHILDREN] < 0)
        //    mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_child_nSharedChildren",
        //                      "parent.child.nSharedChildren field not found");

        std::size_t nChilds;
        for (std::size_t j = 0; j < elemNodes.n_cols; j++)
        {
            fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum);
            if (fieldPtr)
            {
                nChilds = mxGetNumberOfElements(fieldPtr);
                parents(j).children.set_size(nChilds);
            }
            else
                nChilds = 0;
            for (std::size_t k = 0; k < nChilds; k++)
            {
                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[ISTRIANGLE]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                {
                    if (mxIsLogicalScalarTrue(chFieldPtr))
                        parents(j).children(k).shape = TRIANGLE;
                    else
                        parents(j).children(k).shape = QUADRILATERAL;
                }
                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[LOCPANODES]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                    parents(j).children(k).locPaNodes = arma::conv_to<arma::uvec>::from(armaGetPr(chFieldPtr)-1);
                
                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[LOCPAENNODES]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                    parents(j).children(k).locPaEnNodes = arma::conv_to<arma::uvec>::from(armaGetPr(chFieldPtr)-1);

                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[LOCENNODES]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                    parents(j).children(k).locEnNodes = arma::conv_to<arma::uvec>::from(armaGetPr(chFieldPtr)-1);

                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[CHANNELNUM]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                    parents(j).children(k).channelNum = arma::conv_to<arma::uvec>::from(armaGetPr(chFieldPtr)-1);                

                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[CHANNELNODES]);
                //if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                   // parents(j).children(k).channelNodes = arma::conv_to<arma::umat>::from(armaGetPr(chFieldPtr)-1);
              
                
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                {
                    std::size_t nChannels = mxGetNumberOfElements(chFieldPtr);
                    mxArray *cellElemPtr;
                    parents(j).children(k).channelNodes.set_size(2,nChannels);
                    arma::uvec channelNodes;
                    for (std::size_t i = 0; i < nChannels; i++)
                    {
                        cellElemPtr = mxGetCell(chFieldPtr,i);
                        if (cellElemPtr && !mxIsEmpty(cellElemPtr))
                        {
                            channelNodes = arma::conv_to<arma::uvec>::from(armaGetPr(cellElemPtr)-1);
                            parents(j).children(k).channelNodes(0,i) = channelNodes(0);
                            parents(j).children(k).channelNodes(1,i) = channelNodes(channelNodes.n_elem-1);
                        }
                    }
                }

                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[CHANNEL_LOC_NODES]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                {
                    parents(j).children(k).channelLocNodes 
                        = arma::conv_to<arma::umat>::from(armaGetPr(chFieldPtr)-1);
                }

                chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[CHANNELNURBSPARAM]);
                if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                    parents(j).children(k).channelNurbsParam = arma::conv_to<arma::mat>::from(armaGetPr(chFieldPtr));

                //chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[CONDUCTIVITY]);
                //if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                //    parents(j).children(k).conductivity = mxGetScalar(chFieldPtr);
                
                //chFieldPtr = mxGetFieldByNumber(fieldPtr,k,chFieldNum[NSHAREDCHILDREN]);
                //if(chFieldPtr && !mxIsEmpty(chFieldPtr)) 
                //    parents(j).children(k).nSharedChildren = arma::conv_to<arma::uvec>::from(armaGetPr(chFieldPtr));

            }
        }
    } // get child info if field is available
    
 
        
    // validating channel data
    mexPrintf("validating and getting channel data \n");
    fieldNum = mxGetFieldNumber(prhs[6],"mcf");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_channel_mcf",
                           "channel mcf field not found");
    fieldPtr = mxGetFieldByNumber(prhs[6],0,fieldNum);
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:channel_mcf_empty",
                             "channel mcf field empty");   
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:channel_mcf_not_double",
                            "channel must be of type double"); 
 
    chanNetwork channels;
    channels.mcf = armaGetPr(fieldPtr); 
    
    fieldNum = mxGetFieldNumber(prhs[6],"model");
    if (fieldNum < 0)
         mexErrMsgIdAndTxt("mx_assemble_sparse:no_channel_model",
                           "channel model field not found");
    fieldPtr = mxGetFieldByNumber(prhs[6],0,fieldNum);  
    if(!fieldPtr || mxIsEmpty(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:channel_model_empty",
                             "channel model field empty");   
    if (!mxIsDouble(fieldPtr))
          mexErrMsgIdAndTxt("mx_assemble_sparse:channel_model_not_double",
                            "channel model must be of type double"); 

    
    int modelNum = mxGetScalar(fieldPtr); 
    if (modelNum  == 1)
        channels.model = MEAN_TEMP;
    else if (modelNum == 2)
        channels.model = CONST_HEAT;
     
      // channelNurbsParam
    if (channels.model == CONST_HEAT)
    {
        mexPrintf("validating and getting parent channel NURBS parameters \n");
        fieldNum = mxGetFieldNumber(prhs[5],"channelNurbsParam");
        if (fieldNum < 0)
             mexErrMsgIdAndTxt("mx_assemble_sparse:no_parent_channelNurbsParam",
                               "parent(i).channelNurbsParam must be provided for constant heat flux model");


        for (std::size_t j = 0; j < elemNodes.n_cols; j++)
        {
            fieldPtr = mxGetFieldByNumber(prhs[5],j,fieldNum);
            if (fieldPtr && !mxIsEmpty(fieldPtr))
            {
                parents(j).channelNurbsParam =  armaGetPr(fieldPtr);
                if(parents(j).channelNurbsParam.n_rows != 2)
                    mexErrMsgIdAndTxt("mx_assemble_sparse:parent_channelNodes_no_2_rows",
                                "parent(i).channelNodes must be a 2 x number of channelNodes matrix");
            }
        }
        
        mexPrintf("validating and getting const heat flux model parameters \n");
        fieldNum = mxGetFieldNumber(prhs[6],"Tin");
        if (fieldNum < 0)
             mexErrMsgIdAndTxt("mx_assemble_sparse:no_channel_Tin",
                               "channel Tin field not found");
        fieldPtr = mxGetFieldByNumber(prhs[6],0,fieldNum);
        if(!fieldPtr || mxIsEmpty(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_Tin_empty",
                                 "channel Tin field empty");            
        if (!mxIsDouble(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_Tin_not_double",
                                "channel Tin must be of type double"); 
        channels.Tin = mxGetScalar(fieldPtr);

        fieldNum = mxGetFieldNumber(prhs[6],"kapf");
        if (fieldNum < 0)
             mexErrMsgIdAndTxt("mx_assemble_sparse:no_channel_kapf",
                               "channel kapf field not found");
        fieldPtr = mxGetFieldByNumber(prhs[6],0,fieldNum);  
        if(!fieldPtr || mxIsEmpty(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_kapf_empty",
                                 "channel kapf field empty"); 
        if (!mxIsDouble(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_kapf_not_double",
                                "channel kapf must be of type double"); 
        channels.kapf = armaGetPr(fieldPtr);
        
        fieldNum = mxGetFieldNumber(prhs[6],"length");
        if (fieldNum < 0)
             mexErrMsgIdAndTxt("mx_assemble_sparse:no_channel_length",
                               "channel length field not found");
        fieldPtr = mxGetFieldByNumber(prhs[6],0,fieldNum);  
        if(!fieldPtr || mxIsEmpty(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_modelType_empty",
                                 "channel length field empty"); 
        if (!mxIsDouble(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_length_not_double",
                                "channel length must be of type double"); 
        channels.lengths = armaGetPr(fieldPtr);
        
        fieldNum = mxGetFieldNumber(prhs[6],"eigvalsq");
        if (fieldNum < 0)
             mexErrMsgIdAndTxt("mx_assemble_sparse:no_channel_eigvalsq",
                               "channel eigvalsq field not found");
        fieldPtr = mxGetFieldByNumber(prhs[6],0,fieldNum);  
        if(!fieldPtr || mxIsEmpty(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_eigvalsq_empty",
                                 "channel eigvalsq field empty"); 
        if (!mxIsDouble(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_eigvalsq_not_double",
                                "channel eigvalsq must be of type double"); 
   
        channels.eigvalsq = armaGetPr(fieldPtr);
        
        fieldNum = mxGetFieldNumber(prhs[6],"CR1s");
        if (fieldNum < 0)
             mexErrMsgIdAndTxt("mx_assemble_sparse:no_channel_CR1s",
                               "channel CR1s field not found");
        fieldPtr = mxGetFieldByNumber(prhs[6],0,fieldNum);  
        if(!fieldPtr || mxIsEmpty(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_CR1s_empty",
                                 "channel CR1s field empty"); 
        if (!mxIsDouble(fieldPtr))
              mexErrMsgIdAndTxt("mx_assemble_sparse:channel_CR1s_not_double",
                                "channel CR1s must be of type double"); 
   
        channels.CR1s = armaGetPr(fieldPtr);
        
    }    
    
    // get Dirichlet nodes
   
    mexPrintf("validating and getting Dirichlet nodes and values \n");
    //Dirichlet dirichlet;
    std::size_t nDirichletNodes = 0;
    fieldNum = mxGetFieldNumber(prhs[7],"temp_node");
    if (fieldNum < 0)
         mexWarnMsgIdAndTxt("mx_assemble_sparse:no_Dirichlet_node",
                            "Dricihlet temp_node field not found");
    else
    {
        fieldPtr = mxGetFieldByNumber(prhs[7],0,fieldNum);        
        if (fieldPtr && !mxIsEmpty(fieldPtr))
            //dirichlet.nodes = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
            nDirichletNodes = mxGetNumberOfElements(fieldPtr);
        else
            mexWarnMsgIdAndTxt("mx_assemble_sparse:no_Dirichlet_node",
                               "Diricihlet temp_node field not found");
    }

    /* don't need this since we are not solving the system of equations    
    fieldNum = mxGetFieldNumber(prhs[7],"temp_value");
    if (fieldNum < 0)
         mexWarnMsgIdAndTxt("mx_assemble_sparse:no_Dirichlet_value",
                            "Diricihlet temp_value field not found");
    else
    {
        fieldPtr = mxGetFieldByNumber(prhs[7],0,fieldNum);      
        if (fieldPtr && !mxIsEmpty(fieldPtr) )
            dirichlet.vals = armaGetPr(fieldPtr);
        else
            mexWarnMsgIdAndTxt("mx_assemble_sparse:no_Dirichlet_value",
                               "Diricihlet temp_value field not found");
    }    
    */

    mexPrintf("validating and getting Neumann info \n");
    Neumann neumann;
    fieldNum = mxGetFieldNumber(prhs[8],"heatFlux_elem");
    if (fieldNum > 0)
    {
        // heatFlux_elem
        fieldPtr = mxGetFieldByNumber(prhs[8],0,fieldNum);        
        if (fieldPtr && !mxIsEmpty(fieldPtr))
        {
            neumann.elem = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
            
            // heatFlux_surface
            fieldNum = mxGetFieldNumber(prhs[8],"heatFlux_surface");
            if (fieldNum < 0)
                mexErrMsgIdAndTxt("mx_assemble_sparse:heatFlux_surface",
                                   "neumann heatFlux_surface field not found");
            else
            {
                fieldPtr = mxGetFieldByNumber(prhs[8],0,fieldNum);
                if (fieldPtr && !mxIsEmpty(fieldPtr))
                    neumann.surf = arma::conv_to<arma::uvec>::from(armaGetPr(fieldPtr)-1);
                else
                    mexErrMsgIdAndTxt("mx_assemble_sparse:heatFlux_surface",
                                      "neumann heatFlux_surface field not found");
            }

            // heatFlux_value
            fieldNum = mxGetFieldNumber(prhs[8],"heatFlux_value");
            if (fieldNum < 0)
                mexErrMsgIdAndTxt("mx_assemble_sparse:heatFlux_value",
                                   "neumann heatFlux_value field not found");
            else
            {
                fieldPtr = mxGetFieldByNumber(prhs[8],0,fieldNum);
                if (fieldPtr && !mxIsEmpty(fieldPtr))
                    neumann.vals = armaGetPr(fieldPtr);
                else
                    mexErrMsgIdAndTxt("mx_assemble_sparse:heatFlux_value",
                                      "neumann heatFlux_value field not found");
            }
        }
    }


    arma::ivec ii;
    arma::ivec jj;
    arma::vec Kval;
    arma::vec Pglo;
    
    mexPrintf("estimated size of ii, jj or Kval = %i \n", estSize);
    mexPrintf("assemble \n");
    assemble(ii,
             jj,
             Kval,
             Pglo,
             estSize,
             nodeCoords,   
             elemNodes,
             elemHeatSource,                 
             eqNum,
             gauss1,
             parents,
             channels,
             neumann,
             supg);
    

    arma::umat ffInd, fpInd, pfInd, ppInd;
    arma::vec ffKval, fpKval, pfKval, ppKval, Pf, Pp;
    
    std::size_t iinelem = ii.n_elem;
    //std::size_t estNff = elemNodes.n_cols*5*5;
    //estNff = std::min(estNff,iinelem);
    std::size_t estNff = iinelem;
    std::size_t estNfp = nDirichletNodes*8*(elemNodes.n_rows-1);
    estNfp = std::min(estNfp,iinelem);
    std::size_t estNpf = estNfp;
    std::size_t estNpp = ceil(1.4*estNfp);
    estNpp = std::min(estNpp,iinelem);

    mexPrintf("estimated nff = %i, nfp = %i, npf = %i, npp = %i \n",
              estNff,estNfp,estNpf,estNpp);
    partition_stiffness_mat_n_load_vec(ffInd,
                                       ffKval,
                                       fpInd,
                                       fpKval,
                                       pfInd,
                                       pfKval,
                                       ppInd,
                                       ppKval,
                                       Pf,
                                       Pp,
                                       ii,
                                       jj,
                                       Kval,
                                       Pglo,
                                       eqNum,
                                       estNff,
                                       estNfp,
                                       estNpf,
                                       estNpp);

    mexPrintf("actual nff = %i, nfp  = %i, npf = %i, npp = %i \n",
              ffKval.n_elem,fpKval.n_elem,pfKval.n_elem,ppKval.n_elem);
   
    
    ii.clear();
    jj.clear();
    Kval.clear();
    Pglo.clear();
    
    /*       
    arma::mat ffIndK = arma::join_rows(arma::conv_to<arma::mat>::from(ffInd.t()+1),ffKval);
    arma::mat fpIndK = arma::join_rows(arma::conv_to<arma::mat>::from(fpInd.t()+1),fpKval);
    arma::mat pfIndK = arma::join_rows(arma::conv_to<arma::mat>::from(pfInd.t()+1),pfKval);
    arma::mat ppIndK = arma::join_rows(arma::conv_to<arma::mat>::from(ppInd.t()+1),ppKval); 
    */
    
    std::size_t nDof = eqNum.n_elem - nDirichletNodes; 
    arma::sp_mat Kff = sp_mat(true,ffInd,ffKval,nDof,nDof);
    arma::sp_mat Kfp = sp_mat(true,fpInd,fpKval,nDof,nDirichletNodes);
    arma::sp_mat Kpf = sp_mat(true,pfInd,pfKval,nDirichletNodes,nDof);
    arma::sp_mat Kpp = sp_mat(true,ppInd,ppKval,nDirichletNodes,
                                                nDirichletNodes);
    
    ffInd.clear();
    ffKval.clear();
    fpInd.clear();
    fpKval.clear();
    pfInd.clear();
    pfKval.clear();
    ppInd.clear();
    ppKval.clear();
    
    plhs[0] = armaCreateMxSparseMatrix(nDof,nDof,Kff.n_nonzero);
    plhs[1] = armaCreateMxSparseMatrix(nDof,nDirichletNodes,Kfp.n_nonzero);
    plhs[2] = armaCreateMxSparseMatrix(nDirichletNodes,nDof,Kpf.n_nonzero);
    plhs[3] = armaCreateMxSparseMatrix(nDirichletNodes,nDirichletNodes,
                                       Kpp.n_nonzero);
    plhs[4] = armaCreateMxMatrix(Pf.n_elem,1);
    plhs[5] = armaCreateMxMatrix(Pp.n_elem,1);
    
    armaSetSparsePr(plhs[0],Kff);
    armaSetSparsePr(plhs[1],Kfp);
    armaSetSparsePr(plhs[2],Kpf);
    armaSetSparsePr(plhs[3],Kpp);
    armaSetPr(plhs[4],Pf);
    armaSetPr(plhs[5],Pp);
    
    
    
    /*
    plhs[6] = armaCreateMxMatrix(ffIndK.n_rows,ffIndK.n_cols);
    plhs[7] = armaCreateMxMatrix(fpIndK.n_rows,fpIndK.n_cols);
    plhs[8] = armaCreateMxMatrix(pfIndK.n_rows,pfIndK.n_cols);
    plhs[9] = armaCreateMxMatrix(ppIndK.n_rows,ppIndK.n_cols);
    //plhs[4] = armaCreateMxMatrix(Pf.n_elem,1);
    //plhs[5] = armaCreateMxMatrix(Pp.n_elem,1);
    armaSetPr(plhs[6],ffIndK);
    armaSetPr(plhs[7],fpIndK);
    armaSetPr(plhs[8],pfIndK);
    armaSetPr(plhs[9],ppIndK);
    //armaSetPr(plhs[4],Pf);
    //armaSetPr(plhs[5],Pp);
    */
}
    
