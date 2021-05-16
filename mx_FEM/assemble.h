/*
Created by Marcus Tan on 8/7/2014
Updated on 1/29/2016
Copyright 2014 University of Illinois 
Purpose: header file for the assembly of stiffness matrix
*/
#ifndef ASSEMBLE_H_
#define ASSEMBLE_H_
#include <cstddef> // NULL, std::size_t
#include "armadillo"

#define JACTOL 1e-13 // tolerance below which Jacobian is considered zero
#define DENOMTOL 1e-13 // tolerance below which denominator is considered 0
#define MCFTOL 1e-10 // tolerance below which mcf is taken to be 0

namespace igfem
{



// types of elements
enum shapes {TRIANGLE,QUADRILATERAL};
struct child
{
    shapes shape; // child element type
    arma::uvec locPaNodes; //local numbers of child element nodes relative to parent nodes
    arma::uvec locPaEnNodes; //local numbers of child element enrichment nodes rel to parent nodes
    arma::uvec locEnNodes; //local numbers of child element enrichment nodes rel to itself
    arma::uvec channelNum;
    arma::umat channelNodes; // global node number of channel end points
    arma::umat channelLocNodes; //local node number of channel end points wrt child 
    arma::mat channelNurbsParam;    
    //double conductivity;
};

enum parentType {REGULAR,IGFEM};
struct parent
{
    parent() : type(REGULAR){};
    parentType type;
    arma::uvec nodes; 
    arma::uvec channelNum; // label of channels. leave empty if no channel.
    arma::umat channelNodes; // 2 x number of channel matrix of node labels of channel end points
    //arma::field<arma::uvec> channelNodes; 
    arma::mat channelNurbsParam; // 2 x number of channel nurbs parameters of channel end points
    arma::uvec nSharedElems; // number of elements sharing each channel 
    arma::uvec cstrLocNodes;
    arma::vec cstrVals;
    arma::uvec cstrRows;
    arma::field<child> children;
    double conductivity;
};

struct Dirichlet
{
    arma::uvec nodes;
    arma::vec vals;
};


struct Neumann
{
    Neumann()
    {
        surfLocNodes.set_size(2,3);
        surfLocNodes << 0 << 1 << 2 << arma::endr
                     << 1 << 2 << 0 << arma::endr;
    };
    arma::umat surfLocNodes;
    arma::uvec elem;
    arma::uvec surf;
    arma::vec vals;
};


// mean temperature and constant heat flux models
enum modelType {MEAN_TEMP, CONST_HEAT};

// function used only for CONST_HEAT model
double conductance(double x, double kapf, double mcf, const arma::vec& eigvalsq, const arma::vec& CR1s); 

struct chanNetwork
{
    modelType model;
    arma::vec mcf; // vector of mass flow rate x heat capacity of all channels
    // following variables are only used for the constant heat flux model
    double Tin; // channel inlet temperature
    arma::vec kapf; // thermal conductivity of fluid
    arma::vec eigvalsq;
    arma::vec CR1s;
    arma::vec lengths; // lengths of all channels
};

struct gauss
{
    arma::mat elem; // a (number of dimensions+1)x number of gauss points matrix for integration over elements
    arma::mat quadElem;
    arma::mat line;  // a 2 x number of gauss points matrix for integration over lines. 
};

struct Convection
{
    double coef;
    double Tref;
};

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
              bool supg);

template <typename T> 
inline int sgn(T val) {
        return (T(0) < val) - (val < T(0));
}

void compute_regular_element(arma::mat& Kel, 
                             arma::vec& Pel, 
                             const arma::mat& nodeCoords,
                             const arma::uvec& elemNodes, 
                             double elemHeatSource,
                             const Convection& convect,
                             const gauss& gauss1,
                             const parent& parent1,
                             const chanNetwork& channels,
                             bool supg);

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
                           bool supg);

// for supg calculation
void streamwise_elem_length(double& he, 
                            arma::vec& Bsw1, 
                            const arma::vec& channelVec,
                            const arma::mat& B);

void supg_weighing_function(arma::vec& W,
                            const arma::mat& channelUnitVecs,
                            const arma::vec& N,
                            const arma::mat& B,
                            const arma::vec& mcf);

void shape_function_3D(arma::vec& N, arma::mat& DN, const arma::vec& locCoord); // tetrahedral shape function

void shape_function_2D(arma::vec& N, arma::mat& DN, const arma::vec& locCoord); // triangular shape function

void shape_function_2D(arma::vec& N, arma::mat& DN, const arma::vec& locCoord, 
                       shapes shape); // triangular or quadrilateral shape functions

void shape_function_1D(arma::vec& N, arma::mat& DN, double locCoord); // 1D shape function

void child_element_NBJ(arma::vec& N, 
                       arma::mat& B, 
                       double& detJ,
                       arma::vec& X,
                       bool& calcBJch,
                       arma::mat& Bch,
                       bool& calcBel,
                       arma::mat& Bel,
                       std::size_t nNodes,
                       const arma::uvec& paLocOrigNodes,
                       const arma::uvec& chLocPaEnrichNodes,
                       const arma::uvec& chLocEnrichNodes,
                       const arma::mat& Xpa,
                       const arma::mat& Xch,
                       const arma::vec& locCoord,
                       shapes shape);

arma::vec local_coord_2D(const arma::vec& X, const arma::mat& Xel);

arma::vec local_coord_2D(const arma::vec& X, const arma::mat& Xel, shapes shape);

arma::vec local_coord_2D_along_edge(const arma::vec& X, const arma::mat& Xel, shapes shape, int edge);

int edge_number(const arma::uvec& twoNodes, shapes shape);

template<class T>                             
void combination_two_vectors(T& uout, 
                             T& vout, 
                             const T& uin,
                             const T& vin);       
#include "combination_two_vectors.tcc"

void impose_Neumann_BC(arma::vec& Pf,
                       const Neumann& neumann,
                       const arma::umat& elemNodes,
                       const arma::mat& nodeCoords);

//void impose_Dirichlet_BC(arma::vec& Up,
//                         const arma::uvec& dirichletNodes,
//                         const arma::vec& dirichletVals,
//                         const arma::ivec);

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
                                        std::size_t estNpp);


void solve_mat_eqn(arma::vec UUR, arma::vec PUR, 
                   const arma::sp_mat& Kff,
                   const arma::sp_mat& Kfp,
                   const arma::sp_mat& Kpf,
                   const arma::sp_mat& Kpp,
                   const arma::vec& Pf,
                   const arma::vec& Pp,
                   const arma::vec& Up,
                   const arma::ivec& eqNum);

double body_source_straight_channel(const arma::vec& X);
double body_source_semicircular_channel(const arma::vec& X);
double body_source(const arma::vec& X);
}
#endif
