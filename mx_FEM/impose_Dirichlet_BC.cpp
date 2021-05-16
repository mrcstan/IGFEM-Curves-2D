#include "assemble.h"
#include "armadillo"

namespace igfem
{
    void impose_Dirichlet_BC(arma::vec& Up, 
                            const arma::uvec& dirichletNodes,
                            const arma::vec& dirichletVals,
                            const arma::ivec& eqNum)
    {
        Up.set_size(dirichletNodes.n_elem);
        Up(arma::conv_to<arma::uvec>::from(-eqNum(dirichletNodes)-1)) = dirichletVals;
    }
}
