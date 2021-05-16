#include "assemble.h"
#include "armadillo"

namespace igfem
{
double conductance(double x, double kapf, double mcf, const arma::vec& eigvalsq, const arma::vec& CR1s)
{
    const double frac = 11.0/24.0;
    double factor = 0.5*x*arma::datum::pi*kapf/mcf;
    return 2.0*arma::datum::pi*kapf/(4*factor + frac + arma::sum(CR1s%arma::exp(-factor*eigvalsq))); 
}
}
