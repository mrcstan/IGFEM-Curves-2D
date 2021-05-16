/*
Created by Marcus Tan on 1/11/2014
Modified on 1/26/2016
Copyright 2014 University of Illinois 
Purpose: this function calculates the body source at a give position
*/
#include "assemble.h"
#include "armadillo"
#include <iostream>
#include <math.h>

#define X0Y0 0.2489
#define X1Y0 23.42
#define X0Y1 21.87
#define X2Y0 -485.7
#define X1Y1 -16.96
#define X0Y2 -356.8
#define X3Y0 4317
#define X2Y1 170.4
#define X1Y2 276.7
#define X0Y3 2422
#define X4Y0 -14060
#define X3Y1 -791
#define X2Y2 -490
#define X1Y3 -972.2
#define X0Y4 -5849

// only uncomment one of these to activate
// the needed distributed heat source
//#define SINGLE_CHANNEL
//#define SEMICIRCULAR_CHANNEL
//#define CROSS_CHANNELS
#define FOR_VALIDATION
//#define GRIDLIKE_CHANNELS

// To use any of the distributed heat sources in this source code,
// the directive SOURCE_FUNCTION
// must be uncommented in the following source codes
//      compute_regular_element.cpp
//      compute_IGFEM_element.cpp
namespace igfem
{
// distributed heat source for single channel manufactured solution 
#ifdef SINGLE_CHANNEL
double body_source(const arma::vec& X)
{
   const double L = 0.1;
   const double xo = 0.05;
   const double mcf = 1.0;
   const double lambda = L/(mcf*xo*(L-xo));
   const double A1 = 1.0/xo;
   const double A2 = 1.0/(xo-L);
   
   if (X(0) <= xo)
       return -lambda*lambda*A1*X(0)*exp(-lambda*X(1));
   else
       return -lambda*lambda*A2*(X(0)-L)*exp(-lambda*X(1));
}
#endif

// distributed heat source for semicircular channel manufactured solution
// Note: See  "A NURBS-based IGFEM scheme for thermal analysis
//       and design of microvascular composites" by Marcus H. Y. Tan et al for solution
#ifdef SEMICIRCULAR_CHANNEL
void cart2pol(double& theta, double& radius, double x, double y)
{
    radius = sqrt(x*x + y*y);
    theta = atan2(y,x);
}

double body_source(const arma::vec& X)
{
    const double xc = 0.05, yc = 0.0, ro = 0.04, n = 3, mcf = 90, conductivity = 1.0;
    const double c1 = 1.0/pow(ro,n), c2 = pow(ro,n), 
                alpha = mcf/conductivity, lambda = 2.0*n/alpha;
    double theta, radius;
    cart2pol(theta,radius,X(0)-xc,X(1)-yc);
    if (radius < ro)
        return -(n*n+lambda*lambda)*exp(-lambda*theta)*c1*pow(radius,n-2)/conductivity;
    else
        return -(n*n+lambda*lambda)*exp(-lambda*theta)*c2*pow(radius,-n-2)/conductivity; 

}
#endif

// distributed heat source for manufactured solution with a cross channel
// Note: incomplete (see create_channels.m or "A NURBS-based IGFEM scheme for thermal analysis
//       and design of microvascular composites" by Marcus H. Y. Tan et al for solution
#ifdef CROSS_CHANNELS
double body_source(const arma::vec& X)
{
   const double lam1 = 4.0, lam2 = 20.0, lam3 = 20.0, lam4 = 100.0, xo = 0.05, yo = 0.05;
   const double C = exp(-lam2*xo);
    
   if (X(0) <= xo && X(1) <= yo)
       return -lambda*lambda*A1*X(0)*exp(-lambda*X(1));
   else if (X(0) <= xo && X(1) > yo)
   else if (X(0) > xo && X(1) <= yo)
   else
       return -lambda*lambda*A2*(X(0)-L)*exp(-lambda*X(1));
}
#endif

// uneven heating for Stephen's channels
// See "Gradient-based design of actively-cooled microvascular composite panels" 
// by Marcus Hwai Yik Tan et al for more information

#ifdef FOR_VALIDATION
double body_source(const arma::vec& X)
{

    //return 16000*X(0);
    // Stephen's experiment, non-uniform heat source
	return 500*(X0Y0 + X1Y0*X(0) + X0Y1*X(1) + X2Y0*X(0)*X(0) + X1Y1*X(0)*X(1) + X0Y2*X(1)*X(1)
             + X3Y0*pow(X(0),3.0) + X2Y1*X(0)*X(0)*X(1) + X1Y2*X(0)*X(1)*X(1) + X0Y3*pow(X(1),3.0)
             + X4Y0*pow(X(0),4.0) + X3Y1*pow(X(0),3.0)*X(1) + X2Y2*X(0)*X(0)*X(1)*X(1)
             + X1Y3*X(0)*pow(X(1),3.0) + X0Y4*pow(X(1),4.0));    
   
    // Double localized heat sources
    /*
    const double ro = 0.015;
    const double xo = 0.04;
    const double yo = 0.04;
    const double Qo = 250*0.15*0.2*pow(15/(16*ro),2);
    const double xmin = xo - ro;
    const double xmax = xo + ro;
    const double ymin = yo - ro;
    const double ymax = yo + ro;
    double Q = 0;

    if (X(0) >= xmin && X(0) <= xmax && X(1) >= ymin && X(1) <= ymax)
    {
        Q = Qo*pow(1-pow((X(0)-xo)/ro,2),2)*pow(1-pow((X(1)-yo)/ro,2),2);
    }
    
    const double r1 = 0.015;
    const double x1 = 0.11;
    const double y1 = 0.16;
    const double Q1 = 250*0.15*0.2*pow(15/(16*r1),2);
    const double xmin1 = x1 - r1;
    const double xmax1 = x1 + r1;
    const double ymin1 = y1 - r1;
    const double ymax1 = y1 + r1;
    
    if (X(0) >= xmin1 && X(0) <= xmax1 && X(1) >= ymin1 && X(1) <= ymax1)
    {
        Q += Q1*pow(1-pow((X(0)-x1)/r1,2),2)*pow(1-pow((X(1)-y1)/r1,2),2);
    }
   
    return Q;
    */

    /*
    return 500*(0.2489 + 21.87*X(1) + 23.42*X(0) 
                - 356.8*X(1)*X(1) - 16.96*X(0)*X(1) - 485.7*X(0)*X(0)
                + 2422*pow(X(1),3.0) + 276.7*X(1)*X(1)*X(0) + 170.4*X(1)*X(0)*X(0) + 4317*pow(X(0),3.0)
                - 5849*pow(X(1),4.0) - 972.2*pow(X(1),3.0)*X(0) - 490*X(0)*X(0)*X(1)*X(1) 
                - 791*X(1)*pow(X(0),3.0) - 14060*pow(X(0),4.0));
	*/
}
#endif

#ifdef GRIDLIKE_CHANNELS
double body_source(const arma::vec& X)
{
    // For damage-resilence/redundancy project
	return 6000*(0.188625 + 127*X(0) + 99.54*X(1) 
                - 7583*X(0)*X(0) - 4038*X(0)*X(1) - 5169*X(1)*X(1)
                + 1.919E5*pow(X(0),3.0) + 2.205E5*X(0)*X(0)*X(1) + 6.490E4*X(0)*X(1)*X(1) 
                + 1.228E5*pow(X(1),3.0) - 1.700E6*pow(X(0),4.0) - 5.539E6*pow(X(0),3.0)*X(1) 
                - 1.921E6*X(0)*X(0)*X(1)*X(1) + 1.263E6*X(0)*pow(X(1),3.0) - 2.156E6*pow(X(1),4.0)
                - 2.348E6*pow(X(0),5.0) + 5.443E7*pow(X(0),4.0)*X(1) - 8.962E6*pow(X(0),3.0)*pow(X(1),2.0)
                + 4.017E7*pow(X(0),2.0)*pow(X(1),3.0) - 4.462E7*X(0)*pow(X(1),4.0) + 1.982E7*pow(X(1),5.0));    
   
}
#endif

}
