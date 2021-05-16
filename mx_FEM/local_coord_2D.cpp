// this function finds the local_coord of a global coordinate
// wrt a triangle
#include "assemble.h"
#include "armadillo"
#include <iostream>
#include <math.h> // fabs

#define LOCAL_COORD_QUAD_TOL 2e-16 // tolerance for quad elem global to local coordinates calculation
#define LOC_EDGE_TOL 1e-13 // tolerance below which denominator for local_coord_2D_along_edge is considered zero

namespace igfem
{
arma::vec local_coord_2D(const arma::vec& X, const arma::mat& Xel)
{
    if (X.n_elem != 2)
    {
        std::cerr << "local_coord_2D: X must have 2 elements" << std::endl;
        return arma::zeros<arma::vec>(2);

    }
    if (Xel.n_rows != 2)
    {
        std::cerr << "local_coord_2D: Xel must have 2 rows" << std::endl;
        return arma::zeros<arma::vec>(2);
    }
    if (Xel.n_cols != 3)
    {
        std::cerr << "local_coord_2D: Xel must have 3 columns" << std::endl;
        return arma::zeros<arma::vec>(2);
    }
        
    arma::mat A =
    arma::join_horiz(Xel.col(1)-Xel.col(0),Xel.col(2)-Xel.col(0));

    return arma::solve(A,X-Xel.col(0));
}

// CAUTION: this function is not reliable for finding the local coordinate of a point
// on the edge of a quadrilateral element
arma::vec local_coord_2D(const arma::vec& X, const arma::mat& Xel, shapes shape)
{
    if (X.n_elem != 2)
    {
        std::cerr << "local_coord_2D: X must have 2 elements" << std::endl;
        return arma::zeros<arma::vec>(2);

    }
    if (Xel.n_rows != 2)
    {
        std::cerr << "local_coord_2D: Xel must have 2 rows" << std::endl;
        return arma::zeros<arma::vec>(2);
    }

    if (shape == TRIANGLE)
    {
        if (Xel.n_cols != 3)
        {
            std::cerr << "local_coord_2D: Xel must have 3 columns for tri" << std::endl;
            return arma::zeros<arma::vec>(2);
        }
        arma::mat A =
            arma::join_horiz(Xel.col(1)-Xel.col(0),Xel.col(2)-Xel.col(0));
        return arma::solve(A,X-Xel.col(0));
    }
    else if (shape == QUADRILATERAL)
    {
        if (Xel.n_cols != 4)
        {
            std::cerr << "local_coord_2D: Xel must have 4 columns for quad" << std::endl;
            return arma::zeros<arma::vec>(2);
        }
     
        double a1 = Xel(0,0)-Xel(0,1)+Xel(0,2)-Xel(0,3);
        double a2 = Xel(1,0)-Xel(1,1)+Xel(1,2)-Xel(1,3);

        double b1 = -Xel(0,0)+Xel(0,1)+Xel(0,2)-Xel(0,3);
        double b2 = -Xel(1,0)+Xel(1,1)+Xel(1,2)-Xel(1,3);

        double c1 = -Xel(0,0)-Xel(0,1)+Xel(0,2)+Xel(0,3);
        double c2 = -Xel(1,0)-Xel(1,1)+Xel(1,2)+Xel(1,3);

        double d1 = 4.0*X(0) - arma::sum(Xel.row(0));
        double d2 = 4.0*X(1) - arma::sum(Xel.row(1));

        double ab = a1*b2 - a2*b1;
        double ac = a1*c2 - a2*c1;
        double ad = a1*d2 - a2*d1;
        double bc = b1*c2 - b2*c1;
        double bd = b1*d2 - b2*d1;
        double cd = c1*d2 - c2*d1;
        
        arma::vec Xloc(2);
        if (    fabs(a1*a2*ab*ac) > LOCAL_COORD_QUAD_TOL
            || (fabs(a1) < LOCAL_COORD_QUAD_TOL && fabs(a2*c1) > LOCAL_COORD_QUAD_TOL)
            || (fabs(a2) < LOCAL_COORD_QUAD_TOL && fabs(a1*b2) > LOCAL_COORD_QUAD_TOL))
        {
            double B = -bc-ad;
            double C = -cd;
            double sqrtTerm = sqrt(B*B - 4*ab*C);
            double rt1 = 0.5*(-B + sqrtTerm)/ab;
            double rt2 = 0.5*(-B - sqrtTerm)/ab;
            if (rt1 >= -1.0 && rt1 <= 1.0)
                Xloc(0) = rt1;
            else if (rt2 >= -1.0 && rt2 <= 1.0)
                Xloc(0) = rt2;

             Xloc(1) = (ad - ab*Xloc(0))/ac;   
        }
        else if (fabs(a1*a2) > LOCAL_COORD_QUAD_TOL && fabs(ab) < LOCAL_COORD_QUAD_TOL)
        {
            Xloc(0) = -a1*cd/(b1*ac+a1*ad);
            Xloc(1) = ad/ac;
        }
        else if (fabs(a1*a2) > LOCAL_COORD_QUAD_TOL && fabs(ac) < LOCAL_COORD_QUAD_TOL)
        {
            Xloc(0) = ad/ab;
            Xloc(1) = -a1*bd/(c1*ab+a1*ad);
        }
        else
        {
            Xloc(0) = -cd/(a1*d2+bc);
            Xloc(1) = bd/(a2*d1+bc);
        }
        return Xloc;
    }
    else
    {
        std::cerr << "local_coord_2D: unknown shape" << std::endl;
        return arma::zeros<arma::vec>(2);
    }
    
}

arma::vec local_coord_2D_along_edge(const arma::vec& X, const arma::mat& Xel, shapes shape, int edge)
{
    if (X.n_elem != 2)
    {
        std::cerr << "local_coord_2D_along_edge: X must have 2 elements" << std::endl;
        return arma::zeros<arma::vec>(2);

    }
    if (Xel.n_rows != 2)
    {
        std::cerr << "local_coord_2D_along_edge: Xel must have 2 rows" << std::endl;
        return arma::zeros<arma::vec>(2);
    }

    double denom1=LOC_EDGE_TOL+1, denom2=LOC_EDGE_TOL+1;
    arma::vec Xloc(2);
    if (shape == TRIANGLE)
    {
        if (Xel.n_cols != 3)
        {
            std::cerr << "local_coord_2D_along_edge: Xel must have 3 columns for tri" << std::endl;
            return arma::zeros<arma::vec>(2);
        }
        switch (edge)
        {
            case 0:
            {
                denom1 = Xel(0,1) - Xel(0,0);
                denom2 = Xel(1,1) - Xel(1,0);
                if (fabs(denom1) > fabs(denom2))
                    Xloc(0) = (X(0) - Xel(0,0))/denom1;                  
                else
                    Xloc(0) = (X(1) - Xel(1,0))/denom2;
                Xloc(1) = 0;  
                break;
            }
            case 1:
            {
                denom1 = Xel(0,1) - Xel(0,2);
                denom2 = Xel(1,1) - Xel(1,2);
                 if (fabs(denom1) > fabs(denom2))
                    Xloc(0) = (X(0) - Xel(0,2))/denom1;                  
                else
                    Xloc(0) = (X(1) - Xel(1,2))/denom2;
                Xloc(1) = 1 - Xloc(0);
                break;
            }
            case 2:
            {
                denom1 = Xel(0,2) - Xel(0,0);
                denom2 = Xel(1,2) - Xel(1,0);
                if (fabs(denom1) > fabs(denom2))
                    Xloc(1) = (X(0) - Xel(0,0))/denom1;                  
                else
                    Xloc(1) = (X(1) - Xel(1,0))/denom2;
                Xloc(0) = 0; 
                break;

            }
            default:
            {
                std::cerr << "local_coord_2D_along_edge: tri: unknown edge" << std::endl;
                break;
            }
        }
        if (fabs(denom1) < LOC_EDGE_TOL && fabs(denom2) < LOC_EDGE_TOL)
        {
            std::cerr << "local_coord_2D_along_edge: tri: both denominators are smaller than tolerance" << std::endl;
        }
    }
    else if (shape == QUADRILATERAL)
    {
        if (Xel.n_cols != 4)
        {
            std::cerr << "local_coord_2D_along_edge: Xel must have 4 columns for quad" << std::endl;
            return arma::zeros<arma::vec>(2);
        }
        switch(edge)
        {
            case 0:
            {
                denom1 = Xel(0,1) - Xel(0,0);
                denom2 = Xel(1,1) - Xel(1,0);
                if (fabs(denom1) > fabs(denom2))
                    Xloc(0) = (2*X(0) - Xel(0,0) - Xel(0,1))/denom1;
                else
                    Xloc(0) = (2*X(1) - Xel(1,0) - Xel(1,1))/denom2;
                Xloc(1) = -1;
                break;
            }
            case 1:
            {
                denom1 = Xel(0,2) - Xel(0,1);
                denom2 = Xel(1,2) - Xel(1,1);
                if (fabs(denom1) > fabs(denom2))
                    Xloc(1) = (2*X(0) - Xel(0,1) - Xel(0,2))/denom1;
                else
                    Xloc(1) = (2*X(1) - Xel(1,1) - Xel(1,2))/denom2;
                Xloc(0) = 1;
                break;
            }
            case 2:
            {
                denom1 = Xel(0,2) - Xel(0,3);
                denom2 = Xel(1,2) - Xel(1,3);
                if (fabs(denom1) > fabs(denom2))
                    Xloc(0) = (2*X(0) - Xel(0,2) - Xel(0,3))/denom1;
                else
                    Xloc(0) = (2*X(1) - Xel(1,2) - Xel(1,3))/denom2;
                Xloc(1) = 1;
                break;
            }
            case 3:
            {
                denom1 = Xel(0,3) - Xel(0,0);
                denom2 = Xel(1,3) - Xel(1,0);
                if (fabs(denom1) > fabs(denom2))
                    Xloc(1) = (2*X(0) - Xel(0,0) - Xel(0,3))/denom1;
                else
                    Xloc(1) = (2*X(1) - Xel(1,0) - Xel(1,3))/denom2;
                Xloc(0) = -1;
                break;
            }
             default:
            {
                std::cerr << "local_coord_2D_along_edge: quad: unknown edge" << std::endl;
                break;
            }
        }
        if (fabs(denom1) < LOC_EDGE_TOL && fabs(denom2) < LOC_EDGE_TOL)
        {
            std::cerr << "local_coord_2D_along_edge: quad: both denominators are smaller than tolerance" << std::endl;
        }
    }

    return Xloc;

}

int edge_number(const arma::uvec& twoNodes, shapes shape)
{
    if (twoNodes.n_elem != 2)
    {
        std::cerr << "edge_number: first argument must have two nodes" << std::endl;
        return 0;
    }
    std::size_t i1 = twoNodes(0), i2 = twoNodes(1);
    if (i1 > i2)
        std::swap(i1,i2);
    
    if (shape == TRIANGLE)
    {
        if (i1 == 0 && i2 == 1)
            return 0;
        else if (i1 == 1 && i2 == 2)
            return 1;
        else if (i1 == 0 && i2 == 2)
            return 2;
        else
        {
            std::cerr << "edge_number: cannot find edge number for tri, returning 0" << std::endl;
            return 0;
        }

    }
    else if (shape == QUADRILATERAL)
    {
        if (i1 == 0 && i2 == 1)
            return 0;
        else if (i1 == 1 && i2 == 2)
            return 1;
        else if (i1 == 2 && i2 == 3)
            return 2;
        else if (i1 == 0 && i2 == 3)
            return 3;
        else
        {
            std::cerr << "edge_number: cannot find edge number for quad, returning 0" << std::endl;
            return 0;
        }

    }
    else
    {
        std::cerr << "edge_number: unknown shape" << std::endl;
        return 0;
    }

}
/*
int edge_number(const arma::uvec& twoNodes, const arma::uvec& refNodes)
{
    if (twoNodes.n_elem != 2)
    {
        std::cerr << "edge_number: first argument must have two nodes" << std::endl;
        return 0;
    }
    if (refNodes.n_elem == 3)
    {
        arma::uvec i1 = arma::find(refNodes == twoNodes(0),1);
        arma::uvec i2 = arma::find(refNodes == twoNodes(1),1);
        if (i1(0) > i2(0))
            i1.swap(i2);

        if (i1(0) == 0 && i2(0) == 1)
            return 0;
        else if (i1(0) == 1 && i2(0) == 2)
            return 1;
        else if (i1(0) == 0 && i2(0) == 2)
            return 2;
        else
        {
            std::cerr << "edge_number: cannot find edge number for tri, returning 0" << std::endl;
            return 0;
        }

    }
    else if (refNodes.n_elem == 4)
    {
        arma::uvec i1 = arma::find(refNodes == twoNodes(0),1);
        arma::uvec i2 = arma::find(refNodes == twoNodes(1),1);
        if (i1(0) > i2(0))
            i1.swap(i2);

        if (i1(0) == 0 && i2(0) == 1)
            return 0;
        else if (i1(0) == 1 && i2(0) == 2)
            return 1;
        else if (i1(0) == 2 && i2(0) == 3)
            return 2;
        else if (i1(0) == 0 && i2(0) == 3)
            return 3;
        else
        {
            std::cerr << "edge_number: cannot find edge number for quad, returning 0" << std::endl;
            return 0;
        }

    }
    else
    {
        std::cerr << "edge_number: 2nd argument must have three or four nodes" << std::endl;
        return 0;
    }

}
*/
}
