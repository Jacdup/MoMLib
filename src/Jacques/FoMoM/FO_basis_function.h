#ifndef basis_function
#define basis_function

#include "matrix_vector_functions.h"
#include "math.h"
#include "triangle_functions.h"


double BF(double Zeta[3],double points[][3], int edge, int sign, double Delta[]){
 // Computes the basis function from the simplex coordinates
//     edge = (edge+2)%3;
    double N; // Normalisation constant
    int ii;
    double t[3][3];
    
    double Area = TriangleArea(points);
    double normal[3];
    
   for (ii=0; ii<3; ii++){
        t[0][ii] = (points[2][ii] - points[1][ii]); // This is just the side length
        t[1][ii] = (points[0][ii] - points[2][ii]);
        t[2][ii] = (points[1][ii] - points[0][ii]);
    }
    double t_len = Distance(points[edge], points[(edge+1)%3]);
//     double t_len = VectorSize(t[(edge)%3]);
    int t_index;
    int l_index;
    
    int first_order = -1;
    double sf = 1;
    double div = 1;
    double sign_sf = sign;
       t_index = (edge+2)%3;
        l_index = (edge+1)%3;
    //VectorCross(t[1],t[2],normal);
    //double ns[3];
    
    if (abs(sign) == 2){ // (if sign is 2 or -2)
        first_order = 1;
        div =0;
        sign_sf = 0.5*sign;
    }
    
    for (ii=0;ii<3;ii++){
        //ns[ii] = 2*Area*normal[ii];
        Delta[ii] = (t_len/(2*Area))*sign_sf*(first_order*(t[l_index][ii]*Zeta[t_index]) + ((t[t_index][ii]*Zeta[l_index])));
    }
    
    return div*((t_len*sign_sf/(Area)));
}

#endif
