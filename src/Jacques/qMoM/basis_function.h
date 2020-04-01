#ifndef basis_function
#define basis_function

#include "matrix_vector_functions.h"
#include "math.h"

double BF(double uv[2], double drdu[3], double drdv[3], int edge, int sign, double Delta[]){
    /* This function computes the basis function on the reference coordinates
     * then multiplies it with the directional vectors lu and lv, to convert
     * it back to global coordinate system. It returns the normalisation constant.
     * See paper by Graglia, Peterson*/
    
    double Zeta;
    double n_cross_l[3];
    int zu, zv;
    double n[3]; // Unit normal to quad
    double J[3]; // Jacobian
    VectorCross(drdu,drdv,J);
    double J_det = sqrt(pow(J[0],2) + pow(J[1],2) +pow(J[2],2));
    n[0] = J[0]/J_det;
    n[1] = J[1]/J_det;
    n[2] = J[2]/J_det;
    
   switch(edge){
        case 0:
            Zeta = (1-uv[0]); // Zeta_3 * l_4
            zu = 1;
            zv = 0;
            VectorCross(n, drdv, n_cross_l); 
            break;
        case 1:
            Zeta = (1-uv[1]); // Zeta_4 * l_1
            zu = 0;
            zv = 1;
            VectorCross(n, drdu, n_cross_l);
            break;
        case 2:
            Zeta = uv[0]; // Zeta_1 * l_2
            zu = -1;
            zv = 0;
            VectorCross(n, drdv, n_cross_l);
            break;
        case 3:
            Zeta = uv[1]; // Zeta_2 * l_3
            zu = 0;
            zv = -1;
            VectorCross(n, drdu, n_cross_l);
            break;
    }

   // Normalisation constant
    
    double N = VectorSize(n_cross_l); // If reference limits are -1<u,v<1 then divide this by 2
    
    int coord;
    for (coord = 0; coord<3; coord++)
        Delta[coord] = sign * Zeta * ((zv *drdv[coord]) + (zu*drdu[coord]));
    
    return ((1/J_det)*N);
}

#endif
