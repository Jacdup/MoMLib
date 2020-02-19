/*=========================================================
 * Name        : qfarfield.c
 * Author      : JT du Plessis, feat Robey Beswick
 * Version     : 
 * Copyright   :
 * Description : 
 *
 *
 *Author: JT du Plessis. Original RWG by Robey Beswick
 *=======================================================*/
#include "matrix.h"
#include "matrix_vector_functions.h"
#include "quadrature.h"
#include "basis_function.h"
#include "quadCentroid.h"
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
// #include <time.h>


double quadArea(double points[4][3])
{
    // This is the same formula for TriangleArea, but done twice
    // First triangle area
    double A1, A2;
    double AB[3], AC[3], AD[3], cp[3], cd[3];
    VectorSubtract(points[0], points[1], AB);
    VectorSubtract(points[0], points[2], AC);
    VectorCross(AB, AC, cp);
    A1 = 0.5*VectorSize(cp);
    
    // Second triangle area
    VectorSubtract(points[0], points[3], AD);
    VectorCross(AD, AC, cd);
    A2 = 0.5*VectorSize(cd);
    
    return A1 +  A2;
}

double TriangleArea(double points[][3])
{
    double AB[3], AC[3], cp[3];
    VectorSubtract(points[0], points[1], AB);
    VectorSubtract(points[0], points[2], AC);
    VectorCross(AB, AC, cp);
    return 0.5*VectorSize(cp);
}

void BilinearInt(double points[][3], double uv[], double ruv[3], double drdu[3], double drdv[3]){
    
    double b10[3],b11[3],b00[3], b01[3];
    int coord;
    
    
    for (coord=0; coord<3; coord++){
        b00[coord] =  points[0][coord];
        b10[coord] = -points[0][coord] + points[1][coord];
        b01[coord] = -points[0][coord] + points[3][coord] ;
        b11[coord] = (points[0][coord] - points[1][coord]- points[3][coord] + points[2][coord]);
        drdv[coord] = (b01[coord] + b11[coord]*uv[0]);
        drdu[coord] = (b10[coord] + b11[coord]*uv[1]);//b10*u
        ruv[coord] =   (b00[coord] + b01[coord]*uv[1] + b10[coord]*uv[0] + b11[coord]*uv[1]*uv[0]);
    }
    
}

//conversion routines
//-------------------------NEW ROUTINES--------------------------------
void CartesianToSphericalCoords(double cart[3], double sphere[3])
{	//tested and correct
    double r = VectorSize(cart);
    //mexPrintf("r = %f\n", r);
    double theta = acos(cart[2]/r);
    //mexPrintf("theta = %f\n", theta);
    double phi = atan(cart[1]/cart[0]);
    //mexPrintf("phi = %f\n", phi);
    
    sphere[0] = r;
    sphere[1] = theta;
    sphere[2] = phi;
    
    
}


void SphericalToCartesian(double cart[3], double sphere[3])
{	 //tested and correct
    cart[0] = sphere[0]*sin(sphere[1])*cos(sphere[2]);
    //mexPrintf("x = %f\n", cart[0]);
    cart[1] = sphere[0]*sin(sphere[1])*sin(sphere[2]);
    //mexPrintf("y = %f\n", cart[1]);
    cart[2] = sphere[0]*cos(sphere[1]);
    //mexPrintf("phi = %f\n", cart[2]);
    
    
}


void VectorCartesianToSpherical(double r_hat_cart[3], complex double cart_components[3] ,complex double spherical_components[3])
{
    double r = VectorSize(r_hat_cart);
    //mexPrintf("r = %f\n", r);
    double theta = acos(r_hat_cart[2]/r);
    //mexPrintf("theta = %f\n", theta);
    double phi = atan(r_hat_cart[1]/r_hat_cart[0]);
    //mexPrintf("phi = %f\n", phi);
    //
    spherical_components[0] = (sin(theta)*cos(phi)*cart_components[0]) + (sin(theta)*sin(phi)*cart_components[1]) + (cos(theta)*cart_components[2]);
    //mexPrintf("sphere[0] = %f\n", sphere[0]);
    spherical_components[1] = (cos(theta)*cos(phi)*cart_components[0]) + (cos(theta)*sin(phi)*cart_components[1]) - (sin(theta))*cart_components[2];
    //mexPrintf("sphere[1] = %f\n", sphere[1]);
    spherical_components[2] = (-sin(phi))*cart_components[0] + (cos(phi))*cart_components[1] + (0)*cart_components[2];
    //mexPrintf("sphere[2] = %f\n", sphere[2]);
    
}

void setupEpoints(char *input, int amnt_E_points, double *theta_phi_Etheta_Ephi)
{
    //this function populates E_points_spherical with the theta and phi components for farpoints
    //r is always 1
    //Do XY
    int i = 0;
    if (strcmp(input, "XY") == 0)
    {
        //mexPrintf("aweh\n");
        //theta always 90 degrees
        //r always 1
        //increment phi by 2 degrees every time
        //-180 to 180
        for (i = 0; i < amnt_E_points; i++)
        {
            //populate phi
            //remember linear indexing
            //[Col_size*col_index + row_index]
            //theta is col 1
            theta_phi_Etheta_Ephi[amnt_E_points*0 + i] = M_PI/2;
            theta_phi_Etheta_Ephi[amnt_E_points*1 + i] = (360/amnt_E_points)*i*(M_PI/180);
            theta_phi_Etheta_Ephi[amnt_E_points*2 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*3 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*4 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*5 + i] = 0.0;
        }
    }else if (strcmp(input, "XZ") == 0)
    {
        for (i = 0; i < amnt_E_points; i++)
        {
            //populate phi
            //remember linear indexing
            //[Col_size*col_index + row_index]
            //theta is col 1
            theta_phi_Etheta_Ephi[amnt_E_points*0 + i] = (360/amnt_E_points)*i*(M_PI/180);
            theta_phi_Etheta_Ephi[amnt_E_points*1 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*2 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*3 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*4 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*5 + i] = 0.0;
        }
        
    }else if (strcmp(input, "YZ") == 0)
    {
        for (i = 0; i < amnt_E_points; i++)
        {
            //populate phi
            //remember linear indexing
            //[Col_size*col_index + row_index]
            //theta is col 1
            theta_phi_Etheta_Ephi[amnt_E_points*0 + i] = (360/amnt_E_points)*i*(M_PI/180);
            theta_phi_Etheta_Ephi[amnt_E_points*1 + i] = M_PI/2;
            theta_phi_Etheta_Ephi[amnt_E_points*2 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*3 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*4 + i] = 0.0;
            theta_phi_Etheta_Ephi[amnt_E_points*5 + i] = 0.0;
        }
        
    }
}

double Find_pRho_phi(double vec[3], double obsPoint[3])
{
    /*recieves the cartesian pRho and observer point
     * it then returns a scalar value that is the component of pRho in the phi direction*/
    double val =  -vec[0]*sin(obsPoint[2]) + vec[1]*cos(obsPoint[2]);
    return val;
}
double Find_pRho_theta(double vec[3], double obsPoint[3])
{
    /*recieves the cartesian pRho and observer point
     * it then returns a scalar value that is the component of pRho in the theta direction*/
    double val =  vec[0]*cos(obsPoint[1])*cos(obsPoint[2]) + vec[1]*cos(obsPoint[1])*sin(obsPoint[2]) - vec[2]*sin(obsPoint[1]);
    return val;
}

//farfield canculator
void farfield(double freq, int P, int T, double *points, double  *triangles, int amnt_E_points, complex double *I_M, double *theta_phi_Etheta_Ephi)
{
    //define constants
    //double freq = 300000000.0;
    double c0 = 299792456.2;
    double k = (2*M_PI*freq)/c0;
    double mu = 1.25663706143592e-06;
    double w = 2*M_PI*freq;
    double eps = 8.85418781761e-12;
    
    int i, ii, ip, iter;
    int e_iter;
    //can choose accuracy of E vector
    int NumIntPoints = 16;
    double **IntPoints = NULL;
    
    //make a loop for the E field points
    for (e_iter = 0; e_iter < amnt_E_points; e_iter++)//
    {	//r_hat is specific to each observer point
        double r_hat[3] = {1.0,theta_phi_Etheta_Ephi[amnt_E_points*0 + e_iter],theta_phi_Etheta_Ephi[amnt_E_points*1 + e_iter]}; //in this case the same for simplicity
        //mexPrintf("r_hat %d = %f %f %f\n",e_iter, r_hat[0], r_hat[1], r_hat[2] );
        //r_hat will be defined in terms of spherical coordinates so coonvert to cartesian for dot product
        
        //define a total E, the E field due to all the triangles at a SINGLE R
        complex double E_theta_R = {0.0+0.0*I};
        complex double E_phi_R = {0.0+0.0*I};
        
        for (i=0; i < T; i++) //do triangle by triangle again
        {
            
            double pPoints[4][3] = {{points[P*0+(int)(triangles[T*0 + i]-1)],points[P*1+(int)(triangles[T*0 + i]-1)],points[P*2+(int)(triangles[T*0 + i]-1)]},
            {points[P*0+(int)(triangles[T*1 + i]-1)],points[P*1+(int)(triangles[T*1 + i]-1)],points[P*2+(int)(triangles[T*1 + i]-1)]},
            {points[P*0+(int)(triangles[T*2 + i]-1)],points[P*1+(int)(triangles[T*2 + i]-1)],points[P*2+(int)(triangles[T*2 + i]-1)]},
            {points[P*0+(int)(triangles[T*3 + i]-1)],points[P*1+(int)(triangles[T*3 + i]-1)],points[P*2+(int)(triangles[T*3 + i]-1)]}};
            //mexPrintf("%f %f %f\n%f %f %f\n%f %f %f\n ",pPoints[0][0],pPoints[0][1],pPoints[0][2], pPoints[1][0],pPoints[1][1],pPoints[1][2],pPoints[2][0],pPoints[2][1],pPoints[2][2]);
            //allocate the memory for the integration points
            IntPoints = mxMalloc(sizeof(double *) * NumIntPoints);
            
            for (iter=0; iter < NumIntPoints; iter++)
                IntPoints[iter] = mxMalloc(sizeof(double) * 4);
            //get and return the gaussian quadrature points for the triangle
            GaussianQuadrature(pPoints, NumIntPoints, IntPoints);

            //Outer triangle
            int pTri[12] = {(int)(triangles[T*0 + i]-1), (int)(triangles[T*1 + i]-1), (int)(triangles[T*2 + i]-1),(int)(triangles[T*3 + i]), (int)(triangles[T*4 + i]),
            (int)(triangles[T*5 + i]), (int)(triangles[T*6 + i]), (int)(triangles[T*7 + i]),(int)(triangles[T*8 + i]),(int)(triangles[T*9 + i]),
            (int)(triangles[T*10 + i]),(int)(triangles[T*11 + i])};
            //mexPrintf("pTri = %d %d %d %d %d %d %d %d %d\n",pTri[0],pTri[1],pTri[2],pTri[3],pTri[4],pTri[5],pTri[6],pTri[7],pTri[8]);
            
            double pArea = quadArea(pPoints);
            
            //new method
            complex double E_per_tri_theta = {0.0+0.0*I};
            complex double E_per_tri_phi = {0.0+0.0*I};
            
            //now I have a triangle,pick an edge and see if the edge is a basis function edge
            for (ii=0; ii<4; ii++) //Edges of observation triangle
            {
                
                //The basis function number is the solution in this case
                //int pEdgeIndex = pTri[6+(ii+2)%3];//this will also give the current coefficient value
                int pEdgeIndex = pTri[(8 + (4-ii)%4)];
                if (pEdgeIndex + 1)
                {
                    //if there is a basis function on the edge it will have an effect
                    //I am the champion
                    //retrieve the current coefficient value
                    complex double I_m = I_M[pEdgeIndex-1];
                  
                    //now you will iterate through the inner integration points
                    
                    //new method
                    complex double theta_holder = {0.0+0.0*I};
                    complex double phi_holder = {0.0+0.0*I};
                    
                    for (ip = 0; ip < NumIntPoints; ip++)
                    {
                        //find pRho in terms of cartesian coordinates
                        double pRho_cart[3];

                        double ruv[3],drdu[3],drdv[3];
                        BilinearInt(pPoints, IntPoints[ip], ruv, drdu,drdv); // This function call is actually unnecessary, only need drdu/drdv
                        double n = BF(IntPoints[ip],drdu, drdv, ii,pTri[4+(4-ii)%4], pRho_cart); // Get Basis Function
                        
                        //dot product without conversion
                        double r_dot_r_hat = (ruv[0]*cos(r_hat[2]) + ruv[1]*sin(r_hat[2]))*sin(r_hat[1]) + ruv[2]*cos(r_hat[1]);

                        //now new edit we find the i_theta and I_phi components of the fields first
                        double pRho_theta_comp = Find_pRho_theta(pRho_cart, r_hat);
                        double pRho_phi_comp = Find_pRho_phi(pRho_cart, r_hat);
                        
                        //now I have the theta and phi components of pRho
                        //must be multiplied by weight and e component and added to the holder
                        complex double G = (IntPoints[ip][3])*cexp(I*k*r_dot_r_hat); // Green's function (minus constant)
                        
                        theta_holder = theta_holder + ((-I*w*mu*I_m*n*pRho_theta_comp)*G)/(4*M_PI);
                        phi_holder = phi_holder + ((-I*w*mu*I_m*n*pRho_phi_comp)*G)/(4*M_PI);
                        
                        
                    }
                    //will happen a max of three times, adding the basis function effects on the triangle
                    E_per_tri_theta = E_per_tri_theta + theta_holder;
                    E_per_tri_phi = E_per_tri_phi + phi_holder;
                    
                    
                    
                }//ip
                //once we are out of this if statement holder is now the effect of the one basis function on the actual E field, must get added to all the other triangles
                //we want to sum per triangle so need to add up the effect of the triangle as well
      
            }//ii
            
            //now add the effect of the triangle to the total field at R
            E_theta_R = E_theta_R + E_per_tri_theta;
            E_phi_R = E_phi_R + E_per_tri_phi;
            //place value at position in E_field
            for (iter=0; iter < NumIntPoints; iter++)
                mxFree(IntPoints[iter]);
            mxFree(IntPoints);
   
        }//T
        
        //Etheta(returns real and imaginary)
        //theta_phi_Etheta_Ephi[amnt_E_points*2 + e_iter] = creal(E_theta_R);
        //theta_phi_Etheta_Ephi[amnt_E_points*3 + e_iter] = cimag(E_theta_R);
        //Etheta(returns magnitude and phase)
        theta_phi_Etheta_Ephi[amnt_E_points*2 + e_iter] = cabs(E_theta_R);
        theta_phi_Etheta_Ephi[amnt_E_points*3 + e_iter] = carg(E_theta_R);
        //phi(real and imaginery)
        //theta_phi_Etheta_Ephi[amnt_E_points*4 + e_iter] = creal(E_phi_R);
        //theta_phi_Etheta_Ephi[amnt_E_points*5 + e_iter] = cimag(E_phi_R);
        //phi(returns magnitude and phase)
        theta_phi_Etheta_Ephi[amnt_E_points*4 + e_iter] = cabs(E_phi_R);
        theta_phi_Etheta_Ephi[amnt_E_points*5 + e_iter] = carg(E_phi_R);
        
    }//field point
    //for (e_iter = 0; e_iter < amnt_E_points; e_iter++)
    //mexPrintf("\nE_field_theta [%lf] = %lf j%lf E_field_phi[%lf] = %lf j%lf\n",theta_phi_Etheta_Ephi[amnt_E_points*0 + e_iter]*180/M_PI,theta_phi_Etheta_Ephi[amnt_E_points*2 + e_iter], theta_phi_Etheta_Ephi[amnt_E_points*3 + e_iter], theta_phi_Etheta_Ephi[amnt_E_points*1 + e_iter]*180/M_PI, theta_phi_Etheta_Ephi[amnt_E_points*4 + e_iter], theta_phi_Etheta_Ephi[amnt_E_points*5 + e_iter]);
    
}
/*====================================================================================================
 * ========================================Main function================================================
 * ====================================================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variable definitions */
    double *points = NULL;
    double *triangles = NULL;
    double freq;
    double *farpoint = NULL;
    
    
    complex double *I_M = NULL;
    complex double *E_theta = NULL;
    complex double *E_phi = NULL;
    
    
    double *theta_phi_Etheta_Ephi = NULL;
    //clock stuff
//     clock_t start, end;
//     double cpu_time_used;
    
    size_t N,P,T;
    int i,j,k;
    int iter, MODE;
    int blah;//indexin
    int step;
    int amnt_E_points;
    char *plane;
    /*Populate points array*/
    //get the number of rows in the points array
    P = mxGetM(prhs[0]);
    //mexPrintf("%d Points\n",P);
    //set the pointer to recieve points data
    points = mxGetDoubles(prhs[0]);
    //for(iter = 0; iter < P;iter++)
    // mexPrintf("%lf %lf %lf\n",points[P*0 + iter], points[P*1 + iter], points[P*2 + iter]);
    
    
    /*Populate triangles array*/
    T = mxGetM(prhs[1]);
    // mexPrintf("%d Triangles\n",T);
    //mexgetdoubles recieves doubles, this sets them to integers
    triangles =  mxGetDoubles(prhs[1]);
    // mexPrintf("%lf \n",triangles[0]);
    //for(iter = 0; iter < T;iter++)
    //mexPrintf("%f %f %f\n",triangles[T*0 + iter],triangles[T*1 + iter],triangles[T*2 + iter]);
    
    // i= 1;
    //int pTri[9] = {(int)triangles[T*0 + i], (int)triangles[T*1 + i],(int)triangles[T*2 + i],(int)triangles[T*3 + i], (int)triangles[T*4 + i],
    //(int)triangles[T*5 + i],(int)triangles[T*6 + i],(int) triangles[T*7 + i],(int)triangles[T*8 + i]};
    //mexPrintf("%d %d %d %d %d %d %d %d %d\n",pTri[0],pTri[1],pTri[2],pTri[3],pTri[4],pTri[5],pTri[6],pTri[7],pTri[8]);
    
    /*Get the number of edges*/
    N = mxGetScalar(prhs[2]);
    //mexPrintf("%d \n",N);
    
    /*Fetch the frequency*/
    freq = mxGetScalar(prhs[3]);
    // mexPrintf("%lf \n",freq);
    
    
    /*get farpoint*/
    I_M  =  mxGetComplexDoubles(prhs[4]);
    //mexPrintf("%f j%f\n",creal(I_M[0]), cimag(I_M[0]));
    
    step = mxGetScalar(prhs[5]);
    //mexPrintf("%d\n",step);
    /*-----------------make internal pointer E matrix-------------------------------*/
    //retrieve the selected plane
    //mxGetString(prhs[6], plane, 5);
    plane =  mxArrayToString(prhs[6]);
    //mexPrintf("%s\n", plane);
    
    /*Populate triangles array*/
    
    amnt_E_points = (int)360/step;
    //mexPrintf("total E points = %d\n", amnt_E_points);
    
    //this allocates the memory for the return data
    theta_phi_Etheta_Ephi = mxMalloc(sizeof(double) * 6 * amnt_E_points);//the entire memory allocation of
    setupEpoints(plane, amnt_E_points, theta_phi_Etheta_Ephi);
    //for (i = 0; i < amnt_E_points; i++)
    //mexPrintf("%f %f %f %f %f %f\n", theta_phi_Etheta_Ephi[amnt_E_points*0 + i], theta_phi_Etheta_Ephi[amnt_E_points*1 + i], theta_phi_Etheta_Ephi[amnt_E_points*2 + i], theta_phi_Etheta_Ephi[amnt_E_points*3 + i], theta_phi_Etheta_Ephi[amnt_E_points*4 + i], theta_phi_Etheta_Ephi[amnt_E_points*5 + i]);
    /*==========================================================================================*/
    /*==========================================================================================*/
    farfield(freq, P, T, points, triangles, amnt_E_points, I_M, theta_phi_Etheta_Ephi);
    /*==========================================================================================*/
    /*==========================================================================================*/
    
    //end = clock();
    
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    //mexPrintf("the array fill took %f seconds to execute \n", cpu_time_used);
    
    //you don't have to mxFree anything as the prhs and plhs are mxArray types which matlab
    //handles itself
    plhs[0]= mxCreateDoubleMatrix(amnt_E_points, 6, mxREAL);
    mxSetDoubles(plhs[0], theta_phi_Etheta_Ephi);
    
    //mxFree(Zmat);
    
    mexPrintf("\nIf it were done when 'tis done, then 'twere well It were done quickly\n");
    
    
}


