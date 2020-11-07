/*=========================================================
 * Name        : qfillPlane.c
 * Author      : JT du Plessis
 * Version     : 
 * Copyright   :
 * Description : 
 *
 *
 *Author: JT du Plessis
 *=======================================================*/
#include "matrix_vector_functions.h"
#include "quadrature.h"
#include "basis_function.h"
#include "quadCentroid.h"
#include "distance.h"

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


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

void FillPlane(double freq, int P, double *points, int T, double *triangles, int N, complex double *Vvec, int *obs_map, int num_obs, int MODE, double theta_inc, double phi_inc, double eta_pol)
{
    //This code choses the default number of integration points per triangle
    int NumOuterIntPoints = 0;
    //defining a proximity by hand is unwise as this is dependant on triangle diameter
    //double prox = 0.33;
    
    int dir[3][3] = {{2,0},{0,1},{1,2}}; // Incident direction iterator
    
    //define the integration domain vector
//     int mode_select[6][4] = {
//         {3,3,3,6},
//         {3,4,6,7},
//         {4,6,7,12},
//         {6,7,12,16},
//         {7,13,16,25}};
            int mode_select[6][4] = {{4,4,9,4},
    {4,9,9,4},
    {9,9,16,9},
    {9,16,25,9},
    {16,25,36,16}};
    //the RAR level is dependant oon the mode as well.
    int RAR_level[4] = {7,7,7,7};
    
    double unitvec_theta[3] = {cos(theta_inc)*cos(phi_inc), cos(theta_inc)*sin(phi_inc), -sin(theta_inc)};
    double unitvec_phi[3]   = { -sin(phi_inc) ,cos(phi_inc) ,0  };
    double unitvec_beta[3]  = {  -sin(theta_inc)*cos(phi_inc),  -sin(theta_inc)*sin(phi_inc) ,-cos(theta_inc)  };
    double unitvec_eta[3];
    int coord;
    for (coord =0;coord<3; coord++){     
        unitvec_eta[coord] = (-cos(eta_pol)*unitvec_theta[coord])+ (sin(eta_pol)*unitvec_phi[coord]);
    }
    
    double c0 = 299792456.2;
    double k = (2*M_PI*freq)/c0;
    double mu = 1.25663706143592e-06;
    double w = 2*M_PI*freq;
    double eps = 8.85418781761e-12;
    
    int i, j, ii, jj, iter;
    int oip;
//      FILE *fp;
//      FILE *fp1;
// 	
//                 fp = fopen("Delta.txt","w");
//                 fp1 = fopen("Ruv.txt","w");
//     
    for (i=0; i< T; i++) //Observation triangle
    {
        
        //create a 3 by 3 matrix of a single triangles points each row is the x y z coordinate
        //pPoints is the outer triangle
        double pPoints[4][3] = {{points[P*0+(int)(triangles[T*0 + i]-1)],points[P*1+(int)(triangles[T*0 + i]-1)],points[P*2+(int)(triangles[T*0 + i]-1)]},
            {points[P*0+(int)(triangles[T*1 + i]-1)],points[P*1+(int)(triangles[T*1 + i]-1)],points[P*2+(int)(triangles[T*1 + i]-1)]},
            {points[P*0+(int)(triangles[T*2 + i]-1)],points[P*1+(int)(triangles[T*2 + i]-1)],points[P*2+(int)(triangles[T*2 + i]-1)]},
            {points[P*0+(int)(triangles[T*3 + i]-1)],points[P*1+(int)(triangles[T*3 + i]-1)],points[P*2+(int)(triangles[T*3 + i]-1)]}};
        //Outer triangle
        int pTri[12] = {(int)(triangles[T*0 + i]-1), (int)(triangles[T*1 + i]-1), (int)(triangles[T*2 + i]-1),(int)(triangles[T*3 + i]-1), (int)(triangles[T*4 + i]),
            (int)(triangles[T*5 + i]), (int)(triangles[T*6 + i]), (int)(triangles[T*7 + i]),(int)(triangles[T*8 + i]),(int)(triangles[T*9 + i]),
            (int)(triangles[T*10 + i]),(int)(triangles[T*11 + i])};
        int pbasisindex[4] = {pTri[8],pTri[9],pTri[10],pTri[11]};

        if (ObserverOnElement(pbasisindex, obs_map))
        {
            double pArea = quadArea(pPoints);
            double **OuterIntPoints = NULL;
            
            //at this point I know my integration domains and I will allocate memory
            //allocate memory for the integration points, do outer and inner here if far interaction
            NumOuterIntPoints = mode_select[4][MODE]; // Choose highest int points
            OuterIntPoints = mxMalloc(sizeof(double *) * NumOuterIntPoints);
            
            for (iter=0; iter < NumOuterIntPoints; iter++)
                OuterIntPoints[iter] = mxMalloc(sizeof(double) * 4);
            //find the integration points for the outer triangle
            GaussianQuadrature(pPoints, NumOuterIntPoints, OuterIntPoints);
            
            for (ii=0; ii<4; ii++) //Edges of observation triangle
            {	// This checks if the edge is a DOF and returns the edges index in the edge list
                // pEdgeIndex will determine the row of this interaction on the Zmat
                
                //int pEdgeIndex = pTri[6+(ii+2)%3];
                int pEdgeIndex = pTri[(8 + (4-ii)%4)];
                int p_obs_index = obs_map[(pTri[8 + (4-ii)%4]-1)];
                //int p_obs_index = obs_map[(pTri[6+(ii+2)%3]-1)];
                
                // below just checks that the edge index is any number other the -1, meaning it will be a DOF
                // this if statement asks two questions, is it a basis function? and is it an observer?
                if ((pEdgeIndex + 1) && (p_obs_index + 1))
                {
                    //outer integral
                    for (oip=0; oip<NumOuterIntPoints; oip++) //outer integral points
                    { 
                        
                        double pRho[4];
                        
                        //Edit 01/10/2019
                        double ruv[3],drdu[3],drdv[3];
                        BilinearInt(pPoints, OuterIntPoints[oip], ruv, drdu,drdv); // This function call is actually unnecessary, only need drdu/drdv

                        double n = BF(OuterIntPoints[oip], drdu, drdv, ii, pTri[4+(4-ii)%4],pRho); // Get Basis Function

                        /* Uncomment this to visualise basis functions */
//                         if ((pEdgeIndex % 2) != 0){
//                         if ( (pEdgeIndex == 2) ||  (pEdgeIndex == 4) || (pEdgeIndex == 6) || (pEdgeIndex == 8) || (pEdgeIndex == 10) || (pEdgeIndex == 12) ){//|| (pEdgeIndex == 10)){// || (pEdgeIndex == 12)){
//                        if ( (pEdgeIndex == 1) ||  (pEdgeIndex == 3) || (pEdgeIndex == 5) || (pEdgeIndex == 7) || (pEdgeIndex == 9) || (pEdgeIndex == 11) || (pEdgeIndex >= 13)){
//                             if ((pEdgeIndex == 2)){
//                             fprintf(fp,"%f,%f,%f\n", pRho[0], pRho[1], pRho[2]);
//                          fprintf(fp1,"%f,%f,%f\n", ruv[0],ruv[1],ruv[2]);
//                          }


                       // For incident direction:
                        // z: OuterIntPoints[oip][2]*pRho[0]
                        // x: OuterIntPoints[oip][0]*pRho[1]
                        // y: OuterIntPoints[oip][1]*pRho[2]
                        // incident_dir is 0 -> z
                        //                 1 -> x
                        //                 2 -> y
                        int coord;
                        double Ruv_eval;
                        complex double Einc_eval[3];
//                         complex double int_eval[3];
                        double temp[3][1] = {{1}, {1}, {1}};
                        Ruv_eval = dotProduct(ruv,unitvec_beta);
                        for (coord = 0; coord < 3; coord++){
                            Einc_eval[coord] = cexp(-I*k*Ruv_eval) * unitvec_eta[coord]*  pRho[coord];
//                             Einc_eval[coord] = Einc_eval[coord]
//                             int_eval[coord] = Einc_eval[coord] *  pRho[coord];
                        }
//                         complex double int_eval  = dotProduct(Einc_eval, temp);
                        complex double int_eval = 0.0 + 0.0j;
                        
                        for (coord = 0;coord < 3;coord++){
                            int_eval += Einc_eval[coord]*temp[coord][1];
                        }
                        
                        Vvec[p_obs_index-1] += -(OuterIntPoints[oip][3])*int_eval*n;
//                         mxFree(Einc_eval);
//                         Vvec[p_obs_index-1] += -(OuterIntPoints[oip][3])*(cexp(I*k*ruv[1])*n*pRho[2]); // This integral is not in reference coordinates
                    }
                }
            }
            // Free memory
            for (iter=0; iter < NumOuterIntPoints; iter++)
                mxFree(OuterIntPoints[iter]);
            mxFree(OuterIntPoints);
            
            
        }//this brace is new if statement
    }
//                 fclose(fp);
//                 fclose(fp1);
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
    int *obs_map =NULL;
    complex double *Vvec = NULL;
    size_t N,P,T,num_obs ;
    int i,j,k;
    int iter,incident_dir, MODE;
    
    /*Populate points array*/
    //get the number of rows in the points array
    P = mxGetM(prhs[0]);

    //set the pointer to recieve points data
    points = mxGetDoubles(prhs[0]);
    
    /*Populate triangles array*/
    T = mxGetM(prhs[1]);

    //mexgetdoubles recieves doubles, this sets them to integers
    triangles =  mxGetDoubles(prhs[1]);
    
    /*Get the number of edges*/
    N = mxGetScalar(prhs[2]);
    
    /*Fetch the frequency*/
    freq = mxGetScalar(prhs[3]);
    
    /*Get the observer and source map*/
    obs_map = mxGetInt32s(prhs[4]);
    num_obs = mxGetScalar(prhs[5]);
//     incident_dir = mxGetScalar(prhs[6]);
    MODE = mxGetScalar(prhs[6]);
    double theta_inc = mxGetScalar(prhs[7]);
    double phi_inc = mxGetScalar(prhs[8]);
    double eta_pol = mxGetScalar(prhs[9]);
//     mexPrintf("MODE = %d\n", MODE);

    Vvec = mxMalloc(sizeof(complex double) * N);

    for (i = 0; i < N; i++)
    {     
        Vvec[i] = 0.000 + 0.000*I;
    }
    
    /*==========================================================================================*/
    /*==========================================================================================*/
    FillPlane(freq, P, points, T, triangles, N, Vvec, obs_map, num_obs, MODE,theta_inc, phi_inc, eta_pol);
    /*==========================================================================================*/
    /*==========================================================================================*/
    
    plhs[0] =  mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    mxSetComplexDoubles(plhs[0], Vvec);
    
}


