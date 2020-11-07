/*=========================================================
 * Name        : mom.c
 * Author      : Robey Beswick
 * Version     : 3.3 (15/07/2019)
 * Copyright   :
 * Description : version 3.3
 *               Is similar to 3.2 but is a lot more optimized.
 *
 *
 *Author: Robey Beswick
 *=======================================================*/
#include "matrix.h"
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "matrix_vector_functions.h"
#include "FO_basis_function.h"
#include "triangle_functions.h"


// return the integration points and weights based on a gaussian rule

void GaussianQuadrature(double points[][3], int numIntPoints, double **returnPoints)
{
    /*
     * Computes quadrature points on a triangle defined by points.
     * Returns points and weight.
     */
    int i, j, ii;
    
    double quadPoints1 [1][4] = {{1.0/3, 1.0/3, 1.0/3, 1}};
    double quadPoints3 [1][4] = {{2.0/3, 1.0/6, 1.0/6, 1.0/3}};
    
    double quadPoints4 [2][4] = {{1.0/3, 1.0/3, 1.0/3, -0.5625},
    {0.6, 0.2, 0.2, 0.520833333333333}};
    double quadPoints6 [2][4] = {{0.816847572980459, 0.091576213509771, 0.091576213509771, 0.109951743655322},
    {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.223381589678011}};
    double quadPoints7 [3][4] = {{1.0/3, 1.0/3, 1.0/3, 0.225},
    {0.797426985353087, 0.101286507323456, 0.101286507323456, 0.125939180544827},
    {0.059715871789770, 0.470142064105115, 0.470142064105115, 0.132394152788506}};
    double quadPoints12 [4][4] = {{0.873821971016996, 0.063089014491502, 0.063089014491502, 0.050844906370207},
    {0.501426509658179, 0.249286745170910, 0.249286745170910, 0.116786275726379},
    {0.636502499121399, 0.310352451033785, 0.053145049844816, 0.082851075618374},
    {0.636502499121399, 0.053145049844816, 0.310352451033785, 0.082851075618374}};
    double quadPoints13 [5][4] = {{1.0/3, 1.0/3, 1.0/3, -0.149570044467682},
    {0.479308067841920, 0.260345966079040, 0.260345966079040, 0.175615257433208},
    {0.869739794195568, 0.065130102902216, 0.065130102902216, 0.053347235608838},
    {0.048690315425316, 0.312865496004874, 0.638444188569810, 0.077113760890257},
    {0.048690315425316, 0.638444188569810, 0.312865496004874, 0.077113760890257}};
    double quadPoints16 [6][4] = {{1.0/3, 1.0/3, 1.0/3, 0.144315607677787},
    {0.081414823414554, 0.459292588292723, 0.459292588292723, 0.095091634267284},
    {0.658861384496478, 0.170569307751761, 0.170569307751761, 0.103217370534718},
    {0.898905543365938, 0.050547228317031, 0.050547228317031, 0.032458497623198},
    {0.008394777409958, 0.263112829634638, 0.728492392955404, 0.027230314174435},
    {0.008394777409958, 0.728492392955404, 0.263112829634638, 0.027230314174435}};
    double quadPoints25 [9][4] = {{0.333333333333333, 0.333333333333333, 0.333333333333333, 0.090817990382754},
    {0.028844733232685, 0.485577633383657, 0.485577633383657, 0.036725957756467},
    {0.781036849029926, 0.109481575485037, 0.109481575485037, 0.045321059435528},
    {0.141707219414880, 0.307939838764121, 0.550352941820999, 0.072757916845420},
    {0.141707219414880, 0.550352941820999, 0.307939838764121, 0.072757916845420},
    {0.025003534762686, 0.246672560639903, 0.728323904597411, 0.028327242531057},
    {0.025003534762686, 0.728323904597411, 0.246672560639903, 0.028327242531057},
    {0.009540815400299, 0.066803251012200, 0.923655933587500, 0.009421666963733},
    {0.009540815400299, 0.923655933587500, 0.066803251012200, 0.009421666963733}};
    
    
    double quadPoints33 [11][4] = {{0.023565220452390, 0.488217389773805, 0.488217389773805, 0.025731066440455},
    {0.120551215411079, 0.439724392294460, 0.439724392294460, 0.043692544538038},
    {0.457579229975768, 0.271210385012116, 0.271210385012116, 0.062858224217885},
    {0.744847708916828, 0.127576145541586, 0.127576145541586, 0.034796112930709},
    {0.957365299093579, 0.021317350453210, 0.021317350453210, 0.006166261051559},
    {0.115343494534698, 0.608943235779788, 0.275713269685514, 0.040371557766381},
    {0.115343494534698, 0.275713269685514, 0.608943235779788, 0.040371557766381},
    {0.022838332222257, 0.695836086787803, 0.281325580989940, 0.022356773202303},
    {0.022838332222257, 0.281325580989940, 0.695836086787803, 0.022356773202303},
    {0.025734050548330, 0.858014033544073, 0.116251915907597, 0.017316231108659},
    {0.025734050548330, 0.116251915907597, 0.858014033544073, 0.017316231108659}};
    
    
    
    double (*qpPointer)[11][4];
    int qpNum;
    
    if (numIntPoints > 16 && numIntPoints < 25)
        numIntPoints = 16;
    
    switch (numIntPoints)
    {
        case 1:
        case 2:
            qpPointer = &quadPoints1;
            qpNum = 1;
            break;
        case 3:
            qpPointer = &quadPoints3;
            qpNum = 1;
            break;
        case 4:
        case 5:
            qpPointer = &quadPoints4;
            qpNum = 2;
            break;
        case 6:
            qpPointer = &quadPoints6;
            qpNum = 2;
            break;
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
            qpPointer = &quadPoints7;
            qpNum = 3;
            break;
        case 12:
            qpPointer = &quadPoints12;
            qpNum = 4;
            break;
        case 13:
        case 14:
        case 15:
            qpPointer = &quadPoints13;
            qpNum = 5;
            break;
        case 16:
            qpPointer = &quadPoints16;
            qpNum = 6;
            break;
        case 25:
            qpPointer = &quadPoints25;
            qpNum = 9;
            break;
        case 33:
        default:
            qpPointer = &quadPoints33;
            qpNum = 11;
            break;
    }
    
    double AB[3] = {points[1][0]-points[0][0],
    points[1][1]-points[0][1],
    points[1][2]-points[0][2]};
    
    double AC[3] = {points[2][0]-points[0][0],
    points[2][1]-points[0][1],
    points[2][2]-points[0][2]};
    
    double Area = pow(AB[1]*AC[2] - AB[2]*AC[1], 2);
    Area += pow(AB[2]*AC[0] - AB[0]*AC[2], 2);
    Area += pow(AB[0]*AC[1] - AB[1]*AC[0], 2);
    Area = 0.5*sqrt(Area);
    
    int place = 0;
    int iter;
    
    for (i=0; i<qpNum; i++)
    {
        if ((*qpPointer)[i][0] != (*qpPointer)[i][1])
            iter = 3;
        else
            iter = 1;
        
        for (j=0; j<iter; j++)
        {
            //  for(ii=0; ii<3; ii++)
            // {
            //returnPoints[place][ii] = (*qpPointer)[i][j]*points[0][ii] + (*qpPointer)[i][(j+1)%3]*points[1][ii] +
            //        (*qpPointer)[i][(j+2)%3]*points[2][ii];
            returnPoints[place][0] = (*qpPointer)[i][j]; //Zeta_1
            returnPoints[place][1] = (*qpPointer)[i][(j+1)%3]; // Zeta_2
            returnPoints[place][2] = (*qpPointer)[i][(j+2)%3]; // Zeta_3
            // }
            returnPoints[place][3] = Area*(*qpPointer)[i][3];
            place++;
        }
    }
}


void PrintMatrix(double **a, int n)
{
    int i,j;
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            printf("%lf ", a[i][j]);
        printf("\n");
    }
}




void FillPlane(double freq, int P, double *points, int T, double *triangles, int N, complex double *Vvec, int *obs_map, int num_obs, int MODE, int order,double theta_inc, double phi_inc, double eta_pol)
{
    //This code choses the default number of integration points per triangle
    int NumOuterIntPoints = 0;
    //defining a proximity by hand is unwise as this is dependant on triangle diameter
    //double prox = 0.33;
    
    int dir[3][3] = {{2,0},{0,1},{1,2}}; // Incident direction iterator
    
    //define the integration domain vector
    int mode_select[6][4] = {
        {3,3,3,6},
        {3,4,6,7},
        {4,6,7,12},
        {6,7,12,16},
        {7,13,16,25}};
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
        
//        FILE *fp;
//       FILE *fp1;
//        fp = fopen("Delta1.txt","w");
//         fp1 = fopen("Ruv1.txt","w");
//           FILE *fp2;
//       FILE *fp3;
//        fp2 = fopen("Delta.txt","w");
//         fp3 = fopen("Ruv.txt","w");
        
        
        double c0 = 299792456.2;
        double k = (2*M_PI*freq)/c0;
        double mu = 1.25663706143592e-06;
        double w = 2*M_PI*freq;
        double eps = 8.85418781761e-12;
        
        
        int i, j, ii, jj, iter;
        int oip;
        
        for (i=0; i< T; i++) //Observation triangle
        {
            
            //create a 3 by 3 matrix of a single triangles points each row is the x y z coordinate
            //pPoints is the outer triangle
            double pPoints[3][3] = {
                {points[P*0+(int)(triangles[T*0 + i]-1)],points[P*1+(int)(triangles[T*0 + i]-1)],points[P*2+(int)(triangles[T*0 + i]-1)]},
                {points[P*0+(int)(triangles[T*1 + i]-1)],points[P*1+(int)(triangles[T*1 + i]-1)],points[P*2+(int)(triangles[T*1 + i]-1)]},
                {points[P*0+(int)(triangles[T*2 + i]-1)],points[P*1+(int)(triangles[T*2 + i]-1)],points[P*2+(int)(triangles[T*2 + i]-1)]}
            };
            //Outer triangle
//         if (order == 0){
//             int pTri[9] = {(int)(triangles[T*0 + i]-1), (int)(triangles[T*1 + i]-1), (int)(triangles[T*2 + i]-1),(int)(triangles[T*3 + i]), (int)(triangles[T*4 + i]),
//                            (int)(triangles[T*5 + i]), (int)(triangles[T*6 + i]), (int)(triangles[T*7 + i]),(int)(triangles[T*8 + i])};
//             int pbasisindex[3] = {pTri[6],pTri[7],pTri[8]};
//         }else{
            int pTri[15] = {(int)(triangles[T*0 + i]-1), (int)(triangles[T*1 + i]-1), (int)(triangles[T*2 + i]-1),(int)(triangles[T*3 + i]), (int)(triangles[T*4 + i]),
            (int)(triangles[T*5 + i]), (int)(triangles[T*6 + i]), (int)(triangles[T*7 + i]),(int)(triangles[T*8 + i]),
            (int)(triangles[T*9 + i]),(int)(triangles[T*10 + i]),(int)(triangles[T*11 + i]),(int)(triangles[T*12 + i]),
            (int)(triangles[T*13 + i]),(int)(triangles[T*14 + i])};
            int pbasisindex[6] = {pTri[6],pTri[7],pTri[8],pTri[12], pTri[13], pTri[14]};
//         }
            
            if ((abs(pTri[10]) != 1) && (abs(pTri[10]) != 2) ){ // Checks if there are extra columns. If not, solve standard RWG.
                order = 0;
            }else{
                order = 1;
//                 mexPrintf("Not correct, no.");
            }
//         int pbasisindex[3] = {pTri[6],pTri[7],pTri[8]};
            
            if (ObserverOnTriangle(pbasisindex, obs_map, order))
            {
                double pArea = TriangleArea(pPoints);
                double **OuterIntPoints = NULL;
                
                //at this point I know my integration domains and I will allocate memory
                //allocate memory for the integration points, do outer and inner here if far interaction
                NumOuterIntPoints = mode_select[4][MODE]; // Choose highest int points
                OuterIntPoints = mxMalloc(sizeof(double *) * NumOuterIntPoints);
                
                for (iter=0; iter < NumOuterIntPoints; iter++)
                    OuterIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                //find the integration points for the outer triangle
                GaussianQuadrature(pPoints, NumOuterIntPoints, OuterIntPoints);
                
                
                //outer integral
                for (oip=0; oip<NumOuterIntPoints; oip++) //outer integral points
                {
                    int pEdgeIndex ;
                    for (ii=0; ii<(3+(3*order)); ii++) //Edges of observation triangle
                    {	// This checks if the edge is a DOF and returns the edges index in the edge list
                        // pEdgeIndex will determine the row of this interaction on the Zmat
                        
                        if (ii<3){
                            pEdgeIndex = pTri[6+(ii+2)%3];
                        }else{
                            pEdgeIndex = pTri[12+(ii+2)%3];
//                             mexPrintf("Not supposed to be in here now");
                        }
                        int p_obs_index = obs_map[pEdgeIndex-1];
                        
                        // below just checks that the edge index is any number other the -1, meaning it will be a DOF
                        // this if statement asks two questions, is it a basis function? and is it an observer?
                        if ((pEdgeIndex + 1) && (p_obs_index + 1))
                        {
                            //mexPrintf("edge = %d\n", pEdgeIndex);
                            
                            // find the length of the outer edge in question
//                     double lp = Distance(pPoints[ii+1], pPoints[(ii+2)%3]);
                            
                            
                            double pRho[3];
                            double ruv[3];
                            // (ii+2)%3 goes 2,0,1
                            ParametricMap(pPoints, OuterIntPoints[oip], ruv);// this is just for plotting purposes
                            //mexPrintf("%f %f %f\n", OuterIntPoints[oip][0], OuterIntPoints[oip][1], OuterIntPoints[oip][2]);
                            
                            int iiSignSelector = 0;
                            if (ii>=3){
                                iiSignSelector = 6;
                            }
                            
                            double N = BF(OuterIntPoints[oip],pPoints, (ii+2)%3, pTri[3+iiSignSelector+(ii+2)%3], pRho);
                            if (isnan(pRho[0])){
                                mexPrintf("error fillPlane\n");
                            }
                            //VectorSubtract(pPoints[(ii+2)%3], ruv, pRho);
                            //MultiplyVectorByConstant(pRho, 3, pTri[3+(ii+2)%3]);
                            // mexPrintf("pRho = %f %f %f", pRho[0],pRho[1],pRho[2]);
                            
// //
//                          if ((pEdgeIndex == 1)|| (pEdgeIndex == 5)){
// //                            if ((i%2)==0){
//                           //   mexPrintf("intPoints = %f %f %f\n pPoints = %f\n", OuterIntPoints[oip][0],OuterIntPoints[oip][1],OuterIntPoints[oip][2], pPoints[0][0]);
//                             // mexPrintf("ruv = %f %f %f\n", pRho[0],pRho[1],pRho[2]);
//                                        fprintf(fp,"%f,%f,%f\n", pRho[0], pRho[1], pRho[2]); //% Delta1
//                                        fprintf(fp1,"%f,%f,%f\n", ruv[0],ruv[1],ruv[2]);
// 
//                                                 }
// //                          }
//                           if ((pEdgeIndex == 10)||(pEdgeIndex == 42)||(pEdgeIndex == 50)){
// //
// //                           //   mexPrintf("intPoints = %f %f %f\n pPoints = %f\n", OuterIntPoints[oip][0],OuterIntPoints[oip][1],OuterIntPoints[oip][2], pPoints[0][0]);
// //                             // mexPrintf("ruv = %f %f %f\n", pRho[0],pRho[1],pRho[2]);
//                                        fprintf(fp2,"%f,%f,%f\n", pRho[0], pRho[1], pRho[2]);
//                                        fprintf(fp3,"%f,%f,%f\n", ruv[0],ruv[1],ruv[2]);
//                            }
// //                             
                            //This is an issue when using matlabs native edges and boundary solver cuz pEdgeIndex
                            // Fill excitation vector
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
                            
                            Vvec[p_obs_index-1] += -(OuterIntPoints[oip][3])*int_eval;
                            
//                          Vvec[p_obs_index-1] += -(OuterIntPoints[oip][3]/(1))*(cexp(I*k*ruv[2])*pRho[0]);
                            //Vvec[p_obs_index-1] += (OuterIntPoints[oip][3]/pArea)*(-lp/2)*(cexp(I*k*OuterIntPoints[oip][2])*pRho[0]);
                            //Vvec[p_obs_index-1] += -(OuterIntPoints[oip][3])*(cexp(I*k*ruv[2])*n*pRho[0]);
//printf("%d\n", pEdgeIndex);
                            
                        }
                    }
                }
                // Free memory
                for (iter=0; iter < NumOuterIntPoints; iter++)
                    mxFree(OuterIntPoints[iter]);
                mxFree(OuterIntPoints);
                
                
            }//this brace is new if statement
        }
//     fclose(fp);
//     fclose(fp1);
//     fclose(fp2);
//     fclose(fp3);
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
    int order = mxGetScalar(prhs[7]);
    double theta_inc = mxGetScalar(prhs[8]);
    double phi_inc = mxGetScalar(prhs[9]);
    double eta_pol = mxGetScalar(prhs[10]);
    mexPrintf("MODE = %d\n", MODE);
    
    Vvec = mxMalloc(sizeof(complex double) * N);
    
    for (i = 0; i < N; i++)
    {
        Vvec[i] = 0.000 + 0.000*I;
    }
    
    /*==========================================================================================*/
    /*==========================================================================================*/
    FillPlane(freq, P, points, T, triangles, N, Vvec, obs_map, num_obs, MODE, order,theta_inc, phi_inc, eta_pol);
    /*==========================================================================================*/
    /*==========================================================================================*/
    
    plhs[0] =  mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    mxSetComplexDoubles(plhs[0], Vvec);
    
    mexPrintf("\nIf it were done when 'tis done, then 'twere well It were done quickly\n");
    
    
}


