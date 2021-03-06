/*=========================================================
 * Name        : qmom.c
 * Author      : Jacques T du Plessis
 * Version     :
 * Copyright   :
 * Description :
 *
 *
 *Author: JT du Plessis. Original RWG by Robey Beswick
 *=======================================================*/
#include "matrix_vector_functions.h"
#include "quadrature.h"
#include "distance.h"
#include "basis_function.h"
#include "quadCentroid.h"

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <time.h>


double quadArea(double points[4][3])
{
    // This is the same formula for TriangleArea, but done twice
    // First quadrilateral area
    double A1, A2;
    double AB[3], AC[3], AD[3], cp[3], cd[3];
    VectorSubtract(points[0], points[1], AB);
    VectorSubtract(points[0], points[2], AC);
    VectorCross(AB, AC, cp);
    A1 = 0.5*VectorSize(cp);
    
    // Second quadrilateral area
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



void invBilinear(double points[][3], double ruv[3], double *uv){
    
    /* See:
     * http://www.iquilezles.org/www/articles/ibilinear/ibilinear.htm
     * for algorithm.
     * **NB** This function has a check which only handles quads parallel
     * to zx and yz plane. So if the cylinder is defined along the z axis,
     * we need another condition.
     */
    
    int CoordSelect = 0; // 2 = x,y
    int i = 0; //x
    int j = 1; //y
    double E[3],F[3],H[3],G[3];
    double temp[3], temp2[3];
    double u, v;
    
    VectorSubtract(points[0], ruv, H); // H = X-A
    VectorSubtract(points[0], points[1], E); // E = B-A
    VectorSubtract(points[0], points[3], F); // F = D-A
    VectorSubtract(points[1],points[0], temp); // A-B
    VectorSubtract(points[3], points[2], temp2); //C-D
    
    G[0] = temp[0] + temp2[0];
    G[1] = temp[1] + temp2[1];
    G[2] = temp[2] + temp2[2]; // G = A-B + C-D
    
    
    
    //Very important: this check only handles quads parallel to zx and yz plane
    // So if the cylinder is defined along the z axis, we need another condition
    
//     if (fabs(E[0])< EPS && E[1] != 0){ // Quad parallel to xy
//         i = 1;//y
//         j = 0;//x
//     }else if (fabs(F[1])< EPS && F[0] != 0)
//     {
// //         mexPrintf("in here\n");
//         i = 0;
//         j = 2;
//     }else if (fabs(H[1]) < EPS && H[0] != 0)
//     {
//         i = 2;//z
//        j = 0;//x
//     }else{
//         i = 2;//z
//         j = 1;//y
//     }
       if (fabs(E[0])< EPS && E[1] != 0 && F[0] != 0){ // Quad parallel to xy
    	// printf("E = %f %f %f, F = %f %f %f\n", E[0], E[1], E[2],F[0], F[1], F[2]);
    	 i = 1;//y
         j = 0;//x
     }else if (fabs(E[0])< EPS && E[1] != 0 && fabs(F[2]) != 0){
    	 //printf("E2\n");
    	 i = 1;
    	 j = 2;
     }
     else if (fabs(F[1])< EPS && F[0] != 0)
     {
    	// printf("F\n");
 //         mexPrintf("in here\n");
         i = 2;
         j = 0;
     }else if (fabs(H[1]) < EPS && H[0] != 0)
     {
    	 //printf("H\n");
         i = 0;//z
        j = 2;//x
     }else if (fabs(F[0]) < EPS && fabs(F[1]) < EPS) 
     {
    	 i = 1;
    	 j = 2;
    }else{
    	 //printf("Else\n");
         i = 0;//z
         j = 1;//y
     }
    
    
    
    
    
    double k2 = Wedge2D(G,F,i,j);
    double k1 = Wedge2D(E,F,i,j) + Wedge2D(H,G,i,j) ;
    double k0 = Wedge2D(H,E,i,j);
    
    double k2u = Wedge2D(E,G,i,j);
    double k1u = Wedge2D(E,F,i,j) + Wedge2D(G,H,i,j);
    double k0u = Wedge2D(H,F,i,j);
    
    double v1, u1, v2, u2;
    if (fabs(k2) < EPS){
        v1 = -k0/k1;
        u1 = (H[i] - F[i]*v1)/(E[i] + G[i]*v1);
    }
    else if (fabs(k2u) < EPS){
        u1 = k0u /k1u;
        v1 = (H[j] - E[j]*u1)/(F[j] + G[j]*u1);
    }
    else{
        
        double w = k1*k1 - 4.0*k0*k2; // Discriminant
        
        if (w<0.0) {
            u = -1;
            v = -1;
            printf("ABORT!!");
        }
        w = sqrt(w);
        
        v1 = (-k1 -w)/(2.0*k2);
        u1 = (-k1u -w)/(2.0*k2u);
        
        v2 = (-k1 +w)/(2.0*k2);
        u2 = (-k1u + w)/(2.0*k2u);
        
    }
    
    bool b1 = ( v1>0.0 && v1<1.0 && u1>0.0 && u1<1.0 );
    bool b2 = ( v2>0.0 && v2<1.0 && u2>0.0 && u2<1.0);
    
    if (  b1 && !b2 ) {u = u1; v = v1;}
    else if ( !b1 &&  b2 ) {u = u2; v = v2;}
    else {u = u1; v = v1;}
    
    if (fabs(u2)< EPS && fabs(v2) < EPS){ // This is sort of a check when the observation point is outside of quad (which is a valid scenario)
        u = u1;
        v = v1;
    }else if (fabs(u1)< EPS && fabs(v1) < EPS){
        u = u2;
        v = v2;
    }else if (fabs(u1) < fabs(u2) && fabs(v1) < fabs(v2)){
        u = u1;
        v = v1;
    }
    else{
//         u = u*2 - 1;
//         v = v*2 - 1;
         u = u2;
    	 v = v2;
//          mexPrintf("u = %f\n", u);
//         u = 1;
//         v = 1;
    }
    
    uv[0] = u;
    uv[1] = v;
    
    // Get the limits to -1 < u,v < 1
    //uv[0] = uv[0]*2 -1;
    //uv[1] = uv[1]*2 -1;
    
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




void MoM(double freq, int P, double *points, int T, double *triangles, int N, complex double *Zmat,int *obs_map, int *src_map, int num_obs, int num_src, int MODE)
{
    
    /*Variables:
     *P         =>  #new_points
     *points    =>  pointer to new_points
     *T         =>  #new_triangles
     *triangles =>  pointer to new_triangles
     *N         =>  new_N (number of edges)
     *Zmat      =>  memory allocated Z-matrix
     */
    
    //This code choses the default number of integration points per quadrilateral and then deviates from it with checks
    int NumOuterIntPoints = 0;
    int NumInnerIntPoints = 0;
    
    //define the integration domain vector
//     int mode_select[6][4] = {{3,3,6,3},
//     {4,6,7,3},
//     {6,7,12,4},
//     {7,12,16,6},
//     {13,16,25,7}};
    int mode_select[6][4] = {{4,4,9,4},
    {4,9,9,4},
    {9,9,16,9},
    {9,16,25,9},
    {16,25,36,16}};
    
    //the RAR level is dependant oon the mode as well.
    int RAR_level[4] = {7,7,7,7};
    
    
    double c0 = 299792456.2;
    double k = (2*M_PI*freq)/c0;
    double mu = 1.25663706143592e-06;
    double w = 2*M_PI*freq;
    double eps = 8.85418781761e-12;
    
    int i, j, ii, jj, iter;
    int oip, iip;
    
    //Used to count occurences of integration types
    int count_16_RAR = 0; //used to count occurances of cancellation quadrature
    int count_4_4 = 0;
    int count_6_6 = 0;
    int count_7_7 = 0;
    int count_12_12 = 0;
    int count_16_16 = 0;
    int count = 0;
    
    double PhiTot = 0;
    
    //clock stuff
    clock_t StartInner, EndInner;
    double Inner_loop_time_used;
    FILE *fp_ruv;
    fp_ruv = fopen("ruv_bilinear.txt","w");
    
    double prog = 0.0; //used to monitor progress
    
    for (i=0; i< T; i++) //Observation quadrilateral
    {
        
        //create a 3 by 3 matrix of a single triangles points each row is the x y z coordinate
        //pPoints is my outer quadrilateral
        //new reference of pPoints
        double pPoints[4][3] = {{points[P*0+(int)(triangles[T*0 + i]-1)],points[P*1+(int)(triangles[T*0 + i]-1)],points[P*2+(int)(triangles[T*0 + i]-1)]},
        {points[P*0+(int)(triangles[T*1 + i]-1)],points[P*1+(int)(triangles[T*1 + i]-1)],points[P*2+(int)(triangles[T*1 + i]-1)]},
        {points[P*0+(int)(triangles[T*2 + i]-1)],points[P*1+(int)(triangles[T*2 + i]-1)],points[P*2+(int)(triangles[T*2 + i]-1)]},
        {points[P*0+(int)(triangles[T*3 + i]-1)],points[P*1+(int)(triangles[T*3 + i]-1)],points[P*2+(int)(triangles[T*3 + i]-1)]}};
        
        //Outer quadrilateral
        // This is the quad nodes as well as directional and DOF data
        int pQuad[12] = {(int)(triangles[T*0 + i]-1), (int)(triangles[T*1 + i]-1), (int)(triangles[T*2 + i]-1),(int)(triangles[T*3 + i]-1), (int)(triangles[T*4 + i]),
        (int)(triangles[T*5 + i]), (int)(triangles[T*6 + i]), (int)(triangles[T*7 + i]),(int)(triangles[T*8 + i]),(int)(triangles[T*9 + i]),
        (int)(triangles[T*10 + i]),(int)(triangles[T*11 + i])};
        
        //only need to do once
        int pbasisindex[4] = {pQuad[8],pQuad[9],pQuad[10],pQuad[11]};
        //VERSION 3.3 ADDITION
        if (ObserverOnElement(pbasisindex, obs_map))
        {
            
            double pArea = quadArea(pPoints);
            
            //double **OuterIntPoints = NULL;
            //StartInner = clock();
            for (j=0; j<T; j++) //Testing quadrilateral
            {
                // qPoints is my inner quadrilateral where I will possibly need to use RAR to integrate
                double qPoints[4][3] = {{points[P*0+(int)(triangles[T*0 + j]-1)],points[P*1+(int)(triangles[T*0 + j]-1)],points[P*2+(int)(triangles[T*0 + j]-1)]},
                {points[P*0+(int)(triangles[T*1 + j]-1)],points[P*1+(int)(triangles[T*1 + j]-1)],points[P*2+(int)(triangles[T*1 + j]-1)]},
                {points[P*0+(int)(triangles[T*2 + j]-1)],points[P*1+(int)(triangles[T*2 + j]-1)],points[P*2+(int)(triangles[T*2 + j]-1)]},
                {points[P*0+(int)(triangles[T*3 + j]-1)],points[P*1+(int)(triangles[T*3 + j]-1)],points[P*2+(int)(triangles[T*3 + j]-1)]}};
                
                //inner quadrilateral
                int qQuad[12] = {(int)(triangles[T*0 + j]-1), (int)(triangles[T*1 + j]-1), (int)(triangles[T*2 + j]-1),(int)(triangles[T*3 + j]-1), (int)(triangles[T*4 + j]),
                (int)(triangles[T*5 + j]), (int)(triangles[T*6 + j]), (int)(triangles[T*7 + j]),(int)(triangles[T*8 + j]),(int)(triangles[T*9 + j]),
                (int)(triangles[T*10 + j]),(int)(triangles[T*11 + j])};
                
                //below checks if any single edge in the quadrilateral has a basis function defined over it
                int qbasisindex[4] = {qQuad[8],qQuad[9],qQuad[10],qQuad[11]};
                if (ObserverOnElement(qbasisindex, src_map))
                {
                    double **OuterIntPoints = NULL;
                    // Find the area and centrepoint of the quad
                    double qArea = quadArea(qPoints);
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //-----------------------------------------------------First round of checks--------------------------------------------------------
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //these are local variables, they change per quadrilateral pair, will only be used
                    //if completely necessery
                    double d_conservative = 0.0;
                    double max_diameter   = 0.0;
                    double pCentroid[3];
                    double qCentroid[3];
                    double d_centrepoints = 0.0;
                    double pDiameter = 0.0;
                    double qDiameter = 0.0;
                    double d_9 = 0.0;
                    double d_36 = 0.0;
                    int further_check_required = 0;
                    int newNumInt = 0;
                    
                    QuadCentre(pPoints,pCentroid);
                    QuadCentre(qPoints,qCentroid);
                    d_centrepoints = Distance(pCentroid,qCentroid);
                    
                    //find the quad diameters and the maximum of them
                    pDiameter = QuadDiameter(pPoints);
                    qDiameter = QuadDiameter(qPoints);
                    
                    
                    if (pDiameter > qDiameter)
                    {
                        max_diameter = pDiameter;
                    }else{
                        
                        max_diameter = qDiameter;
                    }
                    //shortest possible distance between triangles given the centroid
                    d_conservative = d_centrepoints - (2.0/3.0)*(pDiameter + qDiameter);
                    
                    //d conservative is the shortest possible distance between the triangles given the centroid
                    if (0 > d_conservative)
                    {
                        //if this is the case the triangles are either touching or very close,
                        //it will always be negative if the triangles are touching
                        int DontTouch = ElementsDontTouch(pQuad, qQuad);//can only be used locally
                        if (!DontTouch)
                        {
                            //if we are in here it means the triangles touch
                            //instantly set outer to 16 and further_checks = 1
                            //no inner points set yet, further checks handles that
                            //set mode to level 4
                            NumOuterIntPoints = mode_select[4][MODE];
                            
                            //NumOuterIntPoints = 16;
                            further_check_required = 1;
                            
                        }else
                        {
                            //here we will perform a 9 point check to see how far apart they actually are
                            d_9 = ApproxShortDistanceBetweenElements_9(pPoints,qPoints);
                            
                            if(d_9 > (1.5*max_diameter))
                            {
                                //set to 6 point
                                //set mode to level 1
                                NumOuterIntPoints = mode_select[1][MODE];
                                NumInnerIntPoints = mode_select[1][MODE];
                                count_6_6++;
                            }else if(d_9 > (0.75*max_diameter))
                            {
                                //set to 7 point
                                //set mode to level 2
                                NumOuterIntPoints = mode_select[2][MODE];
                                NumInnerIntPoints = mode_select[2][MODE];
                                count_7_7++;
                            }else
                            {
                                //This means the triangles are very close and we will do strict test
                                //do 36 point check
                                d_36 = ApproxShortDistanceBetweenElements_36(pPoints,qPoints);
                                if(d_36 > 0.375*max_diameter)
                                {
                                    //set both to 12 point
                                    //set mode to level 3
                                    NumOuterIntPoints = mode_select[3][MODE];
                                    NumInnerIntPoints = mode_select[3][MODE];
                                    count_12_12++;
                                }else
                                {
                                    //set outer to 16 and further_checks = 1
                                    //set mode to level 4
                                    NumOuterIntPoints = mode_select[4][MODE];
                                    further_check_required = 1;
                                }
                            }
                        }
                        
                    }else
                    {
                        //d_conservative is positive so continue with checks
                        //now we do the big conservative checks
                        if (d_conservative > (3.0*max_diameter))
                        {
                            //set to 4 point
                            NumOuterIntPoints = mode_select[0][MODE];
                            NumInnerIntPoints = mode_select[0][MODE];
                            count_4_4++;
                        }else if(d_conservative > (1.5*max_diameter))
                        {
                            //set mode to level 1
                            NumOuterIntPoints = mode_select[1][MODE];
                            NumInnerIntPoints = mode_select[1][MODE];
                            count_6_6++;
                        }else if(d_conservative > (0.75*max_diameter))
                        {
                            //set mode to level 2
                            NumOuterIntPoints = mode_select[2][MODE];
                            NumInnerIntPoints = mode_select[2][MODE];
                            count_7_7++;
                        }else
                        {	//lower then 0.75*max_d
                            //we perform the 9 point check
                            d_9 = ApproxShortDistanceBetweenElements_9(pPoints,qPoints);
                            if (d_9 > (1.5*max_diameter))
                            {
                                //set mode to level 1
                                NumOuterIntPoints = mode_select[1][MODE];
                                NumInnerIntPoints = mode_select[1][MODE];
                                count_6_6++;
                            }else if(d_9 > (0.75*max_diameter))
                            {
                                //set mode to level 2
                                NumOuterIntPoints = mode_select[2][MODE];
                                NumInnerIntPoints = mode_select[2][MODE];
                                count_7_7++;
                            }else
                            {
                                //so if less then 0.75 multiplied by the max diameter
                                //more stringent tests applied
                                double d_36 = ApproxShortDistanceBetweenElements_36(pPoints,qPoints);
                                if(d_36 > 0.375*max_diameter)
                                {
                                    //set mode to level 3
                                    NumOuterIntPoints = mode_select[3][MODE];
                                    NumInnerIntPoints = mode_select[3][MODE];
                                    count_12_12++;
                                }else
                                {
                                    //set outer to 16 and further_checks = 1
                                    //set mode to level 4, further checks will chose later mode
                                    NumOuterIntPoints = mode_select[4][MODE];
                                    further_check_required = 1;
                                }
                            }
                        }
                    }
                    //at this point I know my integration domains and I will allocate memory
                    //allocate memory for the integration points, do outer and inner here if far interaction
                    int invert = 1;
                    if (!further_check_required)
                    {
                        //far interaction (further checks aren't required)
                        //mxMalloc outer integration points
                        OuterIntPoints = mxMalloc(sizeof(double *) * NumOuterIntPoints);
                        for (iter=0; iter < NumOuterIntPoints; iter++)
                            OuterIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                        //find the integration points for the outer quadrilateral
                        GaussianQuadrature(pPoints, NumOuterIntPoints, OuterIntPoints);
                        
                        invert = 0;
                    }else
                    {
                        //This is a near interaction a check must now be done per outer integration point.
                        //set the outer integration points to a high value(level 4), I will play with this later
                        
                        OuterIntPoints = mxMalloc(sizeof(double *) * NumOuterIntPoints);
                        for (iter=0; iter < NumOuterIntPoints; iter++)
                            OuterIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                        //find the integration points for the outer quadrilateral
                        GaussianQuadrature(pPoints, NumOuterIntPoints, OuterIntPoints);
                        
                    }
                    //once this round of checks is complete we will either know our symmetrical number of
                    //inner and outer integration points or we will know our outer integration points value
                    //is set to level 4 outer integral value
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //------------------------------------------------First round of checks complete----------------------------------------------------
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //----------------------------------------------------------------------------------------------------------------------------------
                    int pQuadShape = 0;
                    int qQuadShape = 0;
                    if (pQuad[2] == pQuad[3]){ // If degenerate
                        pQuadShape = 1;
                    }
                    if (qQuad[2] == qQuad[3]){
                        qQuadShape = 1;
                    }
                    int RAR_flag = 0;
                    
                    //I only want to assign the outer integration points here
                    
                    //outer integral
                    for (oip=0; oip<NumOuterIntPoints; oip++) //outer integral points
                    {
                    
                        double **InnerIntPoints = NULL;
                        if(further_check_required && (i == j))
                        {
                            
                            //if we come in here outer is already set to 16
                            //16 point outside and RAR inside
                            
                            //we will use radial angular for the inner integral
                            
                            //NumInnerIntPoints is recieved by RAR, max of 10 is set
                            //NumInnerIntPoints is defined by mode
                            NumInnerIntPoints = RAR_level[MODE];
                            newNumInt = 4*(int)pow(NumInnerIntPoints,2);
                            if (newNumInt > 400)
                                newNumInt = 400;
                            
                            //allocate the memory for OuterIntPoints
                            InnerIntPoints = mxMalloc(sizeof(double *)*newNumInt);
                            for (iter=0; iter < newNumInt; iter++)
                                InnerIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                            
                            //so the inner quadrilateral is now integrated using radial angular methods
                            //Note that RAR1S also recieves the outer integration point and original NumInnerIntPoints
                            //17-10-2019 Edit
                            double ruv[3],drdu[3],drdv[3];
                            BilinearInt(pPoints, OuterIntPoints[oip], ruv, drdu,drdv); // Get the x,y,z coords for RAR1S
                            
                           if ((pQuadShape == 1) || (qQuadShape == 1)){
                               // if (i == (T-1)){
                                //  fprintf(fp_ruv,"%f,%f,%f\n", ruv[0], ruv[1], ruv[2]);}
                                  //mexPrintf("node = %f, pQuad = %d %d %d %d\n",pPoints[3][0], pQuad[0], pQuad[1],pQuad[2],pQuad[3]);
                    //      fprintf(fp1,"%f,%f,%f\n", ruv[0],ruv[1],ruv[2]);
                                //    mexPrintf("ruv = %f %f %f\n", ruv[0],ruv[1],ruv[2]);}
                               RAR1STri(qPoints, ruv, NumInnerIntPoints, InnerIntPoints);
                              
                               // mexPrintf("inner points = %f %f %f\n", InnerIntPoints[0][0],InnerIntPoints[0][1],InnerIntPoints[0][2]);
                            }
                            else{
                                RAR1S(qPoints, ruv, NumInnerIntPoints, InnerIntPoints);
//                                   if (i == (T-1)){
//                                 fprintf(fp_ruv,"%f,%f,%f\n", InnerIntPoints[2][0], InnerIntPoints[2][1], InnerIntPoints[2][2]);}
                           }
                            //RAR1S(qPoints, OuterIntPoints[oip], NumInnerIntPoints, InnerIntPoints);
                            count_16_RAR++;
                            RAR_flag = 1;
                            invert = 1;
                            
                        }else if(further_check_required)
                        {
                            //check accurate proximity from oip to inner quadrilateral
                            double proximity = StrictDistance(qPoints, OuterIntPoints[oip]);
                            if(proximity > 0.375*max_diameter)
                            {
                                //set both to 16 point
                                
                                //set inner inner integral to 16 point
                                //mxMalloc inner integration points
                                NumInnerIntPoints = mode_select[4][MODE];
                                newNumInt = NumInnerIntPoints;
                                //The pointer to InnerIntPoints was already set in last checks
                                InnerIntPoints = mxMalloc(sizeof(double *) * NumInnerIntPoints);
                                for (iter=0; iter < NumInnerIntPoints; iter++)
                                    InnerIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                                
                                GaussianQuadrature(qPoints, NumInnerIntPoints, InnerIntPoints);
                                count_16_16++;
                                RAR_flag = 0;
                                invert = 0;
                            }else
                            {
                                //if we come in here outer is already set to level 4
                                
                                //level 4 outside and RAR inside
                                //we will use radial angular for the inner integral
                                
                                //NumInnerIntPoints is recieved by RAR, this can be set higher max of 10 is set
                                NumInnerIntPoints = RAR_level[MODE];
                                newNumInt = 4*(int)pow(NumInnerIntPoints,2);
                                if (newNumInt > 400)
                                    newNumInt = 400;
                                
                                //allocate the memory for InnerIntPoints
                                InnerIntPoints = mxMalloc(sizeof(double *)*newNumInt);
                                for (iter=0; iter < newNumInt; iter++)
                                    InnerIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                                
                                //so the inner quadrilateral is now integrated using radial angular methods
                                //Note that RAR1S also recieves the outer integration point and original NumInnerIntPoints
                                //17-10-2019 Edit
                                double ruv[3],drdu[3],drdv[3];
                                BilinearInt(pPoints, OuterIntPoints[oip], ruv, drdu,drdv);
                                if ((pQuadShape == 1)|| (qQuadShape == 1)){
                                    //    mexPrintf("ruv = %f %f %f\n", ruv[0],ruv[1],ruv[2]);}
                                    RAR1STri(qPoints, ruv, NumInnerIntPoints, InnerIntPoints);
                                    // mexPrintf("success\n");}
                                }
                                else{
                                    RAR1S(qPoints, ruv, NumInnerIntPoints, InnerIntPoints);
                               }
                                // RAR1S(qPoints, ruv, NumInnerIntPoints, InnerIntPoints);
                                //RAR1S(qPoints, OuterIntPoints[oip], NumInnerIntPoints, InnerIntPoints);
                                count_16_RAR++;
                                RAR_flag = 1;
                                invert = 1;
                            }
                            
                            
                        }else if(!further_check_required  && (i != j))
                        {
                            newNumInt = NumInnerIntPoints;
                            InnerIntPoints = mxMalloc(sizeof(double *) * NumInnerIntPoints);
                            for (iter=0; iter < NumInnerIntPoints; iter++)
                                InnerIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                            
                            //find the integration points for the inner quadrilateral
                            GaussianQuadrature(qPoints, NumInnerIntPoints, InnerIntPoints);
                            RAR_flag = 0;
                        }
                        //----------------------------------------------------------------------------------------------------------------------------------
                        //----------------------------------------------------------------------------------------------------------------------------------
                        //------------------------------------------------second round of checks complete---------------------------------------------------
                        //----------------------------------------------------------------------------------------------------------------------------------
                        //----------------------------------------------------------------------------------------------------------------------------------
                        
                        //AT THIS POINT MY OUTER AND INNER INTEGRATION POINTS ARE KNOWN
                        for (ii=0; ii<4; ii++) //Edges of observation quadrilateral
                        {	// This checks if the edge is a DOF and returns the edges index in the edge list
                            // pEdgeIndex will determine the row of this interaction on the Zmat
                            
                            int pEdgeIndex = pQuad[(8 + (4-ii)%4)]; // AD -> AB -> BC -> DC
                            int p_obs_index = obs_map[(pQuad[8 + (4-ii)%4]-1)];
                            // 20/02/2020 Edit
                            
                            
                            
                            // below just checks that the edge index is any number other the -1, meaning it will be a DOF
                            // this if statement asks two questions, is it a basis function? and is it an observer?
                            if ((pEdgeIndex + 1) && (p_obs_index + 1))
                            {
                                
                                //------------Outer Basis Function------------//
                                double pRho[3],OuterRuv[3],drdu_outer[3],drdv_outer[3];
                                BilinearInt(pPoints, OuterIntPoints[oip], OuterRuv, drdu_outer,drdv_outer); // This function call is actually unnecessary, only need drdu/drdv
                                double n_outer = BF(OuterIntPoints[oip], drdu_outer, drdv_outer, ii,pQuad[4+(4-ii)%4],pRho); // Get Basis Function
                                
                                for (jj=0; jj<4; jj++) //Edges of testing quadrilateral
                                {
                                    
                                    //VERSION 3.2.2
                                    int qEdgeIndex = qQuad[(8 + (4-jj)%4)];
                                    int q_src_index = src_map[(qQuad[8 + (4-jj)%4]-1)];
                                    
                                    if ((qEdgeIndex + 1) && (q_src_index + 1))
                                    {
                                        
                                        complex double A[3] = {0.0, 0.0, 0.0};
                                        complex double Phi = 0+0*I;
                                        
                                        int index;
                                        
                                        //inner integral
                                        for (iip=0; iip < newNumInt; iip++)
                                        {
                                            double qRho[3];
                                            
                                            //Edit 17/10/2019
                                            double print = 0;
                                            double uv[2];
                                            
                                            /* This just checks if the inner integral was computed with RAR1S or Gaussian Quadrature,
                                             * since RAR1S returns R(u,v) and Gauss returns u,v */
                                            if (RAR_flag){
                                                invBilinear(qPoints, InnerIntPoints[iip], uv); //Get uv, then do the magic with uv, drdu,drdv
                                                 if (isnan(uv[0]) ){
//                                                      mexPrintf("Problem inv\n");
                                                     mexPrintf("qPoints = %f,%f,%f %f,%f,%f %f,%f,%f, %f,%f,%f, Inner = %f,%f,%f\n",qPoints[0][0],qPoints[0][1],qPoints[0][2],qPoints[1][0],qPoints[1][1],qPoints[1][2],qPoints[2][0],qPoints[2][1],qPoints[2][2],qPoints[3][0],qPoints[3][1],qPoints[3][2], InnerIntPoints[iip][0], InnerIntPoints[iip][1],InnerIntPoints[iip][2]);
                                                     //invBilinear(qPoints, InnerIntPoints[iip], uv);
                                                 }
                                                
//                                                 if (i == (T-1)){
//                                                      if (iip == newNumInt-1){
//                                                          fprintf(fp_ruv,"%f,%f,%f\n", InnerIntPoints[iip][0], InnerIntPoints[iip][1],InnerIntPoints[iip][2]);}}
                                            }
                                            else{
                                                uv[0] = InnerIntPoints[iip][0];
                                                uv[1] = InnerIntPoints[iip][1];
                                            }
                                            
                                            //Edit 21/10/2019
                                            //Compute basis function for inner integration
                                            double ruv[3],drdu[3],drdv[3];
                                            BilinearInt(qPoints, uv, ruv, drdu,drdv);
                                            double n = BF(uv, drdu, drdv, jj,qQuad[4+(4-jj)%4], qRho); // Get Basis Function
                                            
                                            double Rm = Distance(OuterRuv, ruv);
                                            complex double temp = ((InnerIntPoints[iip][3])/(Rm))*cexp(-I*k*Rm); // Green's Function
                                            if (isnan(cimag(temp))){
//                                                 mexPrintf("outer, ruv = %f %f %f, %f %f %f\n", OuterRuv[0],OuterRuv[1],OuterRuv[2],ruv[0],ruv[1],ruv[2]);
                                            }
                                            
//                                             if (isnan(cimag(temp))){
//                                                 mexPrintf("imag temp\n");
//                                             }
                                            
                                            A[0] += ((qRho[0])*n*I*mu*w*temp)/(4*M_PI); //x
                                            A[1] += ((qRho[1])*n*I*mu*w*temp)/(4*M_PI); //y
                                            A[2] += ((qRho[2])*n*I*mu*w*temp)/(4*M_PI); //z
                                            
                                            Phi += (n*qQuad[4+(4-jj)%4]*temp)/(4*I*M_PI*w*eps);
                                        }
                                        complex double AdotpRho = A[0]*(pRho[0]) + A[1]*(pRho[1]) + A[2]*(pRho[2]);
                                        
                                        Zmat[(num_obs*(q_src_index-1)) + (p_obs_index-1)] += (n_outer*(OuterIntPoints[oip][3]))*(AdotpRho + (pQuad[4+(4-ii)%4]*Phi));
                                    }
                                }

                                // if (pEdgeIndex ==6){
                                //   PhiTot = PhiTot + (OuterIntPoints[oip][3]*pQuad[4+(4-ii)%4]);
                                //    mexPrintf("Grad1,Grad2 = %f,%f, Phi Total = %f\n",grad1_outer_abs, grad2_outer_abs,PhiTot);}
                            }//oip iterator
                            
                        }
                        for (iter=0; iter < newNumInt; iter++)
                                    mxFree(InnerIntPoints[iter]);
                                mxFree(InnerIntPoints);
                        
                    }
                    for (iter=0; iter < NumOuterIntPoints; iter++)
                        mxFree(OuterIntPoints[iter]);
                    mxFree(OuterIntPoints);
                }//the brace is new inner if statement
            }
            
            
            //EndInner = clock();
            //Inner_loop_time_used = ((double) (EndInner - StartInner)) / CLOCKS_PER_SEC;
            //mexPrintf("inner loop took %f seconds to execute \n", Inner_loop_time_used);
            
        }//this brace is new if statement
    }
    //prints the occurences of how many times each domain was used
    fclose(fp_ruv);
    mexPrintf("\n%d %d: %d\n%d %d: %d\n%d %d: %d\n%d %d: %d\n%d %d: %d\n%d RAR: %d\n",mode_select[0][MODE],mode_select[0][MODE] ,count_4_4, mode_select[1][MODE],mode_select[1][MODE],count_6_6, mode_select[2][MODE],mode_select[2][MODE],count_7_7, mode_select[3][MODE],mode_select[3][MODE],count_12_12, mode_select[4][MODE],mode_select[4][MODE], count_16_16,mode_select[4][MODE], count_16_RAR);
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
    int *src_map = NULL;
    //complex double *Zarray = NULL;
    complex double *Zmat = NULL;
    
    //clock stuff
    clock_t start, end;
    double cpu_time_used;
    
    size_t N,P,T,num_obs, num_src;
    int i,j,k;
    int iter, MODE;
    int blah;//indexin
    
    /*Populate points array*/
    //get the number of rows in the points array
    P = mxGetM(prhs[0]);
    // mexPrintf("%d \n",P);
    //set the pointer to recieve points data
    points = mxGetDoubles(prhs[0]);
    // for(iter = 0; iter < P;iter++)
    // mexPrintf("%lf %lf %lf\n",points[P*0 + iter], points[P*1 + iter], points[P*2 + iter]);
    
    
    /*Populate triangles array*/
    T = mxGetM(prhs[1]);
    // mexPrintf("%d \n",T);
    //mexgetdoubles recieves doubles, this sets them to integers
    triangles =  mxGetDoubles(prhs[1]);
    // mexPrintf("%lf \n",triangles[0]);
    // for(iter = 0; iter < T;iter++)
    //mexPrintf("%f %f %f\n",triangles[T*0 + iter],triangles[T*1 + iter],triangles[T*2 + iter]);
    
    // i= 1;
    //int pQuad[9] = {(int)triangles[T*0 + i], (int)triangles[T*1 + i],(int)triangles[T*2 + i],(int)triangles[T*3 + i], (int)triangles[T*4 + i],
    //(int)triangles[T*5 + i],(int)triangles[T*6 + i],(int) triangles[T*7 + i],(int)triangles[T*8 + i]};
    //mexPrintf("%d %d %d %d %d %d %d %d %d\n",pQuad[0],pQuad[1],pQuad[2],pQuad[3],pQuad[4],pQuad[5],pQuad[6],pQuad[7],pQuad[8]);
    
    /*Get the number of edges*/
    N = mxGetScalar(prhs[2]);
    // mexPrintf("%d \n",N);
    
    /*Fetch the frequency*/
    freq = mxGetScalar(prhs[3]);
    // mexPrintf("%lf \n",freq);
    
    /*Get the observer and source map*/
    obs_map = mxGetInt32s(prhs[4]);
    src_map = mxGetInt32s(prhs[5]);
    num_obs = mxGetScalar(prhs[6]);
    num_src = mxGetScalar(prhs[7]);
    // mexPrintf("%d %d %d\n",obs_map[0],obs_map[1],obs_map[2]);
    //mexPrintf("%d %d %d\n",src_map[0],src_map[1],src_map[2]);
    MODE = mxGetScalar(prhs[8]);
    mexPrintf("MODE = %d\n", MODE);
    //mexPrintf("%d\n", num_src);
    /*-----------------make internal pointer z matrix-------------------------------*/
    Zmat = mxMalloc(sizeof(complex double) * num_obs * num_src);//the entire memory allocation of
    
    /*THERE IS A MORE INTELLIGENT WAY OF DOING THIS*/
    //the speed for this was checked it is not the cause of the problem
    //start = clock();
    for (i = 0; i < num_obs * num_src; i++)
    {
        
        Zmat[i] = 0.000 + 0.000*I;
        
    }
    
    //start = clock();
    /*==========================================================================================*/
    /*==========================================================================================*/
    MoM(freq, P, points, T, triangles, N, Zmat, obs_map, src_map, num_obs, num_src, MODE);
    /*==========================================================================================*/
    /*==========================================================================================*/
    //end = clock();
    
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    //mexPrintf("the array fill took %f seconds to execute \n", cpu_time_used);
    
    //you don't have to mxFree anything as the prhs and plhs are mxArray types which matlab
    //handles itself
    plhs[0]= mxCreateDoubleMatrix(num_obs, num_src, mxCOMPLEX);
    mxSetComplexDoubles(plhs[0], Zmat);
    
    //mxFree(Zmat);
    
    mexPrintf("\nIf it were done when 'tis done, then 'twere well it were done quickly\n");
    
    
}