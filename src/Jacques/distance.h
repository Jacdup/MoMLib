#ifndef distance
#define distance

#include "math.h"
#include "quadCentroid.h"
#include "mex.h"

/*
 * check if the edge defined by [val1, val2] is a non-boundary edge.
 * Return index+1, else return -1
 */
int IsEdgeInMat(int val1, int val2, int **edges, int N){
    int i;
    for (i=0; i < N; i++)
    {
        if (edges[i][0] == (int)fmin(val1, val2) && edges[i][1] == (int)fmax(val1,val2))
            return i;
    }
    return -1;
}

int QuadsDontTouch(int *quadA, int *quadB)
{
    int i, j;
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            if (quadA[i] == quadB[j])
                return 0;
        }
    }
    return 1;
}

int ElementsDontTouch(int *elemA, int *elemB)
{
    // Function for arbitrary order elements
    int i,j;
    int Ni = sizeof(elemA)/sizeof(elemA[0]);
    int Nj = sizeof(elemB)/sizeof(elemB[0]);
    for (i=0; i<Ni; i++)
    {
        for (j=0; j<Nj; j++)
        {
            if (elemA[i] == elemB[j])
                return 0;
        }
    }
    return 1;
    
    
}
double DistanceToQuad(double QuadPoints[4][3], double obsPoint[3])
{
    double dist = 1e9;
    double temp = 0.0;
    int i;
    for (i=0; i<4; i++)
    {
        temp = Distance(QuadPoints[i], obsPoint);
        if (temp < dist)
            dist = temp;
    }
    return dist;
}

//-------------------------NEW ROUTINES--------------------------------
double StrictDistance(double QuadPoints[4][3], double obsPoint[3])
{ // StrictDistance to quadrilateral
    
    double dist = 1e9;
    double temp = 0.0;
    double **element_data = mxMalloc(sizeof(double *) * 9);
    double centroid[3] = {0.0, 0.0, 0.0};
    double T1_mid_1[3]= {0.0, 0.0, 0.0};
    double T1_mid_2[3]= {0.0, 0.0, 0.0};
    double T1_mid_3[3]= {0.0, 0.0, 0.0};
    double T1_mid_4[3]= {0.0, 0.0, 0.0};
    
    
    int i;
    //Vertices
    for (i=0; i<4; i++){
        element_data[i] = QuadPoints[i];
    }
    
    // Edge Midpoints
    for (i=0; i<3; i++){ // x,y,z iterator
        T1_mid_1[i] = (QuadPoints[0][i] + QuadPoints[1][i])/2;
        T1_mid_2[i] = (QuadPoints[0][i] + QuadPoints[2][i])/2;
        T1_mid_3[i] = (QuadPoints[1][i] + QuadPoints[3][i])/2;
        T1_mid_4[i] = (QuadPoints[2][i] + QuadPoints[3][i])/2;
    }
    element_data[4] = T1_mid_1;
    element_data[5] = T1_mid_2;
    element_data[6] = T1_mid_3;
    element_data[7] = T1_mid_4;
    
    //Centroid
    QuadCentre(QuadPoints, centroid);
    element_data[8] = centroid;
    
    
    for (i=0; i<9; i++)
    {
        temp = Distance(element_data[i], obsPoint);
        if (temp < dist)
            dist = temp;
    }
    mxFree(element_data);
    return dist;
    
}

double QuadDiameter(double QuadPoints[4][3])
{
    double diameter = 0.0;
    double temp = 0.0;
    int i;
    int j[4] = {1,3,0,2};
    
    for (i=0; i<4; i++)
    {
        temp = Distance(QuadPoints[i],QuadPoints[j[i]]);
        if (temp > diameter)
        {
            diameter = temp;
        }
    }
    return diameter;
}

double ApproxShortDistanceBetweenElements_36(double QuadPoints_1[4][3], double QuadPoints_2[4][3])
{
    
    double dist = 1e9;
    double temp = 0.0;
    double **element_data_1 = mxMalloc(sizeof(double *) * 8);
    double **element_data_2 = mxMalloc(sizeof(double *) * 8);
    double centroid[3] = {0.0, 0.0, 0.0};
    double T1_mid_1[3]= {0.0, 0.0, 0.0};
    double T1_mid_2[3]= {0.0, 0.0, 0.0};
    double T1_mid_3[3]= {0.0, 0.0, 0.0};
    double T1_mid_4[3]= {0.0, 0.0, 0.0};
    double T2_mid_1[3]= {0.0, 0.0, 0.0};
    double T2_mid_2[3]= {0.0, 0.0, 0.0};
    double T2_mid_3[3]= {0.0, 0.0, 0.0};
    double T2_mid_4[3]= {0.0, 0.0, 0.0};
    
    int i,j;
    
    //Vertices
    for (i=0; i<4; i++)
    {
        element_data_1[i] = QuadPoints_1[i];
        element_data_2[i] = QuadPoints_2[i];
    }
    
    // Edge Midpoints
    for (i=0; i<3; i++){ // x,y,z iterator
        T1_mid_1[i] = (QuadPoints_1[0][i] + QuadPoints_1[1][i])/2;
        T1_mid_2[i] = (QuadPoints_1[0][i] + QuadPoints_1[2][i])/2;
        T1_mid_3[i] = (QuadPoints_1[1][i] + QuadPoints_1[3][i])/2;
        T1_mid_4[i] = (QuadPoints_1[2][i] + QuadPoints_1[3][i])/2;
        
        T2_mid_1[i] = (QuadPoints_2[0][i] + QuadPoints_2[1][i])/2;
        T2_mid_2[i] = (QuadPoints_2[0][i] + QuadPoints_2[2][i])/2;
        T2_mid_3[i] = (QuadPoints_2[1][i] + QuadPoints_2[3][i])/2;
        T2_mid_4[i] = (QuadPoints_2[2][i] + QuadPoints_2[3][i])/2;
    }
    element_data_1[4] = T1_mid_1;
    element_data_1[5] = T1_mid_2;
    element_data_1[6] = T1_mid_3;
    element_data_1[7] = T1_mid_4;
    element_data_2[4] = T2_mid_1;
    element_data_2[5] = T2_mid_2;
    element_data_2[6] = T2_mid_3;
    element_data_2[7] = T2_mid_4;
    //this functions uses the quadrilateral centroids, vertices, and edge midpoints and returns the smallest value
    
    //all the processing of the elements is now complete, we have the 9 points from each
    //remember to check for plates that are perfectly parallel
    for(i = 0; i<8; i++)
    {
        for (j = 0; j<8; j++)
        {
            temp = Distance(element_data_1[i], element_data_2[j]);
            if (temp < dist)
                dist = temp;
        }
    }
    
    mxFree(element_data_1);
    mxFree(element_data_2);
    return dist;
    
}


double ApproxShortDistanceBetweenElements_9(double QuadPoints_1[4][3], double QuadPoints_2[4][3])
{
    double dist = 1e9;
    double temp = 0.0;
    int i;
    int j;
    
    //finds the shortest distance between the 16 vertices
    for(i = 0; i<4; i++)
    {
        for (j = 0; j<4; j++)
        {
            temp = Distance(QuadPoints_1[i], QuadPoints_2[j]);
            if (temp < dist)
                dist = temp;
        }
    }
    
    return dist;
    
}

int ObserverOnElement(int basisindex[4], int *obs_map)
{
    int i;
    for (i=0; i<4; i++)
    {
        //index the obs map using the basis functions
        
        if (basisindex[i] != -1)
        {
            //we get in here cuz there is a basis function referencing a observer location
            if (obs_map[basisindex[i]-1] != -1)//so if any value on the obs map is not equal to -1
            {
                //return 1 because there is an observer defined in the quadrilateral
                return 1;
            }
        }
    }
    //if we get to this point there is no observer on the quadrilateral and we return zero
    return 0;
    
}

#endif