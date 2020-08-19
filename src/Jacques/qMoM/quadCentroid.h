#ifndef quadCentroid
#define quadCentroid

#include "math.h"
#include "matrix_vector_functions.h"
#include "mex.h"


#define EPS 1e-13 // Make this as small as you like. 
#define FALSE 0
#define TRUE 1

void TriangleCentre (double points[3][3], double centre[])
{
    centre[0] = (points[0][0]+points[1][0]+points[2][0])/3.0; //x
    centre[1] = (points[0][1]+points[1][1]+points[2][1])/3.0; //y
    centre[2] = (points[0][2]+points[1][2]+points[2][2])/3.0; //z
}


int LineIntersect(double p1[3],double p2[3],double p3[3],double p4[3],double pa[],double pb[])
{
    /* See:
     * http://paulbourke.net/geometry/pointlineplane/ 
     * for algorithm */
    
    double p13[3];
    double p43[3];
    double p21[3];
    double mua, mub;
    double d1343,d4321,d1321,d4343,d2121;
    double numer,denom;
    
    VectorSubtract(p3,p1,p13);
    VectorSubtract(p3,p4,p43);
    
    
    if (fabs(p43[0]) < EPS && fabs(p43[1]) < EPS && fabs(p43[2]) < EPS)
        return(FALSE);
    VectorSubtract(p1,p2,p21);
    
    if (fabs(p21[0]) < EPS && fabs(p21[1]) < EPS && fabs(p21[2]) < EPS)
        return(FALSE);
    
    d1343 = dotProduct(p13,p43);
    d4321 = dotProduct(p43,p21);
    d1321 = dotProduct(p13,p21);
    d4343 = dotProduct(p43,p43);
    d2121 = dotProduct(p21,p21);

    denom = d2121 * d4343 - d4321 * d4321;
    if (fabs(denom) < EPS){
        mexPrintf("Denom = %f\n", denom);
        // This will go off if mesh is sufficiently small (elements to the order EPS)
        return(FALSE);
    }
    numer = d1343 * d4321 - d1321 * d4343;
    
    mua = numer / denom;
    mub = (d1343 + d4321 * (mua)) / d4343;

    int i;
    for (i=0; i<3; i++){
        pa[i] = p1[i] + mua*p21[i];
        pb[i] = p3[i] + mub*p43[i];
    }

    return(TRUE);
}



void QuadCentre(double points[][3], double centre[])
{
    /* 
      * Construct 2 triangles inside quadrilateral
      * Find centroid of both
      * Find other 2 triangles and centroids
      * Quadrilateral centroid is point where lines intersect
     */

    //Centroids:
    double c1[3] = {0.0, 0.0, 0.0};
    double c2[3] = {0.0, 0.0, 0.0};
    double c3[3] = {0.0, 0.0, 0.0};
    double c4[3] = {0.0, 0.0, 0.0};
    
    //Triangles
    double tri_1[3][3] = {  {points[0][0],points[0][1],points[0][2]},
    {points[1][0],points[1][1],points[1][2]},
    {points[2][0],points[2][1],points[2][2]}};
    double tri_2[3][3] = {  {points[1][0],points[1][1],points[1][2]},
    {points[2][0],points[2][1],points[2][2]},
    {points[3][0],points[3][1],points[3][2]}};
    double tri_3[3][3] = {  {points[0][0],points[0][1],points[0][2]},
    {points[2][0],points[2][1],points[2][2]},
    {points[3][0],points[3][1],points[3][2]}};
    double tri_4[3][3] = {  {points[0][0],points[0][1],points[0][2]},
    {points[1][0],points[1][1],points[1][2]},
    {points[3][0],points[3][1],points[3][2]}};
    
    TriangleCentre(tri_1, c1);
    TriangleCentre(tri_2, c2);
    TriangleCentre(tri_3, c3);
    TriangleCentre(tri_4, c4);
    
    double pa[3] = {0.0, 0.0, 0.0};
    double pb[3] = {0.0, 0.0, 0.0};
    
    int pointExists = 0;
    pointExists = LineIntersect(c1,c3,c2,c4,pa,pb);
    if (pointExists == 0){
        mexPrintf("Cannot find centroid of quad\n");
        // This will go off if mesh is sufficiently small (elements to the order EPS)
    }
   // if (pointExists == 1){
        // mexPrintf("pa = %f,%f,%f\n", pa[0],pa[1],pa[2]);
        // mexPrintf("pb = %f,%f,%f\n", pb[0],pb[1],pb[2]);
   // }
    // pa and pb should be the same point
    
    centre[0] = pa[0];
    centre[1] = pa[1];
    centre[2] = pa[2];
    
}


#endif