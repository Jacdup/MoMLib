#ifndef triangle_functions
#define triangle_functions

#include "matrix_vector_functions.h"
#include "math.h"

double TriangleArea(double points[][3])
{
	double AB[3], AC[3], cp[3];
	VectorSubtract(points[0], points[1], AB);
	VectorSubtract(points[0], points[2], AC);
	VectorCross(AB, AC, cp);
	return 0.5*VectorSize(cp);
}

void TriangleCentre (double points[][3], double centre[])
{
	centre[0] = (points[0][0]+points[1][0]+points[2][0])/3.0;
	centre[1] = (points[0][1]+points[1][1]+points[2][1])/3.0;
	centre[2] = (points[0][2]+points[1][2]+points[2][2])/3.0;
}

void InverseBarycentric(double points[][3], double R[3], double Zeta[])
{
    //Vector v0 = b - a, v1 = c - a, v2 = p - a;
    double v0[3], v1[3],v2[3];
    VectorSubtract(points[0], points[1], v0);
    VectorSubtract(points[0],points[2], v1);
    VectorSubtract(points[0], R, v2);
    
    double d00 = dotProduct(v0,v0);
    double d01 = dotProduct(v0,v1);
    double d11 = dotProduct(v1,v1);
    double d20 = dotProduct(v2, v0);
    double d21 = dotProduct(v2,v1);
    double denom = (d00 * d11) - (d01 * d01);
    
    Zeta[0] = (d11 * d20 - d01 * d21) / denom;
    Zeta[1] = (d00 * d21 - d01 * d20) / denom;
    Zeta[2] = 1.0 - Zeta[0] - Zeta[1];
    
}

void ParametricMap(double points[][3], double Zeta[3], double ruv[]){
    //double ruv[3];
    int ii;

     for (ii=0; ii<3; ii++){
         ruv[ii] = points[0][ii]*Zeta[0] + points[1][ii]*Zeta[1] + points[2][ii]*Zeta[2];
     }
}

int TrianglesDontTouch(int *triA, int *triB)
{
	int i,j;
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			if (triA[i] == triB[j])
				return 0;
		}
	}
	return 1;
}

double DistanceToTriangle(double TriPoints[3][3], double obsPoint[3])
{
	double dist = 1e9;
	double temp = 0.0;
	int i;
	for (i=0; i<3; i++)
	{
		temp = Distance(TriPoints[i], obsPoint);
		if (temp < dist)
			dist = temp;
	}
	return dist;
}

//Accurate distance to triangle uses vertices side midpoints and centrepoint
double StrictDistanceToTriangle(double TriPoints[3][3], double obsPoint[3])
{
	double dist = 1e9;
	double temp = 0.0;

	double TriangleData[7][3] = {
					{0.0, 0.0, 0.0},//vert 1
					{0.0, 0.0, 0.0},//vert 2
					{0.0, 0.0, 0.0},//vert 3
					{0.0, 0.0, 0.0},//edge midpoint1
					{0.0, 0.0, 0.0},//edge midpoint2
					{0.0, 0.0, 0.0},//edge midpoint2
					{0.0, 0.0, 0.0}	//triangle centroid					
				};

	double tri_1_centroid[3] = {0.0, 0.0, 0.0};
	double T1_mid_1[3]= {0.0, 0.0, 0.0};
	double T1_mid_2[3]= {0.0, 0.0, 0.0};
	double T1_mid_3[3]= {0.0, 0.0, 0.0};
	
	int i;

	//Put the triangles in the data
	TriangleData[0][0] = TriPoints[0][0];
	TriangleData[0][1] = TriPoints[0][1];
	TriangleData[0][2] = TriPoints[0][2];


	TriangleData[1][0] = TriPoints[1][0];
	TriangleData[1][1] = TriPoints[1][1];
	TriangleData[1][2] = TriPoints[1][2];


	TriangleData[2][0] = TriPoints[2][0];
	TriangleData[2][1] = TriPoints[2][1];
	TriangleData[2][2] = TriPoints[2][2];
	
	//Triangle edge midpoints
	for (i = 0; i <3;i++) T1_mid_1[i] = (TriangleData[0][i]+TriangleData[1][i])/2;
	TriangleData[3][0] = T1_mid_1[0];
	TriangleData[3][1] = T1_mid_1[1];
	TriangleData[3][2] = T1_mid_1[2];
	
	for (i = 0; i <3;i++) T1_mid_2[i] = (TriangleData[1][i]+TriangleData[2][i])/2;
	TriangleData[4][0] = T1_mid_2[0];
	TriangleData[4][1] = T1_mid_2[1];
	TriangleData[4][2] = T1_mid_2[2];

	for (i = 0; i <3;i++) T1_mid_3[i] = (TriangleData[2][i]+TriangleData[0][i])/2;
	TriangleData[5][0] = T1_mid_3[0];
	TriangleData[5][1] = T1_mid_3[1];
	TriangleData[5][2] = T1_mid_3[2];
	
	//centroid
	TriangleCentre(TriPoints,tri_1_centroid);
	TriangleData[6][0] = tri_1_centroid[0];
	TriangleData[6][1] = tri_1_centroid[1];
	TriangleData[6][2] = tri_1_centroid[2];
	
	for (i=0; i<7; i++)
	{
		temp = Distance(TriangleData[i], obsPoint);
		if (temp < dist)
			dist = temp;
	}
	return dist;
}

//Find the longest edge of a triangle (its diameter)
double TriangleDiameter(double TriPoints[3][3])
{
	double diameter = 0.0;
	double temp = 0.0;
	int i = 0;
	int j[3] = {1,2,0};

	for (i=0; i<3; i++)
	{
		temp = Distance(TriPoints[i],TriPoints[j[i]]);
		//printf("temp = %lf \n",temp);
		if (temp > diameter)
		{
			diameter = temp; 
		}
	}
	return diameter;
}

//ApproxShortDistanceBetweenTriangles returns d_tilda, reduce to 36
//ApproxShortDistanceBetweenTriangles returns d_tilda, reduce to 36
double ApproxShortDistanceBetweenTriangles_36(double TriPoints_1[3][3], double TriPoints_2[3][3])
{
	double dist = 1e9;
	double temp = 0.0;
	double TriangleData[6][3] = {
					{0.0, 0.0, 0.0},//vert 1
					{0.0, 0.0, 0.0},//vert 2
					{0.0, 0.0, 0.0},//vert 3
					{0.0, 0.0, 0.0},//edge midpoint1
					{0.0, 0.0, 0.0},//edge midpoint2
					{0.0, 0.0, 0.0}//edge midpoint2	
				};
	double TriangleData2[6][3] = {
					{0.0, 0.0, 0.0},//vert 1
					{0.0, 0.0, 0.0},//vert 2
					{0.0, 0.0, 0.0},//vert 3
					{0.0, 0.0, 0.0},//edge midpoint1
					{0.0, 0.0, 0.0},//edge midpoint2
					{0.0, 0.0, 0.0}//edge midpoint2						
				};
	double T1_mid_1[3]= {0.0, 0.0, 0.0};
	double T1_mid_2[3]= {0.0, 0.0, 0.0};
	double T1_mid_3[3]= {0.0, 0.0, 0.0};
	double T2_mid_1[3]= {0.0, 0.0, 0.0};
	double T2_mid_2[3]= {0.0, 0.0, 0.0};
	double T2_mid_3[3]= {0.0, 0.0, 0.0};
	
	int i;
	int j;
	
	//this functions uses the triangle vertices, and edge midpoints and returns the smallest value
	

	//Put the triangles in the data
	TriangleData[0][0] = TriPoints_1[0][0];
	TriangleData[0][1] = TriPoints_1[0][1];
	TriangleData[0][2] = TriPoints_1[0][2];


	TriangleData[1][0] = TriPoints_1[1][0];
	TriangleData[1][1] = TriPoints_1[1][1];
	TriangleData[1][2] = TriPoints_1[1][2];


	TriangleData[2][0] = TriPoints_1[2][0];
	TriangleData[2][1] = TriPoints_1[2][1];
	TriangleData[2][2] = TriPoints_1[2][2];
	

	//Put the triangles in the data
	TriangleData2[0][0] = TriPoints_2[0][0];
	TriangleData2[0][1] = TriPoints_2[0][1];
	TriangleData2[0][2] = TriPoints_2[0][2];


	TriangleData2[1][0] = TriPoints_2[1][0];
	TriangleData2[1][1] = TriPoints_2[1][1];
	TriangleData2[1][2] = TriPoints_2[1][2];


	TriangleData2[2][0] = TriPoints_2[2][0];
	TriangleData2[2][1] = TriPoints_2[2][1];
	TriangleData2[2][2] = TriPoints_2[2][2];
	//find edge midpoints
	
	//Triangle 1 edge midpoints
	for (i = 0; i <3;i++) T1_mid_1[i] = (TriangleData[0][i]+TriangleData[1][i])/2;
	TriangleData[3][0] = T1_mid_1[0];
	TriangleData[3][1] = T1_mid_1[1];
	TriangleData[3][2] = T1_mid_1[2];
	
	for (i = 0; i <3;i++) T1_mid_2[i] = (TriangleData[1][i]+TriangleData[2][i])/2;
	TriangleData[4][0] = T1_mid_2[0];
	TriangleData[4][1] = T1_mid_2[1];
	TriangleData[4][2] = T1_mid_2[2];
	
	for (i = 0; i <3;i++) T1_mid_3[i] = (TriangleData[2][i]+TriangleData[0][i])/2;
	TriangleData[5][0] = T1_mid_3[0];
	TriangleData[5][1] = T1_mid_3[1];
	TriangleData[5][2] = T1_mid_3[2];
	
	//Triangle 2 edge midpoints
	for (i = 0; i <3;i++) T2_mid_1[i] = (TriangleData2[0][i]+TriangleData2[1][i])/2;
	TriangleData2[3][0] = T2_mid_1[0];
	TriangleData2[3][1] = T2_mid_1[1];
	TriangleData2[3][2] = T2_mid_1[2];
	
	for (i = 0; i <3;i++) T2_mid_2[i] = (TriangleData2[1][i]+TriangleData2[2][i])/2;
	TriangleData2[4][0] = T2_mid_2[0];
	TriangleData2[4][1] = T2_mid_2[1];
	TriangleData2[4][2] = T2_mid_2[2];
	
	for (i = 0; i <3;i++) T2_mid_3[i] = (TriangleData2[2][i]+TriangleData2[0][i])/2;
	TriangleData2[5][0] = T2_mid_3[0];
	TriangleData2[5][1] = T2_mid_3[1];
	TriangleData2[5][2] = T2_mid_3[2];
	

	
	//all the processing of the triangles is now complete, we have the 7 points from each
	//remember to check for plates that are perfectly parallel
	for(i = 0; i<6; i++)
	{
		
		for (j = 0; j<6; j++)
		{
			temp = Distance(TriangleData[i], TriangleData2[j]);
			if (temp < dist)
				dist = temp;		
		}
			
	}
	


	return dist;
	
}


//dont know why this is called d_9 it finds shortest distance between the 6 corners
double ApproxShortDistanceBetweenTriangles_9(double TriPoints_1[3][3], double TriPoints_2[3][3])
{
	double dist = 1e9;
	double temp = 0.0;
	int i;
	int j;
	
	//finds the shortest distance between the 9 centroids
	for(i = 0; i<3; i++)
	{
		
		for (j = 0; j<3; j++)
		{
			temp = Distance(TriPoints_1[i], TriPoints_2[j]);
			if (temp < dist)
				dist = temp;		
		}
			
	}
	
	return dist;
	
}

//added 15/07/2019
//create a function that checks whether a triangle is using a 
int ObserverOnTriangle(int basisindex[3], int *obs_map)//not a problem because its passed by reference
{
	int i;
	for (i = 0; i <3;i++)
	{
		//index the obs map using the basis functions
		
		if (basisindex[i] != -1)
		{
			//printf("basisindex = %d\n", basisindex[i]);
			//we get in here cuz there is a basis function referencing a observer location
			//printf("pEdgeIndex = %d\n", obs_map[basisindex[i]-1]);
			if (obs_map[basisindex[i]-1] != -1)//so if any value on the obs map is not equal to -1
			{	
				//printf("pEdgeIndex = %d\n", obs_map[basisindex[i]-1]);
				//return 1 because there is an observer defined in the triangle
				return 1;
			}
		}	
	}
	//if we get to this point there is no observer on the triangle and we return zero
	return 0;

}

#endif