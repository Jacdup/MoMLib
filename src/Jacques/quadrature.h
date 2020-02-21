#ifndef quadrature
#define quadrature

#include "matrix_vector_functions.h"
#include "math.h"
#include "mex.h"
#include "quadCentroid.h"

void GaussianQuadrature(double points[][3], int numIntPoints, double **returnPoints)
{

    /*
     * Computes quadrature points on a quadrilateral defined by points.
     * Returns points and weight.
     */
    int i, j, ii,k;
    
    // Define on interval [-1,1]
    // Formatted as {x1,x2,...,xn}
    //              {w1,w2,...,wn}};
    double quadPoints1[2][1] = {{0},
    {2}};
    double quadPoints2[2][2] = {{-0.577350269189626,0.577350269189626},
    {1,1}};
    double quadPoints3[2][3] = {{-0.774596669241483,0,0.774596669241483},
    {0.555555555555556,0.888888888888889,0.555555555555556}};
    double quadPoints4[2][4] = {{-0.861136311594054,-0.339981043584857,0.339981043584856,0.861136311594053},
    {0.347854845137454,0.652145154862545,0.652145154862546,0.347854845137454}};
    double quadPoints5[2][5] = {{-0.906179845938664,-0.538469310105683,0,0.538469310105683,0.906179845938664},
    {0.236926885056189,0.478628670499367,0.568888888888889,0.478628670499366,0.236926885056189}};
    double quadPoints6[2][6] = {{-0.932469514203153,-0.661209386466264,-0.238619186083197,0.238619186083197,0.661209386466264,0.932469514203152},
    {0.171324492379169,0.360761573048141,0.467913934572689,0.467913934572691,0.360761573048138,0.171324492379171}};
    double quadPoints7[2][7] = {{-0.949107912342759,-0.741531185599395,-0.405845151377397,0,0.405845151377397,0.741531185599396,0.949107912342757},
    {0.129484966168869,0.279705391489277,0.381830050505120,0.417959183673469,0.381830050505122,0.279705391489274,0.129484966168871}};
    double quadPoints8[2][8] = {{-0.960289856497537,-0.796666477413630,-0.525532409916328,-0.183434642495650,0.183434642495650,0.525532409916329,0.796666477413627,0.960289856497536},
    {0.101228536290374,0.222381034453376,0.313706645877890,0.362683783378361,0.362683783378360,0.313706645877890,0.222381034453372,0.101228536290377}};
    double quadPoints13[2][13] = {{-0.984183054718600,-0.917598399222968,-0.801578090733310,-0.642349339440349,-0.448492751036441,-0.230458315955136,0,0.230458315955135,0.448492751036442,0.642349339440350,0.801578090733302,0.917598399222974,0.984183054718594},
    {0.0404840047653037,0.0921214998377659,0.138873510219737,0.178145980761987,0.207816047536865,0.226283180262914,0.232551553234957,0.226283180262969,0.207816047536827,0.178145980762005,0.138873510219732,0.0921214998377660,0.0404840047653028}};
    
    
    double (*qpPointer)[2][13];
    int qpNum;
    
    
//	if (numIntPoints > 16 && numIntPoints < 25)
//		numIntPoints = 16;
    
    switch (numIntPoints)
    {
        case 1:qpPointer = &quadPoints1;
        qpNum = 1;
        break;
        case 2:
            qpPointer = &quadPoints1;
            qpNum = 1;
            break;
        case 3:
            //qpPointer = &quadPoints3;
            //qpNum = 3;
            //break;
        case 4:
            qpPointer = &quadPoints2;
            qpNum = 2;
            break;
        case 5:
            //qpPointer = &quadPoints4;
            //qpNum = 4;
            //break;
        case 6:
            //qpPointer = &quadPoints6;
            //qpNum = 6;
            //break;
        case 7:
        case 8:
        case 9:
            qpPointer = &quadPoints3;
            qpNum = 3;
            break;
        case 10:
        case 11:
            //qpPointer = &quadPoints7;
            //qpNum = 7;
            //break;
        case 12:
            //qpPointer = &quadPoints8;
            //qpNum = 8;
            //break;
        case 13:
        case 14:
        case 15:
            //qpPointer = &quadPoints13;
            //qpNum = 13;
            //break;
        case 16:
            qpPointer = &quadPoints4;
            qpNum = 4;
            break;
        case 25:
            qpPointer = &quadPoints5;
            qpNum = 5;
            break;
        case 36:
            qpPointer = &quadPoints6;
            qpNum = 6;
            break;
        case 49:
            qpPointer = &quadPoints7;
            qpNum = 7;
            break;
        default:
            qpPointer = &quadPoints13;
            qpNum = 13;
            
    }
    
    int jj, coord;
    int place = 0;
    int index = 0;
    int iter = pow(qpNum,2);
    double weightProd[iter];
    double quadProd[iter];
    
    index = 0;
    
   // double Area = 0;

    double simplexCoord[2];
    double Q;
    
    for (ii=0; ii<qpNum; ii++){ // u
        for (jj=0; jj<qpNum; jj++){ // v
            index = (ii*qpNum) + jj;
            weightProd[index] = (*qpPointer)[0][ii+qpNum]*(*qpPointer)[0][jj+qpNum]; // Product Rule on the weights
           // Area += weightProd[index];

            double b10[3],b01[3],b11[3],ruv[3],drdv[3],drdu[3];
            double J[3]; //Jacobian
            
            
            for (coord=0; coord<3; coord++){
                
                //  b10[coord] = (-points[0][coord] + points[1][coord]- points[3][coord] + points[2][coord]);
                //   b01[coord] = (-points[0][coord] - points[1][coord]+ points[3][coord] + points[2][coord]);
                //  drdv[coord] = 0.5*b01[coord];
                //  drdu[coord] = 0.5*b10[coord];

                b10[coord] = -points[0][coord] + points[1][coord];
                b01[coord] = -points[0][coord] + points[3][coord] ;
                b11[coord] = (points[0][coord] - points[1][coord]- points[3][coord] + points[2][coord]);
                drdv[coord] = (b01[coord] + b11[coord]*(*qpPointer)[0][ii]);
                drdu[coord] = (b10[coord] + b11[coord]*(*qpPointer)[0][jj]);//b10*u
            }

            VectorCross(drdv,drdu,J);
            //double J_det = sqrt(pow(drdu[1]*drdv[2] - drdu[2]*drdv[1],2) + pow(drdu[2]*drdv[0] - drdu[0]*drdv[2],2) + pow(drdu[0]*drdv[1] - drdu[1]*drdv[0],2));
            double J_det = sqrt(pow(J[0],2) + pow(J[1],2) +pow(J[2],2));
            // mexPrintf("t1 = %f,t2 = %f,t3 = %f\n", t[1],t[2],t[3]);
            // for (coord=0; coord<3; coord++)
            //  {
            //returnPoints[place][coord] = ruv[coord];//1/Q * t[coord];// //
            returnPoints[place][0] = (1+(*qpPointer)[0][ii]) *0.5 ;//ruv[coord];//1/Q * t[coord];// //
            returnPoints[place][1] = (1+(*qpPointer)[0][jj])* 0.5;
            returnPoints[place][2] = 0.0;
            //   }
            
            //J_det = 1;
            returnPoints[place][3] = J_det * 0.25 * weightProd[index] ;//The multiplicative factor because of change of interval
            place++;
        } //jj
        
    }
    // returnPoints[1][1] = Area;
}

/*
 * Functions for Radial Angular R1 Sqrt
 */
double wFromXY (double x, double y)
{
    return asinh(x/y); //p
}
double qFromYWZ (double y, double w, double z)
{
    return sqrt(sqrt(pow(cosh(w),2)*y*y + z*z) - fabs(z));
}
double yFromQWZ (double q, double w, double z)
{
    return sqrt((pow((q*q + fabs(z)),2) - z*z)/pow(cosh(w),2));
}
double xFromYW (double y, double w)
{
    return y*sinh(w);
}

/*
 * Calculates the integral on subtriangles.
 * The assuption is that the third vertice, v3, is the projection of the singularity
 */
void RAR1S_2D(double v1[3], double v2[3], double v3[3], double crux, int numIntPoints, double **returnPoints)
{
    int i,j;
    
    //Gaussian Quadrature rules in 1D
    double quadPoints1[1][2] = {{2,0}};
    double quadPoints2[2][2] = {{1.0000000000000000,    -0.5773502691896257},
    {1.0000000000000000,    0.5773502691896257}};
    double quadPoints3[3][2] = {{0.8888888888888888,    0.0000000000000000},
    {0.5555555555555556,    -0.7745966692414834},
    {0.5555555555555556,    0.7745966692414834}};
    double quadPoints4[4][2] = {{0.6521451548625461,    -0.3399810435848563},
    {0.6521451548625461,    0.3399810435848563},
    {0.3478548451374538,    -0.8611363115940526},
    {0.3478548451374538,    0.8611363115940526}};
    double quadPoints5[5][2] = {{0.5688888888888889,    0.0000000000000000},
    {0.4786286704993665,    -0.5384693101056831},
    {0.4786286704993665,    0.5384693101056831},
    {0.2369268850561891,    -0.9061798459386640},
    {0.2369268850561891,    0.9061798459386640}};
    double quadPoints6[6][2] = {{0.3607615730481386,    0.6612093864662645},
    {0.3607615730481386,    -0.6612093864662645},
    {0.4679139345726910,    -0.2386191860831969},
    {0.4679139345726910,    0.2386191860831969},
    {0.1713244923791704,    -0.9324695142031521},
    {0.1713244923791704,    0.9324695142031521}};
    double quadPoints7[7][2] = {{0.4179591836734694,    0.0000000000000000},
    {0.3818300505051189,    0.4058451513773972},
    {0.3818300505051189,    -0.4058451513773972},
    {0.2797053914892766,    0.7415311855993945},
    {0.2797053914892766,    -0.7415311855993945},
    {0.1294849661688697,    0.9491079123427585},
    {0.1294849661688697,    -0.9491079123427585}};
    double quadPoints8[8][2] = {{0.10122853629038, -0.96028985649754},
    {0.22238103445338, -0.79666647741362},
    {0.31370664587788, -0.52553240991633},
    {0.36268378337836, -0.18343464249565},
    {0.36268378337836, 0.18343464249565},
    {0.31370664587788, 0.52553240991633},
    {0.22238103445338, 0.79666647741362},
    {0.10122853629038, 0.96028985649754}};
    double quadPoints9[9][2] = {{0.08127438836157, -0.96816023950763},
    {0.18064816069487, -0.83603110732663},
    {0.26061069640291, -0.61337143270059},
    {0.31234707704002, -0.32425342340381},
    {0.33023935500125, 0.0},
    {0.31234707704001, 0.32425342340381},
    {0.26061069640292, 0.61337143270059},
    {0.18064816069487, 0.83603110732663},
    {0.08127438836157, 0.96816023950763}};
    
    double quadPoints10[10][2] = {{0.06667134430869, -0.97390652851717},
    {0.14945134915058, -0.86506336668899},
    {0.21908636251599, -0.67940956829902},
    {0.26926671930998, -0.43339539412925},
    {0.29552422471476, -0.14887433898163},
    {0.29552422471475, 0.14887433898163},
    {0.26926671930999, 0.43339539412925},
    {0.21908636251599, 0.67940956829902},
    {0.14945134915058, 0.86506336668899},
    {0.06667134430869, 0.97390652851717}};
    

    //Gaussian Quadrature rules in 1D
    
    int qpNum;
    double (*qpPointer)[10][2];
    
    switch (numIntPoints)
    {
        case 1:
            qpPointer = &quadPoints1;
            qpNum=1;
            break;
        case 2:
            qpPointer = &quadPoints2;
            qpNum=2;
            break;
        case 3:
            qpPointer = &quadPoints3;
            qpNum=3;
            break;
        case 4:
            qpPointer = &quadPoints4;
            qpNum=4;
            break;
        case 5:
            qpPointer = &quadPoints5;
            qpNum=5;
            break;
        case 6:
            qpPointer = &quadPoints6;
            qpNum=6;
            break;
        case 7:
            qpPointer = &quadPoints7;
            qpNum=7;
            break;
        case 8:
            qpPointer = &quadPoints8;
            qpNum=8;
            break;
        case 9:
            qpPointer = &quadPoints9;
            qpNum=9;
            break;
        case 10:
        default:
            qpNum=10;
            qpPointer = &quadPoints10;
            break;
    }
    
    double tempV1[3], newV1[3], tempV2[3], newV2[3], pVec[3];
    
    VectorSubtract(v3, v1, tempV1);		//set origin
    VectorSubtract(v3, v2, tempV2);
    VectorSubtract(tempV1, tempV2, pVec);   //non-origin edge as a vector
    
    double theta;      //angle to get edge parallel to x-axis
    if (pVec[0]==0)
        theta = M_PI*0.5;
    else
        theta = atan(-1.0*pVec[1]/pVec[0]);
    
    double rotationMatrix[3][3] = {{cos(theta), -1.0*sin(theta),0},
    {sin(theta), cos(theta),0},
    {0,0,1}};
    
    MatrixMultiplyVector(rotationMatrix, 3, tempV1, newV1);
    MatrixMultiplyVector(rotationMatrix, 3, tempV2, newV2); 		//rotate non-origin points
    
    
    int flip = 1.0;
    if (newV1[1] < 0)           //flip triangle if upside down
        flip = -1.0;
    newV1[1] *= flip;
    newV2[1] *= flip;
    
    rotationMatrix[0][1] *= -1;
    rotationMatrix[1][0] *= -1;                        //rotationMatrix will now undo rotation
    
    double *leftPoint, *rightPoint;
    if (newV1[0] < newV2[0])
    {
        leftPoint = &newV1[0];
        rightPoint = &newV2[0];
    }
    else
    {
        leftPoint = &newV2[0];              //left most point of edge for integration lower limit
        rightPoint = &newV1[0];
    }
    
    double w_lower = wFromXY(leftPoint[0], leftPoint[1]);
    double w_upper = wFromXY(rightPoint[0], rightPoint[1]);
    double w_range = w_upper - w_lower;
    
    for (i=0; i < qpNum; i++)
    {
        double wPoint = w_lower + 0.5*w_range*(1 + (*qpPointer)[i][1]);
        
        for (j=0; j < qpNum; j++)
        {
            double q_lower = qFromYWZ(0,wPoint,fabs(crux));
            double q_upper = qFromYWZ(leftPoint[1], wPoint, crux);
            double q_range = q_upper - q_lower;
            
            double qPoint = q_lower + 0.5*q_range*(1 + (*qpPointer)[j][1]);
            double yPoint = yFromQWZ(qPoint, wPoint, crux);
            double xPoint = xFromYW(yPoint, wPoint);
            
            double tempNewPoint[3] = {xPoint, flip*yPoint, 0.0};
            MatrixMultiplyVector(rotationMatrix, 3, tempNewPoint, returnPoints[i*qpNum + j]);

            returnPoints[i*qpNum + j][0] += v3[0];
            returnPoints[i*qpNum + j][1] += v3[1];
            returnPoints[i*qpNum + j][2] += v3[2];//rotate and translate back to original triangle coordinate system
            
            returnPoints[i*qpNum + j][3] = (*qpPointer)[j][0]*(*qpPointer)[i][0]*q_range*w_range*0.25*2*qPoint/cosh(wPoint);   //append weighting of quadpoint
        }
    }
}


void RAR1S(double points[][3], double obsPoint[], int numIntPoints, double **returnPoints)
{
    int i, j;
    double triNorm[3], AB[3], AC[3];
    
    VectorSubtract(points[0], points[1], AB);
    VectorSubtract(points[0], points[2], AC);
    VectorCross(AB, AC, triNorm);				//Find norm of triangle
    
    if (triNorm[2] < 0) //z
        MultiplyVectorByConstant(triNorm, 3, -1.0);
    MultiplyVectorByConstant(triNorm, 3, 1/VectorSize(triNorm));
    
    double thetaNorm = acos(triNorm[2]);
    double phiNorm;
    
    /*
     * 		ThetaNorm and PhiNorm are the required rotation angles not the actual angles in polar coordinates of the norm vector
     */
    if (triNorm[1]== 0)
    {
        if (triNorm[0]>0) //x
            phiNorm = M_PI/2;
        else
            phiNorm = -M_PI/2;
    }
    else
        phiNorm = atan(triNorm[0]/triNorm[1]);
    
    if (triNorm[0]<0 && triNorm[1]<0)
        phiNorm += M_PI;
    else if (triNorm[1]<0)
        phiNorm -= M_PI;
    
    double Rz[3][3] = {{cos(phiNorm), -sin(phiNorm),0},
    {sin(phiNorm), cos(phiNorm), 0},
    {0, 0, 1}};
    
    double Rx[3][3] = {{1, 0, 0},
    {0, cos(thetaNorm), -sin(thetaNorm)},
    {0, sin(thetaNorm), cos(thetaNorm)}};
    
    double rotationMatrix[3][3];
    MatrixMultiplyMatrix(Rx, Rz, rotationMatrix);
    double rotationMatrix_inv[3][3];
    MatrixInverse(rotationMatrix, 3, rotationMatrix_inv);
    
    double newPoints[4][3];
    
    for (i=0; i< 4; i++)
        MatrixMultiplyVector(rotationMatrix, 3, points[i], newPoints[i]);  //quadrilateral is now rotated
    
    
    double zOffset = newPoints[0][2];
    for (i=0; i< 4; i++)
        newPoints[i][2]=0.0;
    
    double newObsPoint[3];
    MatrixMultiplyVector(rotationMatrix, 3, obsPoint, newObsPoint);       //rotate observation point
    double zObsPoint = newObsPoint[2];
    newObsPoint[2]=0;
    
    double quadCentre[3];
    QuadCentre(newPoints, quadCentre);

    for (i=0; i<4; i++)          //split quad into four triangles at the projection of the singularity
    {
        // newPoints -> rotated quad
        double AB[3]; // Side AB
        VectorSubtract(newPoints[i], newPoints[(i+1)%4], AB);
        
        double AtoCentre[3];
        VectorSubtract(newPoints[i], quadCentre, AtoCentre);
        
        double AtoObs[3];
        VectorSubtract(newPoints[i], newObsPoint, AtoObs);
        
        double vec1[3], vec2[3];
        VectorCross(AB, AtoCentre, vec1);
        VectorCross(AB, AtoObs, vec2);
        
        int newNumInt = numIntPoints;
        if (newNumInt>10)
            newNumInt=10;
        
        double **returnPoints2D = NULL;
        returnPoints2D = mxMalloc(sizeof(double *) * (int)pow(newNumInt,2));
        for (j=0; j < pow(newNumInt, 2); j++)
            returnPoints2D[j] = (double *)mxMalloc(sizeof(double) * 4);
        
        RAR1S_2D(newPoints[i], newPoints[(i+1)%4], newObsPoint, zObsPoint-zOffset, newNumInt, returnPoints2D);  //get integration points
        
        int posTri = 1;
        if (vec1[2]*vec2[2] < 0)
            posTri = -1; // Subtract triangle that's outside domain
        //mexPrintf("i = %d\n", i);
        newObsPoint[2] += zObsPoint;
        for (j=0; j < pow(newNumInt,2); j++)
        {
            returnPoints2D[j][2] += zOffset;
            double R = Distance(returnPoints2D[j], newObsPoint);
            MatrixMultiplyVector(rotationMatrix_inv, 3, returnPoints2D[j], returnPoints[(int)pow(newNumInt,2)*i + j]);
            returnPoints[(int)pow(newNumInt,2)*i + j][3] = posTri*R*returnPoints2D[j][3];
        }
        newObsPoint[2] -= zObsPoint;
        for (j=0; j < pow(newNumInt,2); j++)
            mxFree(returnPoints2D[j]);
        mxFree(returnPoints2D);
    }
}

void RAR1STri(double points[][3], double obsPoint[], int numIntPoints, double **returnPoints)
{
	int i, j;
    double triNorm[3], AB[3], AC[3];

    VectorSubtract(points[0], points[1], AB);
    VectorSubtract(points[0], points[2], AC);
    VectorCross(AB, AC, triNorm);				//Find norm of triangle

    if (triNorm[2] < 0)
        MultiplyVectorByConstant(triNorm, 3, -1.0);
    MultiplyVectorByConstant(triNorm, 3, 1/VectorSize(triNorm));

    double thetaNorm = acos(triNorm[2]);
    double phiNorm;

    /*
     * 		ThetaNorm and PhiNorm are the required rotation angles not the actual angles in polar coordinates of the norm vector
     */

	if (triNorm[1]== 0)
	{
		if (triNorm[0]>0)
			phiNorm = M_PI/2;
		else
			phiNorm = -M_PI/2;
	}
	else
		phiNorm = atan(triNorm[0]/triNorm[1]);

	if (triNorm[0]<0 && triNorm[1]<0)
		phiNorm += M_PI;
	else if (triNorm[1]<0)
		phiNorm -= M_PI;

    double Rz[3][3] = {{cos(phiNorm), -sin(phiNorm),0},
                   {sin(phiNorm), cos(phiNorm), 0},
                   {0, 0, 1}};

    double Rx[3][3] = {{1, 0, 0},
                   {0, cos(thetaNorm), -sin(thetaNorm)},
                   {0, sin(thetaNorm), cos(thetaNorm)}};

    double rotationMatrix[3][3];
    MatrixMultiplyMatrix(Rx, Rz, rotationMatrix);
    double rotationMatrix_inv[3][3];
    MatrixInverse(rotationMatrix, 3, rotationMatrix_inv);

    double newPoints[3][3];

    for (i=0; i< 3; i++)
    	MatrixMultiplyVector(rotationMatrix, 3, points[i], newPoints[i]);  //triangle is now rotated

    double zOffset = newPoints[0][2];
    for (i=0; i< 3; i++)
    	newPoints[i][2]=0.0;

    double newObsPoint[3];
    MatrixMultiplyVector(rotationMatrix, 3, obsPoint, newObsPoint);       //rotate observation point
    double zObsPoint = newObsPoint[2];
    newObsPoint[2]=0;

    double triCentre[3] = {newPoints[0][0]+newPoints[1][0]+newPoints[2][0],
    						newPoints[0][1]+newPoints[1][1]+newPoints[2][1],
							newPoints[0][2]+newPoints[1][2]+newPoints[2][2]};
    MultiplyVectorByConstant(triCentre, 3, 1.0/3);							//find centre of triangle

	for (i=0; i<3; i++)          //split triangle into three at the projection of the singularity
	{
		double AB[3];
		VectorSubtract(newPoints[i], newPoints[(i+1)%3], AB);

		double AtoCentre[3];
		VectorSubtract(newPoints[i], triCentre, AtoCentre);

		double AtoObs[3];
		VectorSubtract(newPoints[i], newObsPoint, AtoObs);

		double vec1[3], vec2[3];
		VectorCross(AB, AtoCentre, vec1);
		VectorCross(AB, AtoObs, vec2);

		int newNumInt = numIntPoints;
		if (newNumInt>10)
			newNumInt=10;

		double **returnPoints2D = NULL;
		returnPoints2D = mxMalloc(sizeof(double *) * (int)pow(newNumInt,2));
			for (j=0; j < pow(newNumInt, 2); j++)
				returnPoints2D[j] = (double *)mxMalloc(sizeof(double) * 4);

		RAR1S_2D(newPoints[i], newPoints[(i+1)%3], newObsPoint, zObsPoint-zOffset, newNumInt, returnPoints2D);  //get integration points

		int posTri = 1;
		if (vec1[2]*vec2[2] < 0)
			posTri = -1;

		newObsPoint[2] += zObsPoint;
		for (j=0; j < pow(newNumInt,2); j++)
		{
			returnPoints2D[j][2] += zOffset;
			double R = Distance(returnPoints2D[j], newObsPoint);

			MatrixMultiplyVector(rotationMatrix_inv, 3, returnPoints2D[j], returnPoints[(int)pow(newNumInt,2)*i + j]);
			returnPoints[(int)pow(newNumInt,2)*i + j][3] = posTri*R*returnPoints2D[j][3];
		}
		newObsPoint[2] -= zObsPoint;
		for (j=0; j < pow(newNumInt,2); j++)
			mxFree(returnPoints2D[j]);
		mxFree(returnPoints2D);
	}
}


#endif
