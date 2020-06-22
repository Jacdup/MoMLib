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
double MatrixDeterminant(double **a, int n)
{
    int i,j,j1,j2;
    double det = 0;
    double **m = NULL;
    
    if (n < 1) { /* Error */
        
    } else if (n == 1) { /* Shouldn't get used */
        det = a[0][0];
    } else if (n == 2) {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
        det = 0;
        for (j1=0;j1<n;j1++) {
            m = mxMalloc((n-1)*sizeof(double *));
            for (i=0;i<n-1;i++)
                m[i] = mxMalloc((n-1)*sizeof(double));
            for (i=1;i<n;i++) {
                j2 = 0;
                for (j=0;j<n;j++) {
                    if (j == j1)
                        continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                }
            }
            det += pow(-1.0, j1+2.0) * a[0][j1] * MatrixDeterminant(m,n-1);
            for (i=0;i<n-1;i++)
                mxFree(m[i]);
            mxFree(m);
        }
    }
    return(det);
}

/*
 * Recursive definition of determinate using expansion by minors.
 */
complex double ComplexDeterminant(complex double **a, int n)
{
    int i,j,j1,j2;
    complex double det = 0;
    complex double **m = NULL;
    
    if (n < 1)
    { /* Error */
        printf("Huge Error here!!");
    }
    else if (n == 1)
    { /* Shouldn't get used */
        det = a[0][0];
    }
    else if (n == 2)
    {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    }
    else
    {
        det = 0;
        for (j1=0;j1<n;j1++)
        {
            m = mxMalloc((n-1)*sizeof(complex double *));
            for (i=0;i<n-1;i++)
                m[i] = mxMalloc((n-1)*sizeof(complex double));
            for (i=1;i<n;i++)
            {
                j2 = 0;
                for (j=0;j<n;j++)
                {
                    if (j == j1)
                        continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                }
            }
            det += 0.1 + I*0.2;//pow(-1.0, j1+2.0) * a[0][j1] * ComplexDeterminant(m,n-1);
            for (i=0;i<n-1;i++)
                mxFree(m[i]);
            mxFree(m);
        }
    }
    return(det);
}

void MatrixCoFactor(double **a, int n, double **b)
{
    int i,j,ii,jj,i1,j1;
    double det;
    double **c;
    
    c = mxMalloc((n-1)*sizeof(double *));
    for (i=0;i<n-1;i++)
        c[i] = mxMalloc((n-1)*sizeof(double));
    
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) {
            
            /* Form the adjoint a_ij */
            i1 = 0;
            for (ii=0;ii<n;ii++) {
                if (ii == i)
                    continue;
                j1 = 0;
                for (jj=0;jj<n;jj++) {
                    if (jj == j)
                        continue;
                    c[i1][j1] = a[ii][jj];
                    j1++;
                }
                i1++;
            }
            
            /* Calculate the determinate */
            det = MatrixDeterminant(c,n-1);
            
            /* Fill in the elements of the cofactor */
            b[i][j] = pow(-1.0,i+j+2.0) * det;
        }
    }
    for (i=0;i<n-1;i++)
        mxFree(c[i]);
    mxFree(c);
}

/*
 * Find the cofactor matrix of a square matrix
 */
void ComplexCoFactor(complex double **a, int n, complex double **b)
{
    int i,j,ii,jj,i1,j1;
    complex double det;
    complex double **c;
    
    c = mxMalloc((n-1)*sizeof(complex double *));
    for (i=0;i<n-1;i++)
        c[i] = mxMalloc((n-1)*sizeof(complex double));
    
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) {
            
            /* Form the adjoint a_ij */
            i1 = 0;
            for (ii=0;ii<n;ii++) {
                if (ii == i)
                    continue;
                j1 = 0;
                for (jj=0;jj<n;jj++) {
                    if (jj == j)
                        continue;
                    c[i1][j1] = a[ii][jj];
                    j1++;
                }
                i1++;
            }
            
            /* Calculate the determinate */
            det = ComplexDeterminant(c,n-1);
            
            /* Fill in the elements of the cofactor */
            b[i][j] = pow(-1.0,i+j+2.0) * det;
        }
    }
    for (i=0;i<n-1;i++)
        mxFree(c[i]);
    mxFree(c);
}

void MatrixTranspose(double **a, int n)
{
    int i,j;
    double tmp;
    
    for (i=1;i<n;i++) {
        for (j=0;j<i;j++) {
            tmp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = tmp;
        }
    }
}

/*
 * ComplexTranspose of a square matrix, do it in place
 */
void ComplexTranspose(complex double **a, int n)
{
    int i,j;
    complex double tmp;
    
    for (i=1;i<n;i++) {
        for (j=0;j<i;j++) {
            tmp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = tmp;
        }
    }
}

/*
 * Inverse of a square matrix.
 */
void MatrixInverse(double a_in[3][3], int n, double b_in[3][3])
{
    int i,j;
    double **a, **b;
    
    a = mxMalloc(n*sizeof(double *));
    b = mxMalloc(n*sizeof(double *));
    for (i=0;i<n;i++)
    {
        a[i] = mxMalloc(n*sizeof(double));
        b[i] = mxMalloc(n*sizeof(double));
        for (j=0; j<n; j++)
        {
            a[i][j] = a_in[i][j];
            b[i][j] = b_in[i][j];
        }
    }
    
    
    
    double detA = MatrixDeterminant(a, n);
    
    MatrixCoFactor(a, n, b);
    MatrixTranspose(b, n);
    
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            b[i][j] /= detA;
    }
    if (n==1)
        b[0][0] = 1/detA;
    
    for (i=0;i<n;i++)
        for (j=0; j<n; j++)
            b_in[i][j] = b[i][j];
    
    for (i=0;i<n;i++)
    {
        mxFree(a[i]);
        mxFree(b[i]);
    }
    mxFree(a);
    mxFree(b);
    
}

/*
 * Complex Inverse of a square matrix.
 */
void ComplexInverse(complex double **a, int n, complex double **b)
{
    int i,j;
    complex double detA = ComplexDeterminant(a, n);
    
    ComplexCoFactor(a, n, b);
    ComplexTranspose(b, n);
    
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            b[i][j] /= detA;
    }
    if (n==1)
        b[0][0] = 1/detA;
}

/*
 * Multiply square matrix a with vector b
 */
void MatrixMultiplyVector(double a[][3], int n, double b[], double c[])
{
    int i,j;
    for(i=0;i<n;i++)
    {
        c[i]=0;
        for(j=0;j<n;j++)
        {
            c[i] += a[i][j]*b[j];
        }
    }
}

/*
 * Multiply square matrix a with vector b
 */
void ComplexMultiplyVector(complex double **a, int n, complex double *b, complex double *c)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        c[i]=0;
        for(j=0;j<n;j++)
        {
            c[i] += a[i][j]*b[j];
        }
    }
}

/*
 * Multiply square matrix a with square matrix b
 */
void MatrixMultiplyMatrix(double a[3][3], double b[3][3], double c[3][3])
{
    int i,j,k;
    
    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            c[i][j]=0;
    
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            for(k=0; k<3; k++)
                c[i][j] += a[i][k]*b[k][j];
}

/*
 * Solve the matrix equation Ax = b, with A, b known.
 */
void SolveLinearSystem(complex double **a, int n, complex double *b, complex double *x)
{
    int i;
    complex double **invA;
    invA = mxMalloc(n*sizeof(complex double *));
    for (i=0; i<n; i++)
        invA[i] = mxMalloc(n*sizeof(complex double));
    
    ComplexInverse(a,n,invA);
    ComplexMultiplyVector(invA, n, b, x);
    for (i=0;i<n-1;i++)
        mxFree(invA[i]);
    mxFree(invA);
}

void VectorCross(double u[], double v[], double x[])
{
    x[0] = u[1]*v[2] - u[2]*v[1];
    x[1] = -u[0]*v[2] + u[2]*v[0];
    x[2] = u[0]*v[1] - u[1]*v[0];
}

double VectorSize(double a[])
{
    double temp;
    temp = sqrt(pow(a[0],2) + pow(a[1],2) + pow(a[2],2));
    if (isnan(temp))
        return 0.0;
    else
        return temp;
}

/*
 * Writes the vector pointing from a to be in the vector c
 */
void VectorSubtract(double a[], double b[], double c[])
{
    c[0] = b[0]-a[0];
    c[1] = b[1]-a[1];
    c[2] = b[2]-a[2];
}

double Distance(double a[], double b[])
{
    double distVec[3];
    VectorSubtract(a, b, distVec);
    return VectorSize(distVec);
}

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

void MultiplyVectorByConstant(double a[], int n, double val)
{
    int i;
    for (i=0; i<n; i++)
        a[i] *= val;
}

/*
 * return the integration points and weights based on a gaussian rule
 */
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
    
//     double **quadPoints;
//     switch (numIntPoints)
// 	{
// 		case 1:
//             qpNum=1;
//             mxMalloc(sizeof(double) * qpNum);
//             quadPoints[0] = mxMalloc(sizeof(double) *
//             quadPoints = {{2,0}};
// 			//qpPointer = &quadPoints1;
//             //double quadPoints[1][2] = {{2,0}};
//
// 			break;
// 		case 2:
// 			qpPointer = &quadPoints1;
// 			qpNum = 1;
// 			break;
// 		case 3:
// 			qpPointer = &quadPoints3;
// 			qpNum = 1;
// 			break;
// 		case 4:
// 		case 5:
// 			qpPointer = &quadPoints4;
// 			qpNum = 2;
// 			break;
// 		case 6:
// 			qpPointer = &quadPoints6;
// 			qpNum = 2;
// 			break;
// 		case 7:
// 		case 8:
// 		case 9:
// 		case 10:
// 		case 11:
// 			qpPointer = &quadPoints7;
// 			qpNum = 3;
// 			break;
// 		case 12:
// 			qpPointer = &quadPoints12;
// 			qpNum = 4;
// 			break;
// 		case 13:
// 		case 14:
// 		case 15:
// 			qpPointer = &quadPoints13;
// 			qpNum = 5;
// 			break;
// 		case 16:
// 			qpPointer = &quadPoints16;
// 			qpNum = 6;
// 			break;
// 		case 25:
// 			qpPointer = &quadPoints25;
// 			qpNum = 9;
// 			break;
// 		case 33:
// 		default:
// 			qpPointer = &quadPoints33;
// 			qpNum = 11;
// 			break;
// 	}
    
    
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
            for(ii=0; ii<3; ii++)
            {
                returnPoints[place][ii] = (*qpPointer)[i][j]*points[0][ii] + (*qpPointer)[i][(j+1)%3]*points[1][ii] +
                        (*qpPointer)[i][(j+2)%3]*points[2][ii];
            }
            returnPoints[place][3] = Area*(*qpPointer)[i][3];
            place++;
        }
    }
}

/*
 * Functions for Radial Angular R1 Sqrt
 */
double wFromXY (double x, double y)
{
    return asinh(x/y);
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
    
    int qpNum;
    double (*qpPointer)[7][2];
    
    
    //Gaussian Quadrature rules in 1D
    
    //int qpNum;
    //double (*qpPointer)[7][2];
    
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

//-------------------------NEW ROUTINES--------------------------------

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
    
    //finds the shortest distance between the 9 vertices
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

//void MoM(double freq, int P, double **points, int T, int **triangles, int N, int **edges, complex double **Zmat, complex double *Vvec)
void MoM(double freq, int P, double *points, int T, double *triangles, int N, complex double *Zmat,int *obs_map, int *src_map, int num_obs, int num_src, int MODE)
{
    //This code choses the default number of integration points per triangle and then deviates from it with checks
    int NumOuterIntPoints = 0;
    int NumInnerIntPoints = 0;
    //defining a proximity by hand is unwise as this is dependant on triangle diameter
    //double prox = 0.33;
    
    //define the integration domain vector
    int mode_select[6][4] = {{3,3,3,6},
    {3,4,6,7},
    {4,6,7,12},
    {6,7,12,16},
    {7,13,16,25}};
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
    
    // Allocate memory outside loops
    double **OuterIntPoints = NULL;
    OuterIntPoints = mxMalloc(sizeof(double *) * 25);
    for (iter=0; iter < 25; iter++)
        OuterIntPoints[iter] = mxMalloc(sizeof(double) * 4);
    
    //clock stuff
    clock_t StartInner, EndInner;
    double Inner_loop_time_used;
    
    
    double prog = 0.0; //used to monitor progress
    
    for (i=0; i< T; i++) //Observation triangle
    {
        
        //create a 3 by 3 matrix of a single triangles points each row is the x y z coordinate
        //pPoints is my outer triangle
        //new reference of pPoints
        double pPoints[3][3] = {{points[P*0+(int)(triangles[T*0 + i]-1)],points[P*1+(int)(triangles[T*0 + i]-1)],points[P*2+(int)(triangles[T*0 + i]-1)]},
        {points[P*0+(int)(triangles[T*1 + i]-1)],points[P*1+(int)(triangles[T*1 + i]-1)],points[P*2+(int)(triangles[T*1 + i]-1)]},
        {points[P*0+(int)(triangles[T*2 + i]-1)],points[P*1+(int)(triangles[T*2 + i]-1)],points[P*2+(int)(triangles[T*2 + i]-1)]}};
        //Outer triangle
        int pTri[9] = {(int)(triangles[T*0 + i]-1), (int)(triangles[T*1 + i]-1), (int)(triangles[T*2 + i]-1),(int)(triangles[T*3 + i]), (int)(triangles[T*4 + i]),
        (int)(triangles[T*5 + i]), (int)(triangles[T*6 + i]), (int)(triangles[T*7 + i]),(int)(triangles[T*8 + i])};
        
        int pbasisindex[3] = {pTri[6],pTri[7],pTri[8]};
        if (ObserverOnTriangle(pbasisindex, obs_map))
        {
            double pArea = TriangleArea(pPoints);

            //double **OuterIntPoints = NULL;
            //StartInner = clock();
            for (j=0; j<T; j++) //Testing triangle
            {
                // qPoints is my inner triangle where I will possibly need to use RAR to integrate
                double qPoints[3][3] = {{points[P*0+(int)(triangles[T*0 + j]-1)],points[P*1+(int)(triangles[T*0 + j]-1)],points[P*2+(int)(triangles[T*0 + j]-1)]},
                {points[P*0+(int)(triangles[T*1 + j]-1)],points[P*1+(int)(triangles[T*1 + j]-1)],points[P*2+(int)(triangles[T*1 + j]-1)]},
                {points[P*0+(int)(triangles[T*2 + j]-1)],points[P*1+(int)(triangles[T*2 + j]-1)],points[P*2+(int)(triangles[T*2 + j]-1)]}};
                
                //inner triangle
                int qTri[9] = {(int)(triangles[T*0 + j]-1), (int)(triangles[T*1 + j]-1), (int)(triangles[T*2 + j]-1),(int)(triangles[T*3 + j]), (int)(triangles[T*4 + j]),
                (int)(triangles[T*5 + j]), (int)(triangles[T*6 + j]), (int)(triangles[T*7 + j]),(int)(triangles[T*8 + j])};
                
                //below checks if any single edge in the triangle has a basis function defined over it
                int qbasisindex[3] = {qTri[6],qTri[7],qTri[8]};
                //WIll only come in here if there is source and observer, this is the only time actual memory will be allocated for the solutions
                if (ObserverOnTriangle(qbasisindex, src_map))
                {
                    //Because we now know outer and inner integration points will be allocated memory
// 					double **OuterIntPoints = NULL;
                    //double **InnerIntPoints = NULL;
                    // Find the area and centrepoint of the triangle
                    double qArea = TriangleArea(qPoints);
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //-----------------------------------------------------First round of checks--------------------------------------------------------
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //----------------------------------------------------------------------------------------------------------------------------------
                    //these are local variables, they change per triangle pair, will only be used
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
                    
 
                    //find triangle centroids and distance between them
                    TriangleCentre(pPoints,pCentroid);
                    TriangleCentre(qPoints,qCentroid);
                    d_centrepoints = Distance(pCentroid,qCentroid);
 
                    //find the triangle diameters and the maximum of them
                    pDiameter = TriangleDiameter(pPoints);
                    qDiameter = TriangleDiameter(qPoints);
  
                    if (pDiameter > qDiameter)
                    { max_diameter = pDiameter;
                    }else{ max_diameter = qDiameter;}
                    //shortest possible distance between triangles given the centroid
                    d_conservative = d_centrepoints - (2.0/3.0)*(pDiameter + qDiameter);

                    //d conservative is the shortest possible distance between the triangles given the centroid
                    if (0 > d_conservative)
                    {
                        //if this is the case the triangles are either touching or very close,
                        //it will always be negative if the triangles are touching
                        int DontTouch = TrianglesDontTouch(pTri, qTri);//can only be used locally
                        if (!DontTouch)
                        {
                            //if we are in here it means the triangles touch
                            //instantly set outer to 16 and further_checks = 1
                            //no inner points set yet, further checks handles that
                            //set mode to level 4
                            NumOuterIntPoints = mode_select[4][MODE];
                            further_check_required = 1;
     
                        }else
                        {
                            d_9 = ApproxShortDistanceBetweenTriangles_9(pPoints,qPoints);
                            
                            
                            if(d_9 > (1.5*max_diameter))
                            {
                                NumOuterIntPoints = mode_select[1][MODE];
                                NumInnerIntPoints = mode_select[1][MODE];
                                count_6_6++;
                            }else if(d_9 > (0.75*max_diameter))
                            {
                                NumOuterIntPoints = mode_select[2][MODE];
                                NumInnerIntPoints = mode_select[2][MODE];
                                count_7_7++;
                            }else
                            {
                                //This means the triangles are very close and we will do strict test
                                //do 36 point check  
                                d_36 = ApproxShortDistanceBetweenTriangles_36(pPoints,qPoints);
                                
                                if(d_36 > 0.375*max_diameter)
                                {
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
                            NumOuterIntPoints = mode_select[0][MODE];
                            NumInnerIntPoints = mode_select[0][MODE];
                            count_4_4++;
                        }else if(d_conservative > (1.5*max_diameter))
                        {
                            NumOuterIntPoints = mode_select[1][MODE];
                            NumInnerIntPoints = mode_select[1][MODE];
                            count_6_6++;
                        }else if(d_conservative > (0.75*max_diameter))
                        {
                            
                            NumOuterIntPoints = mode_select[2][MODE];
                            NumInnerIntPoints = mode_select[2][MODE];
                            count_7_7++;
                        }else
                        {	//lower then 0.75*max_d
                            //we perform the 9 point check
                            d_9 = ApproxShortDistanceBetweenTriangles_9(pPoints,qPoints);
                            if (d_9 > (1.5*max_diameter))
                            {
                                NumOuterIntPoints = mode_select[1][MODE];
                                NumInnerIntPoints = mode_select[1][MODE];
                                count_6_6++;
                            }else if(d_9 > (0.75*max_diameter))
                            {
                                NumOuterIntPoints = mode_select[2][MODE];
                                NumInnerIntPoints = mode_select[2][MODE];
                                count_7_7++;
                            }else
                            {
                                //so if less then 0.75 multiplied by the max diameter
                                //more stringent tests applied
                                double d_36 = ApproxShortDistanceBetweenTriangles_36(pPoints,qPoints);
                                if(d_36 > 0.375*max_diameter)
                                {
                                    NumOuterIntPoints = mode_select[3][MODE];
                                    NumInnerIntPoints = mode_select[3][MODE];
                                    count_12_12++;
                                }else
                                {
                                    NumOuterIntPoints = mode_select[4][MODE];
                                    further_check_required = 1;
                                }
                                
                            }
                        }     
                    }
                    //at this point I know my integration domains and I will allocate memory
                    //allocate memory for the integration points, do outer and inner here if far interaction
                    //printf("further_check_required = %d\n",further_check_required);
                    if (!further_check_required)
                    {           
                        //far interaction (further checks aren't required)
                        OuterIntPoints = mxMalloc(sizeof(double *) * NumOuterIntPoints);
                        for (iter=0; iter < NumOuterIntPoints; iter++)
                            OuterIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                        //find the integration points for the outer triangle
                        GaussianQuadrature(pPoints, NumOuterIntPoints, OuterIntPoints);
                        
//                         //mxMalloc inner integration points
//                         //double **InnerIntPoints;
//                         newNumInt = NumInnerIntPoints;
//                         InnerIntPoints = mxMalloc(sizeof(double *) * newNumInt);
//                         for (iter=0; iter < newNumInt; iter++)
//                             InnerIntPoints[iter] = mxMalloc(sizeof(double) * 4);
//
//
//                         //find the integration points for the inner triangle
//                         GaussianQuadrature(qPoints, NumInnerIntPoints, InnerIntPoints);
                        
                        
                    }else
                    {
                        //This is a near interaction a check must now be done per outer integration point.
                        //set the outer integration points to a high value(level 4), I will play with this later
                        OuterIntPoints = mxMalloc(sizeof(double *) * NumOuterIntPoints);
                        for (iter=0; iter < NumOuterIntPoints; iter++)
                            OuterIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                        //find the integration points for the outer triangle
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
                            newNumInt = 3*(int)pow(NumInnerIntPoints,2);
                            if (newNumInt > 300)
                                newNumInt = 300;
                            
                            //allocate the memory for OuterIntPoints
                            InnerIntPoints = mxMalloc(sizeof(double *)*newNumInt);
                            for (iter=0; iter < newNumInt; iter++)
                                InnerIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                            
                            //so the inner triangle is now integrated using radial angular methods
                            //Note that RAR1S also recieves the outer integration point and original NumInnerIntPoints
                            RAR1S(qPoints, OuterIntPoints[oip], NumInnerIntPoints, InnerIntPoints);
                            count_16_RAR++;
                            
                        }else if(further_check_required)
                        {
                            //check accurate proximity from oip to inner triangle
                            double proximity = StrictDistanceToTriangle(qPoints, OuterIntPoints[oip]);
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
                                
                            }else
                            {
                                //if we come in here outer is already set to level 4
                                
                                //level 4 outside and RAR inside
                                //we will use radial angular for the inner integral
                                
                                //NumInnerIntPoints is recieved by RAR, this can be set higher max of 10 is set
                                NumInnerIntPoints = RAR_level[MODE];
                                newNumInt = 3*(int)pow(NumInnerIntPoints,2);
                                if (newNumInt > 300)
                                    newNumInt = 300;
                                
                                //allocate the memory for InnerIntPoints
                                InnerIntPoints = mxMalloc(sizeof(double *)*newNumInt);
                                for (iter=0; iter < newNumInt; iter++)
                                    InnerIntPoints[iter] = mxMalloc(sizeof(double) * 4);
                                
                                //so the inner triangle is now integrated using radial angular methods
                                //Note that RAR1S also recieves the outer integration point and original NumInnerIntPoints
                                RAR1S(qPoints, OuterIntPoints[oip], NumInnerIntPoints, InnerIntPoints);
                                count_16_RAR++;
                            }
                            
                            
                        }else if(!further_check_required  && (i != j))
                        {
                            newNumInt = NumInnerIntPoints;
                            InnerIntPoints = mxMalloc(sizeof(double *) * NumInnerIntPoints);
                            for (iter=0; iter < NumInnerIntPoints; iter++)
                                InnerIntPoints[iter] = mxMalloc(sizeof(double) * 4);   
                            
                            //find the integration points for the inner triangle
                            GaussianQuadrature(qPoints, NumInnerIntPoints, InnerIntPoints);
                        }
                        //----------------------------------------------------------------------------------------------------------------------------------
                        //----------------------------------------------------------------------------------------------------------------------------------
                        //------------------------------------------------second round of checks complete---------------------------------------------------
                        //----------------------------------------------------------------------------------------------------------------------------------
                        //----------------------------------------------------------------------------------------------------------------------------------
                           for (ii=0; ii<3; ii++) //Edges of observation triangle
                        {

                            int pEdgeIndex = pTri[6+(ii+2)%3];
                            int p_obs_index = obs_map[(pTri[6+(ii+2)%3]-1)];
                            
                            // below just checks that the edge index is any number other the -1, meaning it will be a DOF
                            // this if statement asks two questions, is it a basis function? and is it an observer?
                            if ((pEdgeIndex + 1) && (p_obs_index + 1))
                            {
                                double lp = Distance(pPoints[ii], pPoints[(ii+1)%3]);
                                
                                double pRho[3];
                                VectorSubtract(pPoints[(ii+2)%3], OuterIntPoints[oip], pRho);
                                MultiplyVectorByConstant(pRho, 3, pTri[3+(ii+2)%3]);

                                for (jj=0; jj<3; jj++) //Edges of testing triangle
                                {

                                    int qEdgeIndex = qTri[6+(jj+2)%3];
                                    int q_src_index = src_map[(qTri[6+(jj+2)%3]-1)];
                                    
                                    if ((qEdgeIndex + 1) && (q_src_index + 1))
                                    {
                                        
                                        complex double A[3] = {0.0, 0.0, 0.0};
                                        complex double Phi = 0+0*I;
                                        
                                        //inner integral
                                        for (iip=0; iip < newNumInt; iip++)
                                        {
                                            double qRho[3];
                                            VectorSubtract(qPoints[(jj+2)%3], InnerIntPoints[iip], qRho);
                                            MultiplyVectorByConstant( qRho, 3, qTri[3+(jj+2)%3]);
                                            
                                            double Rm = Distance(OuterIntPoints[oip], InnerIntPoints[iip]);
                                            double lq = Distance(qPoints[jj], qPoints[(jj+1)%3]);
                                            
                                            complex double temp = (InnerIntPoints[iip][3]/(Rm*qArea))*cexp(-I*k*Rm);
                                            
//                                             if (isnan(real(temp))){
//                                                 mexPrintf( "in here");
//                                             }

                                            A[0] += (qRho[0]*I*lq*mu*w*temp)/(16*M_PI);
                                            A[1] += (qRho[1]*I*lq*mu*w*temp)/(16*M_PI);
                                            A[2] += (qRho[2]*I*lq*mu*w*temp)/(16*M_PI);
                                            Phi += qTri[3+(jj+2)%3]*temp*lq/(4*I*M_PI*w*eps);
                                        }
                                        complex double AdotpRho = A[0]*pRho[0] + A[1]*pRho[1] + A[2]*pRho[2];
                                        
                                        Zmat[(num_obs*(q_src_index-1)) + (p_obs_index-1)] += (OuterIntPoints[oip][3]/pArea)*lp*(AdotpRho + pTri[3+(ii+2)%3]*Phi);
                                        
                                        
                                    }
                                }
                                   
                            }
                        }//ii
                        for (iter=0; iter < newNumInt; iter++)
                                    mxFree(InnerIntPoints[iter]);
                                mxFree(InnerIntPoints);
                     
                    }//oip
                    
                    
                    for (iter=0; iter < NumOuterIntPoints; iter++)
                        mxFree(OuterIntPoints[iter]);
                    mxFree(OuterIntPoints);
                }//the brace is new inner if statement
            }
            //EndInner = clock();
            //Inner_loop_time_used = ((double) (EndInner - StartInner)) / CLOCKS_PER_SEC;
            //mexPrintf("inner loop took %f seconds to execute \n", Inner_loop_time_used);
            
            
            
            //for (iter=0; iter < NumOuterIntPoints; iter++)
            //    mxFree(OuterIntPoints[iter]);
            //mxFree(OuterIntPoints);
            
        }//this brace is new if statement
        
        
    }
    //prints the occurences of how many times each domain was used
    mexPrintf("\n%d %d: %d\n%d %d: %d\n%d %d: %d\n%d %d: %d\n%d %d: %d\n%d RAR: %d\n",mode_select[0][MODE],mode_select[0][MODE] ,count_4_4, mode_select[1][MODE],mode_select[1][MODE],count_6_6, mode_select[2][MODE],mode_select[2][MODE],count_7_7, mode_select[3][MODE],mode_select[3][MODE],count_12_12, mode_select[4][MODE],mode_select[4][MODE], count_16_16,mode_select[4][MODE], count_16_RAR);
// 	for (iter=0; iter < NumOuterIntPoints; iter++)
// 						mxFree(OuterIntPoints[iter]);
// 					mxFree(OuterIntPoints);
    
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
    
//     //clock stuff
// 	clock_t start, end;
//     double cpu_time_used;
    
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
    //for(iter = 0; iter < P;iter++)
    // mexPrintf("%lf %lf %lf\n",points[P*0 + iter], points[P*1 + iter], points[P*2 + iter]);
    
    
    /*Populate triangles array*/
    T = mxGetM(prhs[1]);
    // mexPrintf("%d \n",T);
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
//     Zmat = mxMalloc(sizeof(complex double) * num_obs * num_src);//the entire memory allocation of
    plhs[0]= mxCreateDoubleMatrix(num_obs, num_src, mxCOMPLEX);
    Zmat = mxGetComplexDoubles(plhs[0]);
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
    
    //trying to fix the data free stuff
    
    
    //you don't have to mxFree anything as the prhs and plhs are mxArray types which matlab
    //handles itself
    
//     mxSetComplexDoubles(plhs[0], Zmat);
    
//     mxFree(Zmat);
    
    mexPrintf("\nIf it were done when 'tis done, then 'twere well It were done quickly\n");
    
    
}


