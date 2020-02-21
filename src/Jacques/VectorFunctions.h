#include <math.h>

double VectorSize(double a[])
{
    double temp;
    temp = sqrt(pow(a[0],2) + pow(a[1],2) + pow(a[2],2));
    if (isnan(temp))
        return 0.0;
    else
        return temp;
}

void VectorCross(double u[], double v[], double x[])
{
    x[0] = u[1]*v[2] - u[2]*v[1];
    x[1] = -u[0]*v[2] + u[2]*v[0];
    x[2] = u[0]*v[1] - u[1]*v[0];
}
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
double dotProduct(double a[3], double b[3]){
    
    int i;
    double val = 0.0;
    
    for (i=0; i<3; i++){
        val += a[i]*b[i];
    }
    
    return val;
}