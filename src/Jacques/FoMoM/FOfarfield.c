/*=========================================================
 * Name        : mom.c
 * Author      : Robey Beswick
 * Version     : 3.3 (15/07/2019)
 * Copyright   : 
 * Description : version 3.3  
 *               Is similar to 3.2 but is a lot more optimized. 
 
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
#include <string.h>

#include "matrix_vector_functions.h"
#include "FO_basis_function.h"
#include "triangle_functions.h"


/*
   return the integration points and weights based on a gaussian rule
 */
void GaussianQuadrature(double points[][3], int numIntPoints, double **returnPoints)
{
    /*
    Computes quadrature points on a triangle defined by points.
    Returns points and weight.
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
			//for(ii=0; ii<3; ii++) // x,y,z
			//{//TODO: only return parent coords u,v
				//returnPoints[place][ii] = (*qpPointer)[i][j]*points[0][ii] + (*qpPointer)[i][(j+1)%3]*points[1][ii] +
											//(*qpPointer)[i][(j+2)%3]*points[2][ii];
                returnPoints[place][0] = (*qpPointer)[i][j]; //Zeta_1
                returnPoints[place][1] = (*qpPointer)[i][(j+1)%3]; // Zeta_2
                returnPoints[place][2] = (*qpPointer)[i][(j+2)%3]; // Zeta_3
                //mexPrintf("returnPoints = %f %f %f\n", returnPoints[place][0],returnPoints[place][1],returnPoints[place][2]);
			//}
			//returnPoints[place][3] = 0.5*J_det*(*qpPointer)[i][3];
                returnPoints[place][3] = Area*(*qpPointer)[i][3];
               // returnPoints[place][3] = (*qpPointer)[i][3];
			place++;
		}
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
	  it then returns a scalar value that is the component of pRho in the phi direction*/
	double val =  -vec[0]*sin(obsPoint[2]) + vec[1]*cos(obsPoint[2]);
	return val;
}
double Find_pRho_theta(double vec[3], double obsPoint[3])
{
	/*recieves the cartesian pRho and observer point
	  it then returns a scalar value that is the component of pRho in the theta direction*/
	double val =  vec[0]*cos(obsPoint[1])*cos(obsPoint[2]) + vec[1]*cos(obsPoint[1])*sin(obsPoint[2]) - vec[2]*sin(obsPoint[1]);
	return val;
}


//farfield canculator
void farfield(double freq, int P, int T, double *points, double  *triangles, int amnt_E_points, complex double *I_M, double *theta_phi_Etheta_Ephi, int order)   
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
	int NumIntPoints = 3;
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
			double pPoints[3][3] = {{points[P*0+(int)(triangles[T*0 + i]-1)],points[P*1+(int)(triangles[T*0 + i]-1)],points[P*2+(int)(triangles[T*0 + i]-1)]},
		                        {points[P*0+(int)(triangles[T*1 + i]-1)],points[P*1+(int)(triangles[T*1 + i]-1)],points[P*2+(int)(triangles[T*1 + i]-1)]},
		                        {points[P*0+(int)(triangles[T*2 + i]-1)],points[P*1+(int)(triangles[T*2 + i]-1)],points[P*2+(int)(triangles[T*2 + i]-1)]}};
			
			//mexPrintf("%f %f %f\n%f %f %f\n%f %f %f\n ",pPoints[0][0],pPoints[0][1],pPoints[0][2], pPoints[1][0],pPoints[1][1],pPoints[1][2],pPoints[2][0],pPoints[2][1],pPoints[2][2]);
			//allocate the memory for the integration points
			IntPoints = mxMalloc(sizeof(double *) * NumIntPoints);
            
			for (iter=0; iter < NumIntPoints; iter++)
							IntPoints[iter] = mxMalloc(sizeof(double) * 4);
			//get and return the gaussian quadrature points for the triangle
			GaussianQuadrature(pPoints, NumIntPoints, IntPoints);
			//for(i= 0;i<NumIntPoints;i++)
				//mexPrintf("IntPoints %d = %f %f %f %f\n",i, IntPoints[i][0], IntPoints[i][1], IntPoints[i][2], IntPoints[i][3]);
			//Outer triangle 
// 			if (order == 0){
//                   int pTri[9] = {(int)(triangles[T*0 + i]-1), (int)(triangles[T*1 + i]-1), (int)(triangles[T*2 + i]-1),(int)(triangles[T*3 + i]), (int)(triangles[T*4 + i]),
//                            (int)(triangles[T*5 + i]), (int)(triangles[T*6 + i]), (int)(triangles[T*7 + i]),(int)(triangles[T*8 + i])};
// //                   int pbasisindex[3] = {pTri[6],pTri[7],pTri[8]};
//            }else{
                  int pTri[15] = {(int)(triangles[T*0 + i]-1), (int)(triangles[T*1 + i]-1), (int)(triangles[T*2 + i]-1),(int)(triangles[T*3 + i]), (int)(triangles[T*4 + i]),
                            (int)(triangles[T*5 + i]), (int)(triangles[T*6 + i]), (int)(triangles[T*7 + i]),(int)(triangles[T*8 + i]),
                            (int)(triangles[T*9 + i]),(int)(triangles[T*10 + i]),(int)(triangles[T*11 + i]),(int)(triangles[T*12 + i]),
                            (int)(triangles[T*13 + i]),(int)(triangles[T*14 + i])};
//                   int pbasisindex[6] = {pTri[6],pTri[7],pTri[8],pTri[12], pTri[13], pTri[14]};
//         }
            if ((abs(pTri[10]) != 1) && (abs(pTri[10]) != 2) ){ // Checks if there are extra columns. If not, solve standard RWG.
                order = 0;
            }else{
                order = 1;
            }
			double pArea = TriangleArea(pPoints);
			//mexPrintf("area %d = %f\n", i,pArea);
			

			//new method
			complex double E_per_tri_theta = {0.0+0.0*I};
			complex double E_per_tri_phi = {0.0+0.0*I};
            
            int pEdgeIndex;
			//now I have a triangle,pick an edge and see if the edge is a basis function edge 
			for (ii=0; ii<(3+(3*order)); ii++) //Edges of observation triangle
		    {	
                
				//The basis function number is the solution in this case
		        if (ii<3){
                    pEdgeIndex = pTri[6+(ii+2)%3];
                }else{
                     pEdgeIndex = pTri[12+(ii+2)%3];
                }
		        if (pEdgeIndex + 1)
				{
					//if there is a basis function on the edge it will have an effect
					//I am the champion
					//retrieve the current coefficient value
					complex double I_m = I_M[pEdgeIndex-1];
					//mexPrintf("I_m[%d] = %f j%f", pEdgeIndex, creal(I_M[pEdgeIndex]), cimag(I_M[pEdgeIndex]));
					// find the length of the outer edge in question
// 		            double lp = Distance(pPoints[ii], pPoints[(ii+1)%3]);
					//mexPrintf("selected edge = %d %d the edges length is %f \n",pTri[ii], pTri[(ii+1)%3], lp);
					//mexPrintf("vertex = %d\n",pTri[(ii+2)%3]);
					
					//now you will iterate through the inner integration points

					//new method
					complex double theta_holder = {0.0+0.0*I};
					complex double phi_holder = {0.0+0.0*I};
					
					for (ip = 0; ip < NumIntPoints; ip++)
					{
						//find pRho in terms of cartesian coordinates
						double pRho_cart[3];
		                //VectorSubtract(pPoints[(ii+2)%3], IntPoints[ip], pRho_cart);
		               // MultiplyVectorByConstant(pRho_cart, 3, pTri[3+(ii+2)%3]);
						//mexPrintf("%d\n", pTri[3+(ii+2)%3]);
						//mexPrintf("pRho %d = %f %f %f\n",ip, pRho_cart[0], pRho_cart[1], pRho_cart[2] );
                        double ruv[3];
                        int iiSignSelector = 0;
                        if (ii>=3){
                            iiSignSelector = 6;   
                        }
                        
                        double N = BF(IntPoints[ip],pPoints, (ii+2)%3,pTri[3+iiSignSelector+(ii+2)%3], pRho_cart);

						ParametricMap(pPoints, IntPoints[ip], ruv);
						
						///do dot product without conversion
						
						double r_dot_r_hat = (ruv[0]*cos(r_hat[2]) + ruv[1]*sin(r_hat[2]))*sin(r_hat[1]) + ruv[2]*cos(r_hat[1]);
						//printf("r_dot_r_hat = %f \n", r_dot_r_hat);

						//now new edit we find the i_theta and I_phi components of the fields first 
						double pRho_theta_comp = Find_pRho_theta(pRho_cart, r_hat);
						double pRho_phi_comp = Find_pRho_phi(pRho_cart, r_hat);

						
						//now I have the theta and phi components of pRho
						//must be multiplied by weight and e component and added to the holder
						complex double temp_theta = (IntPoints[ip][3]/(1))*pRho_theta_comp*cexp(I*k*r_dot_r_hat);
						complex double temp_phi = (IntPoints[ip][3]/(1))*pRho_phi_comp*cexp(I*k*r_dot_r_hat);
						
						theta_holder = theta_holder + ((-I*w*mu*I_m*temp_theta)/(4*M_PI));
						phi_holder = phi_holder + ((-I*w*mu*I_m*temp_phi)/(4*M_PI));
						
	 
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
 ========================================Main function================================================
 ====================================================================================================*/
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
	clock_t start, end;
    double cpu_time_used;
    
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
    int order = mxGetScalar(prhs[7]);
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
    farfield(freq, P, T, points, triangles, amnt_E_points, I_M, theta_phi_Etheta_Ephi, order); 
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


