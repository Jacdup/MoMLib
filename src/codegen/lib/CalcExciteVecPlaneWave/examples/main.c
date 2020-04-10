/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "CalcExciteVecPlaneWave.h"
#include "main.h"
#include "CalcExciteVecPlaneWave_terminate.h"
#include "CalcExciteVecPlaneWave_initialize.h"

/* Function Declarations */
static void argInit_180x3_real_T(double result[540]);
static void argInit_1x460_real_T(double result[460]);
static void argInit_320x15_real_T(double result[4800]);
static double argInit_real_T(void);
static void main_CalcExciteVecPlaneWave(void);

/* Function Definitions */
static void argInit_180x3_real_T(double result[540])
{
  int idx0;
  double result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 180; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_real_T();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 180] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 360] = argInit_real_T();
  }
}

static void argInit_1x460_real_T(double result[460])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 460; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
  }
}

static void argInit_320x15_real_T(double result[4800])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 320; idx0++) {
    for (idx1 = 0; idx1 < 15; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 320 * idx1] = argInit_real_T();
    }
  }
}

static double argInit_real_T(void)
{
  return 0.0;
}

static void main_CalcExciteVecPlaneWave(void)
{
  double dv0[4800];
  double dv1[540];
  double dv2[460];
  static creal_T V_vec[211600];

  /* Initialize function 'CalcExciteVecPlaneWave' input arguments. */
  /* Initialize function input argument 'tri_dofs'. */
  /* Initialize function input argument 'node_coords'. */
  /* Initialize function input argument 'observer_map'. */
  /* Call the entry-point 'CalcExciteVecPlaneWave'. */
  argInit_320x15_real_T(dv0);
  argInit_180x3_real_T(dv1);
  argInit_1x460_real_T(dv2);
  CalcExciteVecPlaneWave(dv0, dv1, dv2, argInit_real_T(), argInit_real_T(),
    argInit_real_T(), argInit_real_T(), argInit_real_T(), V_vec);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  CalcExciteVecPlaneWave_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_CalcExciteVecPlaneWave();

  /* Terminate the application.
     You do not need to do this more than one time. */
  CalcExciteVecPlaneWave_terminate();
  return 0;
}

/* End of code generation (main.c) */
