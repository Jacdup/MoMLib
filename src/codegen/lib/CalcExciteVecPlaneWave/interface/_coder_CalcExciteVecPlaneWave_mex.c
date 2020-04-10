/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_CalcExciteVecPlaneWave_mex.c
 *
 * Code generation for function '_coder_CalcExciteVecPlaneWave_mex'
 *
 */

/* Include files */
#include "_coder_CalcExciteVecPlaneWave_api.h"
#include "_coder_CalcExciteVecPlaneWave_mex.h"

/* Function Declarations */
static void c_CalcExciteVecPlaneWave_mexFun(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[8]);

/* Function Definitions */
static void c_CalcExciteVecPlaneWave_mexFun(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[8])
{
  const mxArray *outputs[1];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 8) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 8, 4,
                        22, "CalcExciteVecPlaneWave");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 22,
                        "CalcExciteVecPlaneWave");
  }

  /* Call the function. */
  CalcExciteVecPlaneWave_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, plhs, outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(CalcExciteVecPlaneWave_atexit);

  /* Module initialization. */
  CalcExciteVecPlaneWave_initialize();

  /* Dispatch the entry-point. */
  c_CalcExciteVecPlaneWave_mexFun(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  CalcExciteVecPlaneWave_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_CalcExciteVecPlaneWave_mex.c) */
