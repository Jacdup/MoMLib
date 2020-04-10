/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_CalcExciteVecPlaneWave_api.h
 *
 * Code generation for function '_coder_CalcExciteVecPlaneWave_api'
 *
 */

#ifndef _CODER_CALCEXCITEVECPLANEWAVE_API_H
#define _CODER_CALCEXCITEVECPLANEWAVE_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_CalcExciteVecPlaneWave_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void CalcExciteVecPlaneWave(real_T tri_dofs[4800], real_T node_coords[540],
  real_T observer_map[460], real_T k0, real_T E_scalfac, real_T theta_inc,
  real_T phi_inc, real_T eta_pol, creal_T V_vec[211600]);
extern void CalcExciteVecPlaneWave_api(const mxArray * const prhs[8], int32_T
  nlhs, const mxArray *plhs[1]);
extern void CalcExciteVecPlaneWave_atexit(void);
extern void CalcExciteVecPlaneWave_initialize(void);
extern void CalcExciteVecPlaneWave_terminate(void);
extern void CalcExciteVecPlaneWave_xil_shutdown(void);
extern void CalcExciteVecPlaneWave_xil_terminate(void);

#endif

/* End of code generation (_coder_CalcExciteVecPlaneWave_api.h) */
