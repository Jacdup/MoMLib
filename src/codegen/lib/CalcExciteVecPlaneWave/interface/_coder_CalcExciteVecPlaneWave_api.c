/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_CalcExciteVecPlaneWave_api.c
 *
 * Code generation for function '_coder_CalcExciteVecPlaneWave_api'
 *
 */

/* Include files */
#include "tmwtypes.h"
#include "_coder_CalcExciteVecPlaneWave_api.h"
#include "_coder_CalcExciteVecPlaneWave_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131482U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "CalcExciteVecPlaneWave",            /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4800];
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *node_coords, const char_T *identifier))[540];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[540];
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *observer_map, const char_T *identifier))[460];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *tri_dofs,
  const char_T *identifier))[4800];
static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const creal_T u
  [211600]);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[460];
static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *k0, const
  char_T *identifier);
static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4800];
static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[540];
static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[460];
static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4800]
{
  real_T (*y)[4800];
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *node_coords, const char_T *identifier))[540]
{
  real_T (*y)[540];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(node_coords), &thisId);
  emlrtDestroyArray(&node_coords);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[540]
{
  real_T (*y)[540];
  y = j_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *observer_map, const char_T *identifier))[460]
{
  real_T (*y)[460];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(observer_map), &thisId);
  emlrtDestroyArray(&observer_map);
  return y;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *tri_dofs,
  const char_T *identifier))[4800]
{
  real_T (*y)[4800];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(tri_dofs), &thisId);
  emlrtDestroyArray(&tri_dofs);
  return y;
}
  static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const creal_T u
  [211600])
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[2] = { 460, 460 };

  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxCOMPLEX);
  emlrtExportNumericArrayR2013b(sp, m0, (void *)&u[0], 8);
  emlrtAssign(&y, m0);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[460]
{
  real_T (*y)[460];
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *k0,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(k0), &thisId);
  emlrtDestroyArray(&k0);
  return y;
}

static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = l_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4800]
{
  real_T (*ret)[4800];
  static const int32_T dims[2] = { 320, 15 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[4800])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[540]
{
  real_T (*ret)[540];
  static const int32_T dims[2] = { 180, 3 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[540])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[460]
{
  real_T (*ret)[460];
  static const int32_T dims[2] = { 1, 460 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[460])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void CalcExciteVecPlaneWave_api(const mxArray * const prhs[8], int32_T nlhs,
  const mxArray *plhs[1])
{
  real_T (*tri_dofs)[4800];
  real_T (*node_coords)[540];
  real_T (*observer_map)[460];
  real_T k0;
  real_T E_scalfac;
  real_T theta_inc;
  real_T phi_inc;
  real_T eta_pol;
  static creal_T V_vec[211600];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  (void)nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  tri_dofs = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "tri_dofs");
  node_coords = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "node_coords");
  observer_map = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "observer_map");
  k0 = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "k0");
  E_scalfac = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "E_scalfac");
  theta_inc = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "theta_inc");
  phi_inc = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "phi_inc");
  eta_pol = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "eta_pol");

  /* Invoke the target function */
  CalcExciteVecPlaneWave(*tri_dofs, *node_coords, *observer_map, k0, E_scalfac,
    theta_inc, phi_inc, eta_pol, V_vec);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(&st, V_vec);
}

void CalcExciteVecPlaneWave_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  CalcExciteVecPlaneWave_xil_terminate();
  CalcExciteVecPlaneWave_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void CalcExciteVecPlaneWave_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void CalcExciteVecPlaneWave_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_CalcExciteVecPlaneWave_api.c) */
