/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CalcExciteVecPlaneWave.h
 *
 * Code generation for function 'CalcExciteVecPlaneWave'
 *
 */

#ifndef CALCEXCITEVECPLANEWAVE_H
#define CALCEXCITEVECPLANEWAVE_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "CalcExciteVecPlaneWave_types.h"

/* Function Declarations */
extern void CalcExciteVecPlaneWave(const double tri_dofs[4800], const double
  node_coords[540], const double observer_map[460], double k0, double E_scalfac,
  double theta_inc, double phi_inc, double eta_pol, creal_T V_vec[211600]);

#endif

/* End of code generation (CalcExciteVecPlaneWave.h) */
