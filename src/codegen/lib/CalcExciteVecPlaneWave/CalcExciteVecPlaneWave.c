/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CalcExciteVecPlaneWave.c
 *
 * Code generation for function 'CalcExciteVecPlaneWave'
 *
 */

/* Include files */
#include <math.h>
#include <string.h>
#include "CalcExciteVecPlaneWave.h"

/* Function Definitions */
void CalcExciteVecPlaneWave(const double tri_dofs[4800], const double
  node_coords[540], const double observer_map[460], double k0, double E_scalfac,
  double theta_inc, double phi_inc, double eta_pol, creal_T V_vec[211600])
{
  double qruletri[28];
  double cross_val_idx_1;
  double cross_val_idx_2;
  double unitvec_beta_idx_0;
  double r;
  double unitvec_beta_idx_1;
  double unitvec_beta_idx_2_tmp;
  double cross_val_idx_0;
  double edgevectemp_idx_2;
  double unitvec_eta_idx_0;
  double unitvec_eta_idx_1;
  double unitvec_eta_idx_2;
  int ii;
  double triidx[3];
  int edgevectemp_idx_0_tmp_tmp;
  int edgevectemp_idx_0_tmp;
  double edgevectemp_idx_0;
  int b_idx_0_tmp;
  double b_idx_0;
  double edgevectemp_idx_1;
  int i0;
  double b_node_coords[9];
  double qcoords[21];
  creal_T y[7];
  creal_T Einc_eval[21];
  static const signed char local_edge_nodes_def[6] = { 2, 1, 1, 3, 3, 2 };

  double b_qcoords[21];
  creal_T b_Einc_eval[21];

  /*  The excitation vector is cacluated for the specfied dofs. <observer_map> */
  /*  is -1 for dofs not to be tested with, with ascending 1,2,3,... for the */
  /*  entries that has to be calculated. If only M entries must be calcuated, */
  /*  then <V_vec> only has M entries. */
  /*  */
  /*  The dof and mesh arguments are described here: */
  /*  mesh_data    : output of the function <CreateMeshData.m> */
  /*  dof_data     : output of the function <CreateBasisFunctions.m> */
  /*  */
  /*  2019-12-16: Created. MMB. */
  /*  Init: */
  memset(&V_vec[0], 0, 211600U * sizeof(creal_T));

  /*  V_vec                = uint16(complex(size(ob */
  /*  Assign the 7 point rule data: [Dunavant] degree of precision 5 */
  qruletri[0] = 0.333333333333333;
  qruletri[1] = 0.05971587178977;
  qruletri[2] = 0.470142064105115;
  qruletri[3] = 0.470142064105115;
  qruletri[4] = 0.797426985353087;
  qruletri[5] = 0.101286507323456;
  qruletri[6] = 0.101286507323456;
  qruletri[7] = 0.333333333333333;
  qruletri[8] = 0.470142064105115;
  qruletri[9] = 0.470142064105115;
  qruletri[10] = 0.05971587178977;
  qruletri[11] = 0.101286507323456;
  qruletri[12] = 0.101286507323456;
  qruletri[13] = 0.797426985353087;
  qruletri[14] = 0.333333333333333;
  qruletri[15] = 0.470142064105115;
  qruletri[16] = 0.05971587178977;
  qruletri[17] = 0.470142064105115;
  qruletri[18] = 0.101286507323456;
  qruletri[19] = 0.797426985353087;
  qruletri[20] = 0.101286507323456;
  qruletri[21] = 0.225;
  qruletri[22] = 0.132394152788506;
  qruletri[23] = 0.132394152788506;
  qruletri[24] = 0.132394152788506;
  qruletri[25] = 0.125939180544827;
  qruletri[26] = 0.125939180544827;
  qruletri[27] = 0.125939180544827;
  cross_val_idx_1 = sin(theta_inc);
  cross_val_idx_2 = cos(phi_inc);
  unitvec_beta_idx_0 = -(cross_val_idx_1 * cross_val_idx_2);
  r = sin(phi_inc);
  unitvec_beta_idx_1 = -(cross_val_idx_1 * r);
  unitvec_beta_idx_2_tmp = cos(theta_inc);
  cross_val_idx_0 = -cos(eta_pol);
  edgevectemp_idx_2 = sin(eta_pol);
  unitvec_eta_idx_0 = cross_val_idx_0 * (unitvec_beta_idx_2_tmp *
    cross_val_idx_2) + edgevectemp_idx_2 * -r;
  unitvec_eta_idx_1 = cross_val_idx_0 * (unitvec_beta_idx_2_tmp * r) +
    edgevectemp_idx_2 * cross_val_idx_2;
  unitvec_eta_idx_2 = cross_val_idx_0 * -cross_val_idx_1 + edgevectemp_idx_2 *
    0.0;

  /*  Cycle over all tri_dofs (list of triangles, possibly listing some twice, */
  /*  with dofs associated with each listed triangle). Triangles are listed */
  /*  twice in case of junctions where two basis functions can be associated */
  /*  with the same edge of the same triangle.  */
  /*  Uver each triangle the relevant basis functions dot Einc are integrated */
  /*  and assembled into V-vec. */
  for (ii = 0; ii < 320; ii++) {
    /*  Make a list of the indices into <V_vec> for the three basis functions */
    /*  on this tri: */
    triidx[0] = tri_dofs[ii + 1920];
    triidx[1] = tri_dofs[ii + 2240];
    triidx[2] = tri_dofs[ii + 2560];
    if (tri_dofs[ii + 1920] > 0.0) {
      triidx[0] = observer_map[(int)tri_dofs[ii + 1920] - 1];
    }

    if (tri_dofs[ii + 2240] > 0.0) {
      triidx[1] = observer_map[(int)tri_dofs[ii + 2240] - 1];
    }

    if (tri_dofs[ii + 2560] > 0.0) {
      triidx[2] = observer_map[(int)tri_dofs[ii + 2560] - 1];
    }

    /*  Calculate the contributions from this tri to <V_vec>:     */
    /*  global nodes of current tri */
    edgevectemp_idx_0_tmp_tmp = (int)tri_dofs[ii];
    cross_val_idx_2 = node_coords[edgevectemp_idx_0_tmp_tmp - 1];
    edgevectemp_idx_0_tmp = (int)tri_dofs[320 + ii];
    edgevectemp_idx_0 = node_coords[edgevectemp_idx_0_tmp - 1] - cross_val_idx_2;
    b_idx_0_tmp = (int)tri_dofs[640 + ii];
    b_idx_0 = node_coords[b_idx_0_tmp - 1] - cross_val_idx_2;
    r = node_coords[edgevectemp_idx_0_tmp_tmp + 179];
    edgevectemp_idx_1 = node_coords[edgevectemp_idx_0_tmp + 179] - r;
    r = node_coords[b_idx_0_tmp + 179] - r;
    cross_val_idx_2 = node_coords[edgevectemp_idx_0_tmp_tmp + 359];
    edgevectemp_idx_2 = node_coords[edgevectemp_idx_0_tmp + 359] -
      cross_val_idx_2;
    cross_val_idx_2 = node_coords[b_idx_0_tmp + 359] - cross_val_idx_2;
    cross_val_idx_0 = edgevectemp_idx_1 * cross_val_idx_2 - edgevectemp_idx_2 *
      r;
    cross_val_idx_1 = edgevectemp_idx_2 * b_idx_0 - edgevectemp_idx_0 *
      cross_val_idx_2;
    cross_val_idx_2 = edgevectemp_idx_0 * r - edgevectemp_idx_1 * b_idx_0;
    edgevectemp_idx_1 = 0.5 * sqrt((cross_val_idx_0 * cross_val_idx_0 +
      cross_val_idx_1 * cross_val_idx_1) + cross_val_idx_2 * cross_val_idx_2);
    for (i0 = 0; i0 < 3; i0++) {
      b_node_coords[3 * i0] = node_coords[((int)tri_dofs[ii] + 180 * i0) - 1];
      b_node_coords[1 + 3 * i0] = node_coords[(edgevectemp_idx_0_tmp + 180 * i0)
        - 1];
      b_node_coords[2 + 3 * i0] = node_coords[(b_idx_0_tmp + 180 * i0) - 1];
    }

    for (i0 = 0; i0 < 7; i0++) {
      for (edgevectemp_idx_0_tmp_tmp = 0; edgevectemp_idx_0_tmp_tmp < 3;
           edgevectemp_idx_0_tmp_tmp++) {
        qcoords[i0 + 7 * edgevectemp_idx_0_tmp_tmp] = (qruletri[i0] *
          b_node_coords[3 * edgevectemp_idx_0_tmp_tmp] + qruletri[i0 + 7] *
          b_node_coords[1 + 3 * edgevectemp_idx_0_tmp_tmp]) + qruletri[i0 + 14] *
          b_node_coords[2 + 3 * edgevectemp_idx_0_tmp_tmp];
      }
    }

    cross_val_idx_0 = k0 * 0.0;
    for (edgevectemp_idx_0_tmp_tmp = 0; edgevectemp_idx_0_tmp_tmp < 7;
         edgevectemp_idx_0_tmp_tmp++) {
      edgevectemp_idx_2 = qcoords[edgevectemp_idx_0_tmp_tmp] * cross_val_idx_0;
      cross_val_idx_2 = qcoords[edgevectemp_idx_0_tmp_tmp] * -k0;
      r = edgevectemp_idx_2 * unitvec_beta_idx_0 - cross_val_idx_2 * 0.0;
      b_idx_0 = edgevectemp_idx_2 * 0.0 + cross_val_idx_2 * unitvec_beta_idx_0;
      cross_val_idx_2 = qcoords[edgevectemp_idx_0_tmp_tmp + 7];
      edgevectemp_idx_2 = cross_val_idx_2 * cross_val_idx_0;
      cross_val_idx_2 *= -k0;
      r += edgevectemp_idx_2 * unitvec_beta_idx_1 - cross_val_idx_2 * 0.0;
      b_idx_0 += edgevectemp_idx_2 * 0.0 + cross_val_idx_2 * unitvec_beta_idx_1;
      cross_val_idx_2 = qcoords[edgevectemp_idx_0_tmp_tmp + 14];
      edgevectemp_idx_2 = cross_val_idx_2 * cross_val_idx_0;
      cross_val_idx_2 *= -k0;
      r += edgevectemp_idx_2 * -unitvec_beta_idx_2_tmp - cross_val_idx_2 * 0.0;
      b_idx_0 += edgevectemp_idx_2 * 0.0 + cross_val_idx_2 *
        -unitvec_beta_idx_2_tmp;
      y[edgevectemp_idx_0_tmp_tmp].re = r;
      y[edgevectemp_idx_0_tmp_tmp].im = b_idx_0;
      if (b_idx_0 == 0.0) {
        y[edgevectemp_idx_0_tmp_tmp].re = exp(r);
        y[edgevectemp_idx_0_tmp_tmp].im = 0.0;
      } else {
        r = exp(r / 2.0);
        y[edgevectemp_idx_0_tmp_tmp].re = r * (r * cos(b_idx_0));
        y[edgevectemp_idx_0_tmp_tmp].im = r * (r * sin(b_idx_0));
      }
    }

    for (i0 = 0; i0 < 7; i0++) {
      r = E_scalfac * y[i0].re;
      cross_val_idx_2 = E_scalfac * y[i0].im;
      Einc_eval[i0].re = r * unitvec_eta_idx_0 - cross_val_idx_2 * 0.0;
      Einc_eval[i0].im = r * 0.0 + cross_val_idx_2 * unitvec_eta_idx_0;
      Einc_eval[i0 + 7].re = r * unitvec_eta_idx_1 - cross_val_idx_2 * 0.0;
      Einc_eval[i0 + 7].im = r * 0.0 + cross_val_idx_2 * unitvec_eta_idx_1;
      Einc_eval[i0 + 14].re = r * unitvec_eta_idx_2 - cross_val_idx_2 * 0.0;
      Einc_eval[i0 + 14].im = r * 0.0 + cross_val_idx_2 * unitvec_eta_idx_2;
    }

    /*  each row is the E-field(Ex,Ey,Ez) at each quad point */
    for (edgevectemp_idx_0_tmp = 0; edgevectemp_idx_0_tmp < 3;
         edgevectemp_idx_0_tmp++) {
      /*  cycle over the three edges of this triangle */
      if (triidx[edgevectemp_idx_0_tmp] > 0.0) {
        /*  this basis function has a contribution to <V_vec> */
        /*  Edge length: */
        /*  global nodes of current edge */
        b_idx_0_tmp = (int)tri_dofs[ii + 320 * (local_edge_nodes_def[3 +
          edgevectemp_idx_0_tmp] - 1)];
        edgevectemp_idx_0_tmp_tmp = (int)tri_dofs[ii + 320 *
          (local_edge_nodes_def[edgevectemp_idx_0_tmp] - 1)];

        /*  Eval basis function at the qpoints: */
        cross_val_idx_2 = node_coords[b_idx_0_tmp - 1] -
          node_coords[edgevectemp_idx_0_tmp_tmp - 1];
        r = cross_val_idx_2 * cross_val_idx_2;
        cross_val_idx_2 = node_coords[b_idx_0_tmp + 179] -
          node_coords[edgevectemp_idx_0_tmp_tmp + 179];
        r += cross_val_idx_2 * cross_val_idx_2;
        cross_val_idx_2 = node_coords[b_idx_0_tmp + 359] -
          node_coords[edgevectemp_idx_0_tmp_tmp + 359];
        r += cross_val_idx_2 * cross_val_idx_2;
        cross_val_idx_0 = 0.5 * (sqrt(r) / edgevectemp_idx_1) * tri_dofs[ii +
          320 * (edgevectemp_idx_0_tmp + 3)];

        /*  Calculate the integral result: */
        /*  Assemble into <V_vec>: */
        b_idx_0_tmp = (int)tri_dofs[ii + 320 * edgevectemp_idx_0_tmp];
        for (i0 = 0; i0 < 7; i0++) {
          b_qcoords[i0] = qcoords[i0] - node_coords[b_idx_0_tmp - 1];
          b_qcoords[i0 + 7] = qcoords[i0 + 7] - node_coords[b_idx_0_tmp + 179];
          b_qcoords[i0 + 14] = qcoords[i0 + 14] - node_coords[b_idx_0_tmp + 359];
        }

        for (i0 = 0; i0 < 21; i0++) {
          cross_val_idx_2 = cross_val_idx_0 * b_qcoords[i0];
          b_Einc_eval[i0].re = cross_val_idx_2 * Einc_eval[i0].re;
          b_Einc_eval[i0].im = cross_val_idx_2 * Einc_eval[i0].im;
        }

        cross_val_idx_2 = 0.0;
        cross_val_idx_1 = 0.0;
        for (i0 = 0; i0 < 7; i0++) {
          r = ((b_Einc_eval[i0].re - b_Einc_eval[i0].im * 0.0) + (b_Einc_eval[i0
                + 7].re - b_Einc_eval[i0 + 7].im * 0.0)) + (b_Einc_eval[i0 + 14]
            .re - b_Einc_eval[i0 + 14].im * 0.0);
          b_idx_0 = ((b_Einc_eval[i0].re * 0.0 + b_Einc_eval[i0].im) +
                     (b_Einc_eval[i0 + 7].re * 0.0 + b_Einc_eval[i0 + 7].im)) +
            (b_Einc_eval[i0 + 14].re * 0.0 + b_Einc_eval[i0 + 14].im);
          y[i0].re = r;
          y[i0].im = b_idx_0;
          edgevectemp_idx_2 = edgevectemp_idx_1 * qruletri[21 + i0];
          cross_val_idx_2 += edgevectemp_idx_2 * r - 0.0 * b_idx_0;
          cross_val_idx_1 += edgevectemp_idx_2 * b_idx_0 + 0.0 * r;
        }

        i0 = (int)triidx[edgevectemp_idx_0_tmp] - 1;
        V_vec[i0].re += cross_val_idx_2;
        V_vec[i0].im += cross_val_idx_1;
      }
    }
  }
}

/* End of code generation (CalcExciteVecPlaneWave.c) */
