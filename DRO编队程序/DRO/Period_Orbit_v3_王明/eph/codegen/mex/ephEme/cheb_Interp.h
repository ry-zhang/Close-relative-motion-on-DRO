/*
 * cheb_Interp.h
 *
 * Code generation for function 'cheb_Interp'
 *
 */

#ifndef CHEB_INTERP_H
#define CHEB_INTERP_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "ephEme_types.h"

/* Function Declarations */
extern void b_cheb_Interp(const emlrtStack *sp, const real_T coeff[39], real_T
  tau, real_T rv[6]);
extern void cheb_Interp(const emlrtStack *sp, const real_T coeff_data[], const
  int32_T coeff_size[2], real_T tau, real_T nc, real_T rv[6]);

#endif

/* End of code generation (cheb_Interp.h) */
