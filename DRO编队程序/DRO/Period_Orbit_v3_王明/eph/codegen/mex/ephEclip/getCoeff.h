/*
 * getCoeff.h
 *
 * Code generation for function 'getCoeff'
 *
 */

#ifndef GETCOEFF_H
#define GETCOEFF_H

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
#include "ephEclip_types.h"

/* Function Declarations */
extern void b_getCoeff(const emlrtStack *sp, real_T jIntv, const real_T
  PCtemp_data[], const int32_T PCtemp_size[2], real_T coeff[39]);
extern void c_getCoeff(const emlrtStack *sp, real_T jIntv, const real_T
  PCtemp_data[], const int32_T PCtemp_size[2], real_T coeff[39]);
extern void getCoeff(const emlrtStack *sp, real_T icStart, real_T nc, real_T ns,
                     real_T jIntv, const real_T PCtemp_data[], const int32_T
                     PCtemp_size[2], real_T coeff_data[], int32_T coeff_size[2]);

#endif

/* End of code generation (getCoeff.h) */
