/*
 * ephEclip.h
 *
 * Code generation for function 'ephEclip'
 *
 */

#ifndef EPHECLIP_H
#define EPHECLIP_H

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
extern void ephEclip(const emlrtStack *sp, real_T jd, real_T ntarg, real_T ncent,
                     const real_T C_Mat[233580], real_T rv[6]);

#endif

/* End of code generation (ephEclip.h) */
