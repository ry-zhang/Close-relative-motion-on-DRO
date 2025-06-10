/*
 * ephEme.h
 *
 * Code generation for function 'ephEme'
 *
 */

#ifndef EPHEME_H
#define EPHEME_H

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
extern void ephEme(const emlrtStack *sp, real_T jd, real_T ntarg, real_T ncent,
                   const real_T C_Mat[233580], real_T rv[6]);

#endif

/* End of code generation (ephEme.h) */
