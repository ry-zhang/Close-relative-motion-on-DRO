/*
 * jplEph.h
 *
 * Code generation for function 'jplEph'
 *
 */

#ifndef JPLEPH_H
#define JPLEPH_H

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
extern void jplEph(const emlrtStack *sp, real_T JD, real_T eObj, const real_T
                   CoeffMat[233580], real_T posVel[6]);

#endif

/* End of code generation (jplEph.h) */
