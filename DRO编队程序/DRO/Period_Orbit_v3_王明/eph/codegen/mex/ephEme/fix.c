/*
 * fix.c
 *
 * Code generation for function 'fix'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEme.h"
#include "fix.h"
#include "blas.h"

/* Function Definitions */
void b_fix(real_T *x)
{
  if (*x < 0.0) {
    *x = muDoubleScalarCeil(*x);
  } else {
    *x = muDoubleScalarFloor(*x);
  }
}

/* End of code generation (fix.c) */
