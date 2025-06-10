/*
 * error.c
 *
 * Code generation for function 'error'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEclip.h"
#include "error.h"
#include "blas.h"

/* Variable Definitions */
static emlrtRTEInfo f_emlrtRTEI = { 17, 9, "error",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\eml\\+coder\\+internal\\error.m"
};

/* Function Definitions */
void error(const emlrtStack *sp)
{
  emlrtErrorWithMessageIdR2012b(sp, &f_emlrtRTEI,
    "Coder:toolbox:reshape_emptyReshapeLimit", 0);
}

/* End of code generation (error.c) */
