/*
 * ephEme_initialize.c
 *
 * Code generation for function 'ephEme_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEme.h"
#include "ephEme_initialize.h"
#include "_coder_ephEme_mex.h"
#include "ephEme_data.h"
#include "blas.h"

/* Function Definitions */
void ephEme_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (ephEme_initialize.c) */
