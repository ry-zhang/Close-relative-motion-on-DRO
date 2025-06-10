/*
 * ephEclip_initialize.c
 *
 * Code generation for function 'ephEclip_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEclip.h"
#include "ephEclip_initialize.h"
#include "_coder_ephEclip_mex.h"
#include "ephEclip_data.h"
#include "blas.h"

/* Function Definitions */
void ephEclip_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (ephEclip_initialize.c) */
