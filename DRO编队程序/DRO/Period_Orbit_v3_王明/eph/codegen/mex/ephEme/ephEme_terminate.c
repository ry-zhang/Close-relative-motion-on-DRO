/*
 * ephEme_terminate.c
 *
 * Code generation for function 'ephEme_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEme.h"
#include "ephEme_terminate.h"
#include "_coder_ephEme_mex.h"
#include "ephEme_data.h"
#include "blas.h"

/* Function Definitions */
void ephEme_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void ephEme_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (ephEme_terminate.c) */
