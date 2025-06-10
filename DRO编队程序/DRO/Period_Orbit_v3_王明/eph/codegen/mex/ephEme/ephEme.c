/*
 * ephEme.c
 *
 * Code generation for function 'ephEme'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEme.h"
#include "jplEph.h"
#include "blas.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 8, "ephEme",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\ephEme.m"
};

static emlrtRSInfo b_emlrtRSI = { 9, "ephEme",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\ephEme.m"
};

/* Function Definitions */
void ephEme(const emlrtStack *sp, real_T jd, real_T ntarg, real_T ncent, const
            real_T C_Mat[233580], real_T rv[6])
{
  real_T rvCent[6];
  real_T rvTarg[6];
  int32_T i0;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;

  /*  */
  /*  ¼ÆËãeme2000ÐÇÀú */
  /*   */
  /* %%%%%%%%%%%%%%%%%%%%%% */
  st.site = &emlrtRSI;
  jplEph(&st, jd, ncent, C_Mat, rvCent);
  st.site = &b_emlrtRSI;
  jplEph(&st, jd, ntarg, C_Mat, rvTarg);
  for (i0 = 0; i0 < 3; i0++) {
    rv[i0] = rvTarg[i0] - rvCent[i0];
    rv[i0 + 3] = (rvTarg[3 + i0] - rvCent[3 + i0]) / 86400.0;
  }
}

/* End of code generation (ephEme.c) */
