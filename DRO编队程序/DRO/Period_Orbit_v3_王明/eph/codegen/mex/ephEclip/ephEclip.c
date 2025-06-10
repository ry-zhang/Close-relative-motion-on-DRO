/*
 * ephEclip.c
 *
 * Code generation for function 'ephEclip'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEclip.h"
#include "jplEph.h"
#include "blas.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 13, "ephEclip",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\ephEclip.m"
};

static emlrtRSInfo b_emlrtRSI = { 14, "ephEclip",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\ephEclip.m"
};

/* Function Definitions */
void ephEclip(const emlrtStack *sp, real_T jd, real_T ntarg, real_T ncent, const
              real_T C_Mat[233580], real_T rv[6])
{
  real_T rvCent[6];
  real_T rvTarg[6];
  real_T b_rvTarg[3];
  real_T c_rvTarg[3];
  int32_T i0;
  real_T a[3];
  real_T b_a[3];
  real_T d0;
  int32_T i1;
  static const real_T c_a[9] = { 1.0, 4.4036E-7, -1.90919E-7, -4.79966E-7,
    0.917482137087, -0.397776982902, 0.0, 0.397776982902, 0.917482137087 };

  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;

  /*  */
  /*  ¼ÆËãeclipÐÇÀú */
  /*  */
  /* %%%%%%%%%%%%%%%%%%%%%% */
  /*  eme2000 -> eclip */
  st.site = &emlrtRSI;
  jplEph(&st, jd, ncent, C_Mat, rvCent);
  st.site = &b_emlrtRSI;
  jplEph(&st, jd, ntarg, C_Mat, rvTarg);
  for (i0 = 0; i0 < 3; i0++) {
    b_rvTarg[i0] = rvTarg[i0] - rvCent[i0];
    c_rvTarg[i0] = rvTarg[3 + i0] - rvCent[3 + i0];
  }

  for (i0 = 0; i0 < 3; i0++) {
    a[i0] = 0.0;
    d0 = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      d0 += c_a[i0 + 3 * i1] * c_rvTarg[i1];
      a[i0] += c_a[i0 + 3 * i1] * b_rvTarg[i1];
    }

    b_a[i0] = d0 / 86400.0;
    rv[i0] = a[i0];
    rv[i0 + 3] = b_a[i0];
  }
}

/* End of code generation (ephEclip.c) */
