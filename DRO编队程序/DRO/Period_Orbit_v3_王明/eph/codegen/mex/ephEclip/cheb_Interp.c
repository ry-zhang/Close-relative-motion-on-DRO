/*
 * cheb_Interp.c
 *
 * Code generation for function 'cheb_Interp'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEclip.h"
#include "cheb_Interp.h"
#include "ephEclip_data.h"
#include "blas.h"

/* Variable Definitions */
static emlrtRSInfo o_emlrtRSI = { 18, "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m"
};

static emlrtRSInfo p_emlrtRSI = { 19, "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m"
};

static emlrtRSInfo q_emlrtRSI = { 68, "eml_mtimes_helper",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_helper.m"
};

static emlrtRSInfo r_emlrtRSI = { 21, "eml_mtimes_helper",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_helper.m"
};

static emlrtRTEInfo g_emlrtRTEI = { 98, 23, "eml_mtimes_helper",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_helper.m"
};

static emlrtRTEInfo h_emlrtRTEI = { 103, 23, "eml_mtimes_helper",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_helper.m"
};

static emlrtECInfo emlrtECI = { -1, 16, 1, "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m"
};

static emlrtBCInfo e_emlrtBCI = { -1, -1, 16, 1, "pc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtECInfo b_emlrtECI = { 2, 16, 14, "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m"
};

static emlrtBCInfo f_emlrtBCI = { -1, -1, 16, 42, "pc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 16, 24, "pc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo h_emlrtBCI = { -1, -1, 9, 1, "vc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 7, 1, "pc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo j_emlrtBCI = { -1, -1, 12, 20, "pc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo k_emlrtBCI = { -1, -1, 12, 29, "pc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo l_emlrtBCI = { -1, -1, 12, 5, "pc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo m_emlrtBCI = { -1, -1, 13, 16, "pc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo n_emlrtBCI = { -1, -1, 13, 31, "vc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo o_emlrtBCI = { -1, -1, 13, 40, "vc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

static emlrtBCInfo p_emlrtBCI = { -1, -1, 13, 5, "vc", "cheb_Interp",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\cheb_Interp.m",
  0 };

/* Function Definitions */
void b_cheb_Interp(const emlrtStack *sp, const real_T coeff[39], real_T tau,
                   real_T rv[6])
{
  int32_T i10;
  real_T pc[13];
  real_T vc[13];
  int32_T ii;
  real_T a;
  real_T b_a[11];
  int32_T i11;

  /*  km, 32/ns-day */
  for (i10 = 0; i10 < 6; i10++) {
    rv[i10] = 0.0;
  }

  memset(&pc[0], 0, 13U * sizeof(real_T));
  memset(&vc[0], 0, 13U * sizeof(real_T));
  pc[0] = 1.0;
  pc[1] = tau;
  vc[0] = 0.0;
  vc[1] = 1.0;
  ii = 0;
  while (ii < 11) {
    pc[ii + 2] = 2.0 * tau * pc[ii + 1] - pc[ii];
    vc[ii + 2] = (2.0 * pc[ii + 1] + 2.0 * tau * vc[ii + 1]) - vc[ii];
    ii++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  a = 2.0 * tau;
  for (i10 = 0; i10 < 11; i10++) {
    b_a[i10] = a * pc[1 + i10] - pc[i10];
  }

  memcpy(&pc[2], &b_a[0], 11U * sizeof(real_T));
  for (i10 = 0; i10 < 3; i10++) {
    rv[i10] = 0.0;
    rv[i10] = 0.0;
    rv[3 + i10] = 0.0;
    rv[3 + i10] = 0.0;
    for (i11 = 0; i11 < 13; i11++) {
      rv[i10] += coeff[i10 + 3 * i11] * pc[i11];
      rv[3 + i10] += coeff[i10 + 3 * i11] * vc[i11];
    }
  }
}

void cheb_Interp(const emlrtStack *sp, const real_T coeff_data[], const int32_T
                 coeff_size[2], real_T tau, real_T nc, real_T rv[6])
{
  int32_T i4;
  int32_T loop_ub;
  real_T pc_data[899];
  real_T vc_data[899];
  int32_T ii;
  int32_T i5;
  int32_T i6;
  real_T a;
  int32_T tmp_size[2];
  int32_T iv1[2];
  real_T tmp_data[898];
  int32_T iv2[2];
  int32_T i7;
  real_T b_tmp_data[899];
  real_T b_data[899];
  real_T beta1;
  char_T TRANSB;
  char_T TRANSA;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;

  /*  km, 32/ns-day */
  for (i4 = 0; i4 < 6; i4++) {
    rv[i4] = 0.0;
  }

  loop_ub = (int32_T)nc;
  for (i4 = 0; i4 < loop_ub; i4++) {
    pc_data[i4] = 0.0;
  }

  loop_ub = (int32_T)nc;
  for (i4 = 0; i4 < loop_ub; i4++) {
    vc_data[i4] = 0.0;
  }

  pc_data[0] = 1.0;
  i4 = (int32_T)nc;
  if (!(2 <= i4)) {
    emlrtDynamicBoundsCheckR2012b(2, 1, i4, &i_emlrtBCI, sp);
  }

  pc_data[1] = tau;
  vc_data[0] = 0.0;
  i4 = (int32_T)nc;
  if (!(2 <= i4)) {
    emlrtDynamicBoundsCheckR2012b(2, 1, i4, &h_emlrtBCI, sp);
  }

  vc_data[1] = 1.0;
  ii = 0;
  while (ii <= (int32_T)(nc + -2.0) - 1) {
    i4 = (int32_T)nc;
    i5 = (int32_T)((3.0 + (real_T)ii) - 1.0);
    if (!((i5 >= 1) && (i5 <= i4))) {
      emlrtDynamicBoundsCheckR2012b(i5, 1, i4, &j_emlrtBCI, sp);
    }

    i4 = (int32_T)nc;
    loop_ub = (int32_T)((3.0 + (real_T)ii) - 2.0);
    if (!((loop_ub >= 1) && (loop_ub <= i4))) {
      emlrtDynamicBoundsCheckR2012b(loop_ub, 1, i4, &k_emlrtBCI, sp);
    }

    i4 = (int32_T)nc;
    i6 = (int32_T)(3.0 + (real_T)ii);
    if (!((i6 >= 1) && (i6 <= i4))) {
      emlrtDynamicBoundsCheckR2012b(i6, 1, i4, &l_emlrtBCI, sp);
    }

    pc_data[i6 - 1] = 2.0 * tau * pc_data[i5 - 1] - pc_data[loop_ub - 1];
    i4 = (int32_T)nc;
    i5 = (int32_T)((3.0 + (real_T)ii) - 1.0);
    if (!((i5 >= 1) && (i5 <= i4))) {
      emlrtDynamicBoundsCheckR2012b(i5, 1, i4, &m_emlrtBCI, sp);
    }

    i4 = (int32_T)nc;
    loop_ub = (int32_T)((3.0 + (real_T)ii) - 1.0);
    if (!((loop_ub >= 1) && (loop_ub <= i4))) {
      emlrtDynamicBoundsCheckR2012b(loop_ub, 1, i4, &n_emlrtBCI, sp);
    }

    i4 = (int32_T)nc;
    i6 = (int32_T)((3.0 + (real_T)ii) - 2.0);
    if (!((i6 >= 1) && (i6 <= i4))) {
      emlrtDynamicBoundsCheckR2012b(i6, 1, i4, &o_emlrtBCI, sp);
    }

    i4 = (int32_T)nc;
    i7 = (int32_T)(3.0 + (real_T)ii);
    if (!((i7 >= 1) && (i7 <= i4))) {
      emlrtDynamicBoundsCheckR2012b(i7, 1, i4, &p_emlrtBCI, sp);
    }

    vc_data[i7 - 1] = (2.0 * pc_data[i5 - 1] + 2.0 * tau * vc_data[loop_ub - 1])
      - vc_data[i6 - 1];
    ii++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  if (2.0 > nc - 1.0) {
    i4 = 0;
    loop_ub = 0;
  } else {
    i4 = (int32_T)nc;
    if (!(2 <= i4)) {
      emlrtDynamicBoundsCheckR2012b(2, 1, i4, &g_emlrtBCI, sp);
    }

    i4 = 1;
    i5 = (int32_T)nc;
    loop_ub = (int32_T)(nc - 1.0);
    if (!((loop_ub >= 1) && (loop_ub <= i5))) {
      emlrtDynamicBoundsCheckR2012b(loop_ub, 1, i5, &g_emlrtBCI, sp);
    }
  }

  if (1.0 > nc - 2.0) {
    i6 = 0;
  } else {
    i5 = (int32_T)nc;
    i6 = (int32_T)(nc - 2.0);
    if (!((i6 >= 1) && (i6 <= i5))) {
      emlrtDynamicBoundsCheckR2012b(i6, 1, i5, &f_emlrtBCI, sp);
    }
  }

  a = 2.0 * tau;
  tmp_size[0] = 1;
  tmp_size[1] = loop_ub - i4;
  loop_ub -= i4;
  for (i5 = 0; i5 < loop_ub; i5++) {
    tmp_data[i5] = a * pc_data[i4 + i5];
  }

  for (i4 = 0; i4 < 2; i4++) {
    iv1[i4] = tmp_size[i4];
  }

  iv2[0] = 1;
  iv2[1] = i6;
  if ((iv1[0] != iv2[0]) || (iv1[1] != iv2[1])) {
    emlrtSizeEqCheckNDR2012b(&iv1[0], &iv2[0], &b_emlrtECI, sp);
  }

  if (3.0 > nc) {
    i4 = 0;
    loop_ub = 0;
  } else {
    i4 = (int32_T)nc;
    if (!(3 <= i4)) {
      emlrtDynamicBoundsCheckR2012b(3, 1, i4, &e_emlrtBCI, sp);
    }

    i4 = 2;
    i5 = (int32_T)nc;
    loop_ub = (int32_T)nc;
    if (!((loop_ub >= 1) && (loop_ub <= i5))) {
      emlrtDynamicBoundsCheckR2012b(loop_ub, 1, i5, &e_emlrtBCI, sp);
    }
  }

  i5 = loop_ub - i4;
  if (i5 != tmp_size[1]) {
    emlrtSizeEqCheck1DR2012b(i5, tmp_size[1], &emlrtECI, sp);
  }

  loop_ub = tmp_size[1];
  for (i5 = 0; i5 < loop_ub; i5++) {
    b_tmp_data[i5] = tmp_data[i5] - pc_data[i5];
  }

  loop_ub = tmp_size[1];
  for (i5 = 0; i5 < loop_ub; i5++) {
    pc_data[i4 + i5] = b_tmp_data[i5];
  }

  st.site = &o_emlrtRSI;
  loop_ub = (int32_T)nc;
  for (i4 = 0; i4 < loop_ub; i4++) {
    b_data[i4] = pc_data[i4];
  }

  b_st.site = &r_emlrtRSI;
  if (!(coeff_size[1] == (int32_T)nc)) {
    if ((int32_T)nc == 1) {
      emlrtErrorWithMessageIdR2012b(&b_st, &g_emlrtRTEI,
        "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
    } else {
      emlrtErrorWithMessageIdR2012b(&b_st, &h_emlrtRTEI, "Coder:MATLAB:innerdim",
        0);
    }
  }

  if ((coeff_size[1] == 1) || ((int32_T)nc == 1)) {
    for (i4 = 0; i4 < 3; i4++) {
      rv[i4] = 0.0;
      loop_ub = coeff_size[1];
      for (i5 = 0; i5 < loop_ub; i5++) {
        rv[i4] += coeff_data[i4 + coeff_size[0] * i5] * b_data[i5];
      }
    }
  } else {
    b_st.site = &q_emlrtRSI;
    for (i4 = 0; i4 < 3; i4++) {
      rv[i4] = 0.0;
    }

    a = 1.0;
    beta1 = 0.0;
    TRANSB = 'N';
    TRANSA = 'N';
    m_t = (ptrdiff_t)3;
    n_t = (ptrdiff_t)1;
    k_t = (ptrdiff_t)coeff_size[1];
    lda_t = (ptrdiff_t)3;
    ldb_t = (ptrdiff_t)coeff_size[1];
    ldc_t = (ptrdiff_t)3;
    dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, &a, &coeff_data[0], &lda_t,
          &b_data[0], &ldb_t, &beta1, &rv[0], &ldc_t);
  }

  st.site = &p_emlrtRSI;
  loop_ub = (int32_T)nc;
  for (i4 = 0; i4 < loop_ub; i4++) {
    b_data[i4] = vc_data[i4];
  }

  b_st.site = &r_emlrtRSI;
  if (!(coeff_size[1] == (int32_T)nc)) {
    if ((int32_T)nc == 1) {
      emlrtErrorWithMessageIdR2012b(&b_st, &g_emlrtRTEI,
        "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
    } else {
      emlrtErrorWithMessageIdR2012b(&b_st, &h_emlrtRTEI, "Coder:MATLAB:innerdim",
        0);
    }
  }

  if ((coeff_size[1] == 1) || ((int32_T)nc == 1)) {
    for (i4 = 0; i4 < 3; i4++) {
      rv[3 + i4] = 0.0;
      loop_ub = coeff_size[1];
      for (i5 = 0; i5 < loop_ub; i5++) {
        rv[3 + i4] += coeff_data[i4 + coeff_size[0] * i5] * b_data[i5];
      }
    }
  } else {
    b_st.site = &q_emlrtRSI;
    for (i4 = 0; i4 < 3; i4++) {
      rv[3 + i4] = 0.0;
    }

    a = 1.0;
    beta1 = 0.0;
    TRANSB = 'N';
    TRANSA = 'N';
    m_t = (ptrdiff_t)3;
    n_t = (ptrdiff_t)1;
    k_t = (ptrdiff_t)coeff_size[1];
    lda_t = (ptrdiff_t)3;
    ldb_t = (ptrdiff_t)coeff_size[1];
    ldc_t = (ptrdiff_t)3;
    dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, &a, &coeff_data[0], &lda_t,
          &b_data[0], &ldb_t, &beta1, &rv[3], &ldc_t);
  }
}

/* End of code generation (cheb_Interp.c) */
