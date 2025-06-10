/*
 * getCoeff.c
 *
 * Code generation for function 'getCoeff'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEclip.h"
#include "getCoeff.h"
#include "indexShapeCheck.h"
#include "ephEclip_emxutil.h"
#include "error.h"
#include "ephEclip_data.h"
#include "blas.h"

/* Variable Definitions */
static emlrtRSInfo g_emlrtRSI = { 2, "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m"
};

static emlrtRSInfo h_emlrtRSI = { 3, "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m"
};

static emlrtRSInfo i_emlrtRSI = { 6, "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m"
};

static emlrtRSInfo l_emlrtRSI = { 39, "reshape",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\elmat\\reshape.m"
};

static emlrtRSInfo m_emlrtRSI = { 61, "reshape",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\elmat\\reshape.m"
};

static emlrtRSInfo n_emlrtRSI = { 131, "reshape",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\elmat\\reshape.m"
};

static emlrtRTEInfo emlrtRTEI = { 1, 18, "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m"
};

static emlrtRTEInfo b_emlrtRTEI = { 1, 50, "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m"
};

static emlrtBCInfo c_emlrtBCI = { -1, -1, 2, 10, "PCtemp", "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m",
  0 };

static emlrtDCInfo b_emlrtDCI = { 5, 20, "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m",
  1 };

static emlrtBCInfo d_emlrtBCI = { -1, -1, 5, 20, "PCtemp", "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m",
  0 };

static emlrtRTEInfo c_emlrtRTEI = { 71, 15, "reshape",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\elmat\\reshape.m"
};

static emlrtRTEInfo d_emlrtRTEI = { 53, 23, "assertValidSizeArg",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\eml\\+coder\\+internal\\assertValidSizeArg.m"
};

static emlrtBCInfo q_emlrtBCI = { 1, 8, 5, 20, "PCtemp", "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m",
  0 };

static emlrtBCInfo r_emlrtBCI = { 1, 2, 5, 20, "PCtemp", "getCoeff",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\getCoeff.m",
  0 };

/* Function Definitions */
void b_getCoeff(const emlrtStack *sp, real_T jIntv, const real_T PCtemp_data[],
                const int32_T PCtemp_size[2], real_T coeff[39])
{
  int32_T k;
  int32_T i8;
  int32_T i9;
  real_T PCtemp[312];
  real_T b_PCtemp[312];
  real_T b_coeff[39];
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &g_emlrtRSI;
  b_indexShapeCheck(&st, PCtemp_size);
  k = PCtemp_size[0] * 1020;
  for (i8 = 0; i8 < 312; i8++) {
    i9 = 441 + i8;
    if (!(i9 <= k)) {
      emlrtDynamicBoundsCheckR2012b(i9, 1, k, &c_emlrtBCI, sp);
    }
  }

  /*  all subintervals */
  memcpy(&PCtemp[0], &PCtemp_data[440], 312U * sizeof(real_T));
  for (i8 = 0; i8 < 39; i8++) {
    for (i9 = 0; i9 < 8; i9++) {
      b_PCtemp[i9 + (i8 << 3)] = PCtemp[i8 + 39 * i9];
    }
  }

  if (jIntv != (int32_T)muDoubleScalarFloor(jIntv)) {
    emlrtIntegerCheckR2012b(jIntv, &b_emlrtDCI, sp);
  }

  i8 = (int32_T)jIntv;
  if (!((i8 >= 1) && (i8 <= 8))) {
    emlrtDynamicBoundsCheckR2012b(i8, 1, 8, &q_emlrtBCI, sp);
  }

  /*  get the correct subinterval */
  for (k = 0; k < 39; k++) {
    b_coeff[k] = b_PCtemp[((int32_T)jIntv + (k << 3)) - 1];
  }

  for (i8 = 0; i8 < 13; i8++) {
    for (i9 = 0; i9 < 3; i9++) {
      coeff[i9 + 3 * i8] = b_coeff[i8 + 13 * i9];
    }
  }
}

void c_getCoeff(const emlrtStack *sp, real_T jIntv, const real_T PCtemp_data[],
                const int32_T PCtemp_size[2], real_T coeff[39])
{
  boolean_T nonSingletonDimFound;
  int32_T k;
  int32_T i12;
  int32_T i13;
  real_T PCtemp[78];
  real_T b_PCtemp[78];
  real_T b_coeff[39];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &g_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  nonSingletonDimFound = false;
  if (PCtemp_size[0] != 1) {
    nonSingletonDimFound = true;
  }

  if (nonSingletonDimFound) {
    nonSingletonDimFound = false;
  } else {
    nonSingletonDimFound = true;
  }

  if (nonSingletonDimFound) {
    b_st.site = &j_emlrtRSI;
    c_st.site = &k_emlrtRSI;
    if (!!(PCtemp_size[0] == 1)) {
    } else {
      emlrtErrorWithMessageIdR2012b(&c_st, &e_emlrtRTEI,
        "Coder:FE:PotentialMatrixMatrix", 0);
    }
  }

  k = PCtemp_size[0] * 1020;
  for (i12 = 0; i12 < 78; i12++) {
    i13 = 231 + i12;
    if (!(i13 <= k)) {
      emlrtDynamicBoundsCheckR2012b(i13, 1, k, &c_emlrtBCI, sp);
    }
  }

  /*  all subintervals */
  memcpy(&PCtemp[0], &PCtemp_data[230], 78U * sizeof(real_T));
  for (i12 = 0; i12 < 39; i12++) {
    for (i13 = 0; i13 < 2; i13++) {
      b_PCtemp[i13 + (i12 << 1)] = PCtemp[i12 + 39 * i13];
    }
  }

  if (jIntv != (int32_T)muDoubleScalarFloor(jIntv)) {
    emlrtIntegerCheckR2012b(jIntv, &b_emlrtDCI, sp);
  }

  i12 = (int32_T)jIntv;
  if (!((i12 >= 1) && (i12 <= 2))) {
    emlrtDynamicBoundsCheckR2012b(i12, 1, 2, &r_emlrtBCI, sp);
  }

  /*  get the correct subinterval */
  for (k = 0; k < 39; k++) {
    b_coeff[k] = b_PCtemp[((int32_T)jIntv + (k << 1)) - 1];
  }

  for (i12 = 0; i12 < 13; i12++) {
    for (i13 = 0; i13 < 3; i13++) {
      coeff[i13 + 3 * i12] = b_coeff[i12 + 13 * i13];
    }
  }
}

void getCoeff(const emlrtStack *sp, real_T icStart, real_T nc, real_T ns, real_T
              jIntv, const real_T PCtemp_data[], const int32_T PCtemp_size[2],
              real_T coeff_data[], int32_T coeff_size[2])
{
  real_T varargin_1;
  int32_T i2;
  int32_T i3;
  int32_T iv0[2];
  int32_T maxdimlen;
  int32_T x_size_idx_0;
  int32_T loop_ub;
  real_T x_data[1020];
  boolean_T b0;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  emxArray_real_T *PCtemp;
  real_T coeffTemp_data[1020];
  int32_T b_PCtemp_size[2];
  int16_T b_varargin_1[2];
  real_T b_coeff_data[2697];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  varargin_1 = (icStart - 1.0) + 3.0 * nc * ns;
  if (icStart > varargin_1) {
    i3 = 1;
    i2 = 1;
  } else {
    i2 = PCtemp_size[0] * 1020;
    i3 = (int32_T)icStart;
    if (!((i3 >= 1) && (i3 <= i2))) {
      emlrtDynamicBoundsCheckR2012b(i3, 1, i2, &c_emlrtBCI, sp);
    }

    maxdimlen = (int32_T)varargin_1;
    if (!((maxdimlen >= 1) && (maxdimlen <= i2))) {
      emlrtDynamicBoundsCheckR2012b(maxdimlen, 1, i2, &c_emlrtBCI, sp);
    }

    i2 = maxdimlen + 1;
  }

  iv0[0] = 1;
  iv0[1] = i2 - i3;
  st.site = &g_emlrtRSI;
  indexShapeCheck(&st, PCtemp_size, iv0);

  /*  all subintervals */
  st.site = &h_emlrtRSI;
  x_size_idx_0 = i2 - i3;
  loop_ub = i2 - i3;
  for (i2 = 0; i2 < loop_ub; i2++) {
    x_data[i2] = PCtemp_data[(i3 + i2) - 1];
  }

  varargin_1 = nc * 3.0;
  b_st.site = &l_emlrtRSI;
  c_st.site = &n_emlrtRSI;
  if ((!(varargin_1 != varargin_1)) && (!(-2.147483648E+9 > varargin_1))) {
    b0 = true;
  } else {
    b0 = false;
  }

  if (b0) {
  } else {
    emlrtErrorWithMessageIdR2012b(&c_st, &d_emlrtRTEI,
      "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }

  c_st.site = &n_emlrtRSI;
  if ((!(ns != ns)) && (!(-2.147483648E+9 > ns))) {
    b0 = true;
  } else {
    b0 = false;
  }

  if (b0) {
  } else {
    emlrtErrorWithMessageIdR2012b(&c_st, &d_emlrtRTEI,
      "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }

  maxdimlen = (int16_T)x_size_idx_0;
  if (1 > (int16_T)x_size_idx_0) {
    maxdimlen = 1;
  }

  maxdimlen = muIntScalarMax_sint32(x_size_idx_0, maxdimlen);
  if ((int32_T)varargin_1 > maxdimlen) {
    b_st.site = &m_emlrtRSI;
    error(&b_st);
  }

  if ((int32_T)ns > maxdimlen) {
    b_st.site = &m_emlrtRSI;
    error(&b_st);
  }

  emxInit_real_T(&st, &y, 2, &emlrtRTEI, true);
  i2 = y->size[0] * y->size[1];
  y->size[0] = (int32_T)varargin_1;
  y->size[1] = (int32_T)ns;
  emxEnsureCapacity(&st, (emxArray__common *)y, i2, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  if (x_size_idx_0 == y->size[0] * y->size[1]) {
  } else {
    emlrtErrorWithMessageIdR2012b(&st, &c_emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  for (maxdimlen = 0; maxdimlen + 1 <= x_size_idx_0; maxdimlen++) {
    y->data[maxdimlen] = x_data[maxdimlen];
  }

  emxInit_real_T(&st, &b_y, 2, &emlrtRTEI, true);
  i2 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = y->size[1];
  b_y->size[1] = y->size[0];
  emxEnsureCapacity(sp, (emxArray__common *)b_y, i2, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = y->size[0];
  for (i2 = 0; i2 < loop_ub; i2++) {
    maxdimlen = y->size[1];
    for (i3 = 0; i3 < maxdimlen; i3++) {
      b_y->data[i3 + b_y->size[0] * i2] = y->data[i2 + y->size[0] * i3];
    }
  }

  emxInit_real_T(sp, &PCtemp, 2, &b_emlrtRTEI, true);
  maxdimlen = y->size[1];
  loop_ub = y->size[0];
  i2 = PCtemp->size[0] * PCtemp->size[1];
  PCtemp->size[0] = maxdimlen;
  PCtemp->size[1] = loop_ub;
  emxEnsureCapacity(sp, (emxArray__common *)PCtemp, i2, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  emxFree_real_T(&y);
  for (i2 = 0; i2 < loop_ub; i2++) {
    for (i3 = 0; i3 < maxdimlen; i3++) {
      PCtemp->data[i3 + PCtemp->size[0] * i2] = b_y->data[i3 + maxdimlen * i2];
    }
  }

  emxFree_real_T(&b_y);
  loop_ub = PCtemp->size[1];
  if (jIntv != (int32_T)muDoubleScalarFloor(jIntv)) {
    emlrtIntegerCheckR2012b(jIntv, &b_emlrtDCI, sp);
  }

  i2 = PCtemp->size[0];
  maxdimlen = (int32_T)jIntv;
  if (!((maxdimlen >= 1) && (maxdimlen <= i2))) {
    emlrtDynamicBoundsCheckR2012b(maxdimlen, 1, i2, &d_emlrtBCI, sp);
  }

  for (i2 = 0; i2 < loop_ub; i2++) {
    coeffTemp_data[i2] = PCtemp->data[(maxdimlen + PCtemp->size[0] * i2) - 1];
  }

  /*  get the correct subinterval */
  st.site = &i_emlrtRSI;
  i2 = PCtemp->size[1];
  b_st.site = &l_emlrtRSI;
  c_st.site = &n_emlrtRSI;
  if ((!(nc != nc)) && (!(-2.147483648E+9 > nc))) {
    b0 = true;
  } else {
    b0 = false;
  }

  if (b0) {
  } else {
    emlrtErrorWithMessageIdR2012b(&c_st, &d_emlrtRTEI,
      "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }

  i3 = PCtemp->size[1];
  b_PCtemp_size[0] = 1;
  b_PCtemp_size[1] = i3;
  for (i3 = 0; i3 < 2; i3++) {
    b_varargin_1[i3] = (int16_T)b_PCtemp_size[i3];
  }

  maxdimlen = 1;
  if (b_varargin_1[1] > 1) {
    maxdimlen = b_varargin_1[1];
  }

  if (i2 > maxdimlen) {
    maxdimlen = PCtemp->size[1];
  }

  emxFree_real_T(&PCtemp);
  if ((int32_T)nc > maxdimlen) {
    b_st.site = &m_emlrtRSI;
    error(&b_st);
  }

  if (3 > maxdimlen) {
    b_st.site = &m_emlrtRSI;
    error(&b_st);
  }

  if (i2 == (int32_T)nc * 3) {
  } else {
    emlrtErrorWithMessageIdR2012b(&st, &c_emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  for (maxdimlen = 0; maxdimlen + 1 <= i2; maxdimlen++) {
    b_coeff_data[maxdimlen] = coeffTemp_data[maxdimlen];
  }

  coeff_size[0] = 3;
  coeff_size[1] = (int32_T)nc;
  loop_ub = (int32_T)nc;
  for (i2 = 0; i2 < loop_ub; i2++) {
    for (i3 = 0; i3 < 3; i3++) {
      coeff_data[i3 + 3 * i2] = b_coeff_data[i2 + (int32_T)nc * i3];
    }
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (getCoeff.c) */
