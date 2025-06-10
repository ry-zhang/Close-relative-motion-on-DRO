/*
 * _coder_ephEclip_api.c
 *
 * Code generation for function '_coder_ephEclip_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEclip.h"
#include "_coder_ephEclip_api.h"
#include "ephEclip_data.h"
#include "blas.h"

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *C_Mat,
  const char_T *identifier))[233580];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[233580];
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *jd, const
  char_T *identifier);
static const mxArray *emlrt_marshallOut(const real_T u[6]);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[233580];

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *C_Mat,
  const char_T *identifier))[233580]
{
  real_T (*y)[233580];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(C_Mat), &thisId);
  emlrtDestroyArray(&C_Mat);
  return y;
}
  static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[233580]
{
  real_T (*y)[233580];
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *jd, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(jd), &thisId);
  emlrtDestroyArray(&jd);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u[6])
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv3[1] = { 0 };

  static const int32_T iv4[1] = { 6 };

  y = NULL;
  m0 = emlrtCreateNumericArray(1, iv3, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)u);
  emlrtSetDimensions((mxArray *)m0, iv4, 1);
  emlrtAssign(&y, m0);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[233580]
{
  real_T (*ret)[233580];
  static const int32_T dims[2] = { 229, 1020 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[233580])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  void ephEclip_api(const mxArray * const prhs[4], const mxArray *plhs[1])
{
  real_T (*rv)[6];
  real_T jd;
  real_T ntarg;
  real_T ncent;
  real_T (*C_Mat)[233580];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  rv = (real_T (*)[6])mxMalloc(sizeof(real_T [6]));

  /* Marshall function inputs */
  jd = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "jd");
  ntarg = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "ntarg");
  ncent = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "ncent");
  C_Mat = c_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "C_Mat");

  /* Invoke the target function */
  ephEclip(&st, jd, ntarg, ncent, *C_Mat, *rv);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*rv);
}

/* End of code generation (_coder_ephEclip_api.c) */
