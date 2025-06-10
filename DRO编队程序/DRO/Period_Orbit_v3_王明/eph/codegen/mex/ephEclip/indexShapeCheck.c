/*
 * indexShapeCheck.c
 *
 * Code generation for function 'indexShapeCheck'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEclip.h"
#include "indexShapeCheck.h"
#include "ephEclip_data.h"
#include "blas.h"

/* Function Definitions */
void b_indexShapeCheck(const emlrtStack *sp, const int32_T matrixSize[2])
{
  boolean_T nonSingletonDimFound;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  nonSingletonDimFound = false;
  if (matrixSize[0] != 1) {
    nonSingletonDimFound = true;
  }

  if (matrixSize[1] != 1) {
    if (nonSingletonDimFound) {
      nonSingletonDimFound = false;
    } else {
      nonSingletonDimFound = true;
    }
  }

  if (nonSingletonDimFound) {
    st.site = &j_emlrtRSI;
    b_st.site = &k_emlrtRSI;
    if (!!(matrixSize[0] == 1)) {
    } else {
      emlrtErrorWithMessageIdR2012b(&b_st, &e_emlrtRTEI,
        "Coder:FE:PotentialMatrixMatrix", 0);
    }
  }
}

void indexShapeCheck(const emlrtStack *sp, const int32_T matrixSize[2], const
                     int32_T indexSize[2])
{
  boolean_T nonSingletonDimFound;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  nonSingletonDimFound = false;
  if (matrixSize[0] != 1) {
    nonSingletonDimFound = true;
  }

  if (matrixSize[1] != 1) {
    if (nonSingletonDimFound) {
      nonSingletonDimFound = false;
    } else {
      nonSingletonDimFound = true;
    }
  }

  if (nonSingletonDimFound) {
    nonSingletonDimFound = false;
    if (indexSize[0] != 1) {
      nonSingletonDimFound = true;
    }

    if (indexSize[1] != 1) {
      if (nonSingletonDimFound) {
        nonSingletonDimFound = false;
      } else {
        nonSingletonDimFound = true;
      }
    }

    if (nonSingletonDimFound) {
      st.site = &j_emlrtRSI;
      nonSingletonDimFound = !(matrixSize[0] == 1);
      if (nonSingletonDimFound || (indexSize[1] == 1)) {
        nonSingletonDimFound = true;
      }

      b_st.site = &k_emlrtRSI;
      if (!nonSingletonDimFound) {
      } else {
        emlrtErrorWithMessageIdR2012b(&b_st, &e_emlrtRTEI,
          "Coder:FE:PotentialMatrixMatrix", 0);
      }
    }
  }
}

/* End of code generation (indexShapeCheck.c) */
