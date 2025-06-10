/*
 * jplEph.c
 *
 * Code generation for function 'jplEph'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "ephEme.h"
#include "jplEph.h"
#include "cheb_Interp.h"
#include "getCoeff.h"
#include "fix.h"
#include "blas.h"

/* Variable Definitions */
static emlrtRSInfo c_emlrtRSI = { 54, "jplEph",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\jplEph.m"
};

static emlrtRSInfo d_emlrtRSI = { 55, "jplEph",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\jplEph.m"
};

static emlrtRSInfo e_emlrtRSI = { 67, "jplEph",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\jplEph.m"
};

static emlrtRSInfo f_emlrtRSI = { 84, "jplEph",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\jplEph.m"
};

static emlrtBCInfo emlrtBCI = { -1, -1, 40, 6, "PCtemp", "jplEph",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\jplEph.m",
  0 };

static emlrtDCInfo emlrtDCI = { 46, 21, "jplEph",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\jplEph.m",
  1 };

static emlrtBCInfo b_emlrtBCI = { 1, 13, 46, 21, "header", "jplEph",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\jplEph.m",
  0 };

static emlrtRSInfo t_emlrtRSI = { 85, "jplEph",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\jplEph.m"
};

static emlrtRSInfo u_emlrtRSI = { 68, "jplEph",
  "D:\\\\360\xe5\xae\x89\xe5\x85\xa8\xe4\xba\x91\xe7\x9b\x98\xe5\x90\x8c\xe6\xad\xa5\xe7\x89\x88\\[Code]\\\\jplEph\\\\jplEph.m"
};

/* Function Definitions */
void jplEph(const emlrtStack *sp, real_T JD, real_T eObj, const real_T CoeffMat
            [233580], real_T posVel[6])
{
  boolean_T x[229];
  int32_T idx;
  uint8_T ii_data[1];
  int32_T ii_size_idx_0;
  int32_T ii;
  boolean_T exitg1;
  int32_T PCtemp_size[2];
  real_T PCtemp_data[1020];
  real_T dt;
  real_T interval_length;
  static const int16_T header[52] = { 1, 3, 14, 4, 2, 171, 10, 2, 3, 231, 13, 2,
    4, 309, 11, 1, 5, 342, 8, 1, 6, 366, 7, 1, 7, 387, 6, 1, 8, 405, 6, 1, 9,
    423, 6, 1, 10, 441, 13, 8, 11, 753, 11, 2, 12, 819, 10, 4, 13, 899, 10, 4 };

  real_T d0;
  real_T coeff_data[2697];
  int32_T coeff_size[2];
  real_T coeff[39];
  real_T rvMnGeo[6];
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;

  /*  compute ssb pos and vel */
  /*  ----------input----------- */
  /*  JD -- Julian date */
  /*  eObj -- 1~9=Mercury to Pluto */
  /*  eObj -- 10=moon, 11 = sun */
  /*  */
  /*  ----------output -------- */
  /*  % from ssb to celetrial body */
  /*  unit - km, day */
  /*  3*2 mat */
  /* -------------- */
  /*  Copyright(C) 2019 by Hao Zhang @CSU,CAS */
  /*  hao.zhang.zhr@gmail.com */
  /*  */
  /* % code starts here */
  /*  main */
  /*   DE430 structure */
  /*         1 = mercury           8 = neptune */
  /*         2 = venus              9 = pluto */
  /*         3 = EM                10 = moon(geo) */
  /*         4 = mars             11 = sun */
  /*         5 = jupiter           12 = nutation */
  /*         6 = saturn           13 = moon lib */
  /*         7 = uranus */
  /*  legend */
  /*  start */
  /*  # of coeff */
  /*  # of subinterval */
  /*  DE430  % m_Earth/m_Moon */
  /* % get the coeffs by checking  JD */
  for (idx = 0; idx < 229; idx++) {
    x[idx] = ((CoeffMat[idx] <= JD) && (JD < CoeffMat[229 + idx]));
  }

  idx = 0;
  ii_size_idx_0 = 1;
  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii < 230)) {
    if (x[ii - 1]) {
      idx = 1;
      ii_data[0] = (uint8_T)ii;
      exitg1 = true;
    } else {
      ii++;
    }
  }

  if (idx == 0) {
    ii_size_idx_0 = 0;
  }

  PCtemp_size[0] = ii_size_idx_0;
  PCtemp_size[1] = 1020;
  for (idx = 0; idx < 1020; idx++) {
    for (ii = 0; ii < ii_size_idx_0; ii++) {
      PCtemp_data[ii + ii_size_idx_0 * idx] = CoeffMat[(ii_data[ii] + 229 * idx)
        - 1];
    }
  }

  /*  the entire row */
  idx = ii_size_idx_0 * 1020;
  if (!(1 <= idx)) {
    emlrtDynamicBoundsCheckR2012b(1, 1, idx, &emlrtBCI, sp);
  }

  /*  JD at start of interval */
  dt = JD - CoeffMat[ii_data[0 % ii_size_idx_0] - 1];

  /* % get the coeffs of the planet */
  /*  [nc*3, nc*3, nc*3,...] */
  if (eObj != (int32_T)muDoubleScalarFloor(eObj)) {
    emlrtIntegerCheckR2012b(eObj, &emlrtDCI, sp);
  }

  idx = (int32_T)eObj;
  if (!((idx >= 1) && (idx <= 13))) {
    emlrtDynamicBoundsCheckR2012b(idx, 1, 13, &b_emlrtBCI, sp);
  }

  interval_length = 32.0 / (real_T)header[3 + (((int32_T)eObj - 1) << 2)];
  d0 = dt / interval_length;
  b_fix(&d0);
  st.site = &c_emlrtRSI;
  getCoeff(&st, header[1 + (((int32_T)eObj - 1) << 2)], header[2 + (((int32_T)
             eObj - 1) << 2)], header[3 + (((int32_T)eObj - 1) << 2)], d0 + 1.0,
           PCtemp_data, PCtemp_size, coeff_data, coeff_size);
  st.site = &d_emlrtRSI;
  cheb_Interp(&st, coeff_data, coeff_size, 2.0 * (dt - ((d0 + 1.0) - 1.0) *
    interval_length) / interval_length - 1.0, header[2 + (((int32_T)eObj - 1) <<
    2)], posVel);
  for (idx = 0; idx < 3; idx++) {
    posVel[3 + idx] = posVel[3 + idx] * 2.0 / interval_length;
  }

  if (eObj == 3.0) {
    /*  Earth */
    d0 = dt / 4.0;
    b_fix(&d0);
    st.site = &e_emlrtRSI;
    b_getCoeff(&st, d0 + 1.0, PCtemp_data, PCtemp_size, coeff);
    st.site = &u_emlrtRSI;
    b_cheb_Interp(&st, coeff, 2.0 * (dt - ((d0 + 1.0) - 1.0) * 4.0) / 4.0 - 1.0,
                  rvMnGeo);

    /*  Moon @geo */
    for (idx = 0; idx < 3; idx++) {
      rvMnGeo[3 + idx] = rvMnGeo[3 + idx] * 2.0 / 4.0;
    }

    /*  Earth from Solar Bary */
    for (idx = 0; idx < 6; idx++) {
      posVel[idx] -= 0.012150584269940352 * rvMnGeo[idx];
    }
  }

  if (eObj == 10.0) {
    /*  Moon */
    d0 = dt / 16.0;
    b_fix(&d0);
    st.site = &f_emlrtRSI;
    c_getCoeff(&st, d0 + 1.0, PCtemp_data, PCtemp_size, coeff);
    st.site = &t_emlrtRSI;
    b_cheb_Interp(&st, coeff, 2.0 * (dt - ((d0 + 1.0) - 1.0) * 16.0) / 16.0 -
                  1.0, rvMnGeo);

    /*  EM @ssb */
    for (idx = 0; idx < 3; idx++) {
      rvMnGeo[3 + idx] = rvMnGeo[3 + idx] * 2.0 / 16.0;
    }

    for (idx = 0; idx < 6; idx++) {
      posVel[idx] = 0.98784941573005969 * posVel[idx] + rvMnGeo[idx];
    }
  }
}

/* End of code generation (jplEph.c) */
