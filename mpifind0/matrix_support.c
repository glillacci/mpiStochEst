/* Produced by CVXGEN, 2012-05-16 08:23:25 -0700.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */

#include "solver.h"

void multbymA(double *lhs, double *rhs) {
}

void multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
  lhs[3] = 0;
  lhs[4] = 0;
  lhs[5] = 0;
  lhs[6] = 0;
  lhs[7] = 0;
  lhs[8] = 0;
  lhs[9] = 0;
  lhs[10] = 0;
  lhs[11] = 0;
  lhs[12] = 0;
  lhs[13] = 0;
}

void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1);
  lhs[1] = -rhs[1]*(-1);
  lhs[2] = -rhs[2]*(-1);
  lhs[3] = -rhs[3]*(-1);
  lhs[4] = -rhs[4]*(-1);
  lhs[5] = -rhs[5]*(-1);
  lhs[6] = -rhs[6]*(-1);
  lhs[7] = -rhs[7]*(-1);
  lhs[8] = -rhs[8]*(-1);
  lhs[9] = -rhs[9]*(-1);
  lhs[10] = -rhs[10]*(-1);
  lhs[11] = -rhs[11]*(-1);
  lhs[12] = -rhs[12]*(-1);
  lhs[13] = -rhs[0]*(1);
  lhs[14] = -rhs[1]*(1);
  lhs[15] = -rhs[2]*(1);
  lhs[16] = -rhs[3]*(1);
  lhs[17] = -rhs[4]*(1);
  lhs[18] = -rhs[5]*(1);
  lhs[19] = -rhs[6]*(1);
  lhs[20] = -rhs[7]*(1);
  lhs[21] = -rhs[8]*(1);
  lhs[22] = -rhs[9]*(1);
  lhs[23] = -rhs[10]*(1);
  lhs[24] = -rhs[11]*(1);
  lhs[25] = -rhs[12]*(1);
  lhs[26] = -rhs[13]*(1);
  lhs[27] = -rhs[0]*(1)-rhs[13]*(-1);
  lhs[28] = -rhs[1]*(1)-rhs[13]*(-1);
  lhs[29] = -rhs[2]*(1)-rhs[13]*(-1);
  lhs[30] = -rhs[3]*(1)-rhs[13]*(-1);
  lhs[31] = -rhs[4]*(1)-rhs[13]*(-1);
  lhs[32] = -rhs[5]*(1)-rhs[13]*(-1);
  lhs[33] = -rhs[6]*(1)-rhs[13]*(-1);
  lhs[34] = -rhs[7]*(1)-rhs[13]*(-1);
  lhs[35] = -rhs[8]*(1)-rhs[13]*(-1);
  lhs[36] = -rhs[9]*(1)-rhs[13]*(-1);
  lhs[37] = -rhs[10]*(1)-rhs[13]*(-1);
  lhs[38] = -rhs[11]*(1)-rhs[13]*(-1);
  lhs[39] = -rhs[12]*(1)-rhs[13]*(-1);
  lhs[40] = -rhs[0]*(-1)-rhs[13]*(-1);
  lhs[41] = -rhs[1]*(-1)-rhs[13]*(-1);
  lhs[42] = -rhs[2]*(-1)-rhs[13]*(-1);
  lhs[43] = -rhs[3]*(-1)-rhs[13]*(-1);
  lhs[44] = -rhs[4]*(-1)-rhs[13]*(-1);
  lhs[45] = -rhs[5]*(-1)-rhs[13]*(-1);
  lhs[46] = -rhs[6]*(-1)-rhs[13]*(-1);
  lhs[47] = -rhs[7]*(-1)-rhs[13]*(-1);
  lhs[48] = -rhs[8]*(-1)-rhs[13]*(-1);
  lhs[49] = -rhs[9]*(-1)-rhs[13]*(-1);
  lhs[50] = -rhs[10]*(-1)-rhs[13]*(-1);
  lhs[51] = -rhs[11]*(-1)-rhs[13]*(-1);
  lhs[52] = -rhs[12]*(-1)-rhs[13]*(-1);
}

void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1)-rhs[13]*(1)-rhs[27]*(1)-rhs[40]*(-1);
  lhs[1] = -rhs[1]*(-1)-rhs[14]*(1)-rhs[28]*(1)-rhs[41]*(-1);
  lhs[2] = -rhs[2]*(-1)-rhs[15]*(1)-rhs[29]*(1)-rhs[42]*(-1);
  lhs[3] = -rhs[3]*(-1)-rhs[16]*(1)-rhs[30]*(1)-rhs[43]*(-1);
  lhs[4] = -rhs[4]*(-1)-rhs[17]*(1)-rhs[31]*(1)-rhs[44]*(-1);
  lhs[5] = -rhs[5]*(-1)-rhs[18]*(1)-rhs[32]*(1)-rhs[45]*(-1);
  lhs[6] = -rhs[6]*(-1)-rhs[19]*(1)-rhs[33]*(1)-rhs[46]*(-1);
  lhs[7] = -rhs[7]*(-1)-rhs[20]*(1)-rhs[34]*(1)-rhs[47]*(-1);
  lhs[8] = -rhs[8]*(-1)-rhs[21]*(1)-rhs[35]*(1)-rhs[48]*(-1);
  lhs[9] = -rhs[9]*(-1)-rhs[22]*(1)-rhs[36]*(1)-rhs[49]*(-1);
  lhs[10] = -rhs[10]*(-1)-rhs[23]*(1)-rhs[37]*(1)-rhs[50]*(-1);
  lhs[11] = -rhs[11]*(-1)-rhs[24]*(1)-rhs[38]*(1)-rhs[51]*(-1);
  lhs[12] = -rhs[12]*(-1)-rhs[25]*(1)-rhs[39]*(1)-rhs[52]*(-1);
  lhs[13] = -rhs[26]*(1)-rhs[27]*(-1)-rhs[28]*(-1)-rhs[29]*(-1)-rhs[30]*(-1)-rhs[31]*(-1)-rhs[32]*(-1)-rhs[33]*(-1)-rhs[34]*(-1)-rhs[35]*(-1)-rhs[36]*(-1)-rhs[37]*(-1)-rhs[38]*(-1)-rhs[39]*(-1)-rhs[40]*(-1)-rhs[41]*(-1)-rhs[42]*(-1)-rhs[43]*(-1)-rhs[44]*(-1)-rhs[45]*(-1)-rhs[46]*(-1)-rhs[47]*(-1)-rhs[48]*(-1)-rhs[49]*(-1)-rhs[50]*(-1)-rhs[51]*(-1)-rhs[52]*(-1);
}

void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(params.H[0])+rhs[1]*(params.H[13])+rhs[2]*(params.H[26])+rhs[3]*(params.H[39])+rhs[4]*(params.H[52])+rhs[5]*(params.H[65])+rhs[6]*(params.H[78])+rhs[7]*(params.H[91])+rhs[8]*(params.H[104])+rhs[9]*(params.H[117])+rhs[10]*(params.H[130])+rhs[11]*(params.H[143])+rhs[12]*(params.H[156]);
  lhs[1] = rhs[0]*(params.H[1])+rhs[1]*(params.H[14])+rhs[2]*(params.H[27])+rhs[3]*(params.H[40])+rhs[4]*(params.H[53])+rhs[5]*(params.H[66])+rhs[6]*(params.H[79])+rhs[7]*(params.H[92])+rhs[8]*(params.H[105])+rhs[9]*(params.H[118])+rhs[10]*(params.H[131])+rhs[11]*(params.H[144])+rhs[12]*(params.H[157]);
  lhs[2] = rhs[0]*(params.H[2])+rhs[1]*(params.H[15])+rhs[2]*(params.H[28])+rhs[3]*(params.H[41])+rhs[4]*(params.H[54])+rhs[5]*(params.H[67])+rhs[6]*(params.H[80])+rhs[7]*(params.H[93])+rhs[8]*(params.H[106])+rhs[9]*(params.H[119])+rhs[10]*(params.H[132])+rhs[11]*(params.H[145])+rhs[12]*(params.H[158]);
  lhs[3] = rhs[0]*(params.H[3])+rhs[1]*(params.H[16])+rhs[2]*(params.H[29])+rhs[3]*(params.H[42])+rhs[4]*(params.H[55])+rhs[5]*(params.H[68])+rhs[6]*(params.H[81])+rhs[7]*(params.H[94])+rhs[8]*(params.H[107])+rhs[9]*(params.H[120])+rhs[10]*(params.H[133])+rhs[11]*(params.H[146])+rhs[12]*(params.H[159]);
  lhs[4] = rhs[0]*(params.H[4])+rhs[1]*(params.H[17])+rhs[2]*(params.H[30])+rhs[3]*(params.H[43])+rhs[4]*(params.H[56])+rhs[5]*(params.H[69])+rhs[6]*(params.H[82])+rhs[7]*(params.H[95])+rhs[8]*(params.H[108])+rhs[9]*(params.H[121])+rhs[10]*(params.H[134])+rhs[11]*(params.H[147])+rhs[12]*(params.H[160]);
  lhs[5] = rhs[0]*(params.H[5])+rhs[1]*(params.H[18])+rhs[2]*(params.H[31])+rhs[3]*(params.H[44])+rhs[4]*(params.H[57])+rhs[5]*(params.H[70])+rhs[6]*(params.H[83])+rhs[7]*(params.H[96])+rhs[8]*(params.H[109])+rhs[9]*(params.H[122])+rhs[10]*(params.H[135])+rhs[11]*(params.H[148])+rhs[12]*(params.H[161]);
  lhs[6] = rhs[0]*(params.H[6])+rhs[1]*(params.H[19])+rhs[2]*(params.H[32])+rhs[3]*(params.H[45])+rhs[4]*(params.H[58])+rhs[5]*(params.H[71])+rhs[6]*(params.H[84])+rhs[7]*(params.H[97])+rhs[8]*(params.H[110])+rhs[9]*(params.H[123])+rhs[10]*(params.H[136])+rhs[11]*(params.H[149])+rhs[12]*(params.H[162]);
  lhs[7] = rhs[0]*(params.H[7])+rhs[1]*(params.H[20])+rhs[2]*(params.H[33])+rhs[3]*(params.H[46])+rhs[4]*(params.H[59])+rhs[5]*(params.H[72])+rhs[6]*(params.H[85])+rhs[7]*(params.H[98])+rhs[8]*(params.H[111])+rhs[9]*(params.H[124])+rhs[10]*(params.H[137])+rhs[11]*(params.H[150])+rhs[12]*(params.H[163]);
  lhs[8] = rhs[0]*(params.H[8])+rhs[1]*(params.H[21])+rhs[2]*(params.H[34])+rhs[3]*(params.H[47])+rhs[4]*(params.H[60])+rhs[5]*(params.H[73])+rhs[6]*(params.H[86])+rhs[7]*(params.H[99])+rhs[8]*(params.H[112])+rhs[9]*(params.H[125])+rhs[10]*(params.H[138])+rhs[11]*(params.H[151])+rhs[12]*(params.H[164]);
  lhs[9] = rhs[0]*(params.H[9])+rhs[1]*(params.H[22])+rhs[2]*(params.H[35])+rhs[3]*(params.H[48])+rhs[4]*(params.H[61])+rhs[5]*(params.H[74])+rhs[6]*(params.H[87])+rhs[7]*(params.H[100])+rhs[8]*(params.H[113])+rhs[9]*(params.H[126])+rhs[10]*(params.H[139])+rhs[11]*(params.H[152])+rhs[12]*(params.H[165]);
  lhs[10] = rhs[0]*(params.H[10])+rhs[1]*(params.H[23])+rhs[2]*(params.H[36])+rhs[3]*(params.H[49])+rhs[4]*(params.H[62])+rhs[5]*(params.H[75])+rhs[6]*(params.H[88])+rhs[7]*(params.H[101])+rhs[8]*(params.H[114])+rhs[9]*(params.H[127])+rhs[10]*(params.H[140])+rhs[11]*(params.H[153])+rhs[12]*(params.H[166]);
  lhs[11] = rhs[0]*(params.H[11])+rhs[1]*(params.H[24])+rhs[2]*(params.H[37])+rhs[3]*(params.H[50])+rhs[4]*(params.H[63])+rhs[5]*(params.H[76])+rhs[6]*(params.H[89])+rhs[7]*(params.H[102])+rhs[8]*(params.H[115])+rhs[9]*(params.H[128])+rhs[10]*(params.H[141])+rhs[11]*(params.H[154])+rhs[12]*(params.H[167]);
  lhs[12] = rhs[0]*(params.H[12])+rhs[1]*(params.H[25])+rhs[2]*(params.H[38])+rhs[3]*(params.H[51])+rhs[4]*(params.H[64])+rhs[5]*(params.H[77])+rhs[6]*(params.H[90])+rhs[7]*(params.H[103])+rhs[8]*(params.H[116])+rhs[9]*(params.H[129])+rhs[10]*(params.H[142])+rhs[11]*(params.H[155])+rhs[12]*(params.H[168]);
  lhs[13] = 0;
}

void fillq(void) {
  work.q[0] = params.g[0];
  work.q[1] = params.g[1];
  work.q[2] = params.g[2];
  work.q[3] = params.g[3];
  work.q[4] = params.g[4];
  work.q[5] = params.g[5];
  work.q[6] = params.g[6];
  work.q[7] = params.g[7];
  work.q[8] = params.g[8];
  work.q[9] = params.g[9];
  work.q[10] = params.g[10];
  work.q[11] = params.g[11];
  work.q[12] = params.g[12];
  work.q[13] = 0;
}

void fillh(void) {
  work.h[0] = -params.lb[0];
  work.h[1] = -params.lb[1];
  work.h[2] = -params.lb[2];
  work.h[3] = -params.lb[3];
  work.h[4] = -params.lb[4];
  work.h[5] = -params.lb[5];
  work.h[6] = -params.lb[6];
  work.h[7] = -params.lb[7];
  work.h[8] = -params.lb[8];
  work.h[9] = -params.lb[9];
  work.h[10] = -params.lb[10];
  work.h[11] = -params.lb[11];
  work.h[12] = -params.lb[12];
  work.h[13] = params.ub[0];
  work.h[14] = params.ub[1];
  work.h[15] = params.ub[2];
  work.h[16] = params.ub[3];
  work.h[17] = params.ub[4];
  work.h[18] = params.ub[5];
  work.h[19] = params.ub[6];
  work.h[20] = params.ub[7];
  work.h[21] = params.ub[8];
  work.h[22] = params.ub[9];
  work.h[23] = params.ub[10];
  work.h[24] = params.ub[11];
  work.h[25] = params.ub[12];
  work.h[26] = params.d[0];
  work.h[27] = 0;
  work.h[28] = 0;
  work.h[29] = 0;
  work.h[30] = 0;
  work.h[31] = 0;
  work.h[32] = 0;
  work.h[33] = 0;
  work.h[34] = 0;
  work.h[35] = 0;
  work.h[36] = 0;
  work.h[37] = 0;
  work.h[38] = 0;
  work.h[39] = 0;
  work.h[40] = 0;
  work.h[41] = 0;
  work.h[42] = 0;
  work.h[43] = 0;
  work.h[44] = 0;
  work.h[45] = 0;
  work.h[46] = 0;
  work.h[47] = 0;
  work.h[48] = 0;
  work.h[49] = 0;
  work.h[50] = 0;
  work.h[51] = 0;
  work.h[52] = 0;
}

void fillb(void) {
}

void pre_ops(void) {
}
