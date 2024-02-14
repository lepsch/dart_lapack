import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlasd6.dart';
import 'package:lapack/src/dlasdq.dart';
import 'package:lapack/src/dlasdt.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlasda(
  final int ICOMPQ,
  final int SMLSIZ,
  final int N,
  final int SQRE,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> VT_,
  final Array<int> K_,
  final Matrix<double> DIFL_,
  final Matrix<double> DIFR_,
  final Matrix<double> Z_,
  final Matrix<double> POLES_,
  final Array<int> GIVPTR_,
  final Matrix<int> GIVCOL_,
  final int LDGCOL,
  final Matrix<int> PERM_,
  final Matrix<double> GIVNUM_,
  final Array<double> C_,
  final Array<double> S_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.dim();
  final E = E_.dim();
  final U = U_.dim(LDU);
  final VT = VT_.dim(LDU);
  final K = K_.dim();
  final DIFL = DIFL_.dim(LDU);
  final DIFR = DIFR_.dim(LDU);
  final Z = Z_.dim(LDU);
  final POLES = POLES_.dim(LDU);
  final GIVPTR = GIVPTR_.dim();
  final GIVCOL = GIVCOL_.dim(LDGCOL);
  final PERM = PERM_.dim(LDGCOL);
  final GIVNUM = GIVNUM_.dim(LDU);
  final C = C_.dim();
  final S = S_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  int I,
      I1,
      IC,
      IDXQ,
      IDXQI,
      IM1,
      INODE,
      ITEMP,
      IWK,
      J,
      LF,
      LL,
      LVL,
      LVL2,
      M,
      NCC,
      NDB1,
      NDIML,
      NDIMR,
      NL,
      NLF,
      NLP1,
      NR,
      NRF,
      NRP1,
      NRU,
      NWORK1,
      NWORK2,
      SMLSZP,
      SQREI,
      VF,
      VFI,
      VL,
      VLI;
  final ALPHA = Box(0.0), BETA = Box(0.0);
  final NLVL = Box(0), ND = Box(0);

  // Test the input parameters.

  INFO.value = 0;

  if ((ICOMPQ < 0) || (ICOMPQ > 1)) {
    INFO.value = -1;
  } else if (SMLSIZ < 3) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if ((SQRE < 0) || (SQRE > 1)) {
    INFO.value = -4;
  } else if (LDU < (N + SQRE)) {
    INFO.value = -8;
  } else if (LDGCOL < N) {
    INFO.value = -17;
  }
  if (INFO.value != 0) {
    xerbla('DLASDA', -INFO.value);
    return;
  }

  M = N + SQRE;

  // If the input matrix is too small, call DLASDQ to find the SVD.

  if (N <= SMLSIZ) {
    if (ICOMPQ == 0) {
      dlasdq('U', SQRE, N, 0, 0, 0, D, E, VT, LDU, U, LDU, U, LDU, WORK, INFO);
    } else {
      dlasdq('U', SQRE, N, M, N, 0, D, E, VT, LDU, U, LDU, U, LDU, WORK, INFO);
    }
    return;
  }

  // Book-keeping and  set up the computation tree.

  INODE = 1;
  NDIML = INODE + N;
  NDIMR = NDIML + N;
  IDXQ = NDIMR + N;
  IWK = IDXQ + N;

  NCC = 0;
  NRU = 0;

  SMLSZP = SMLSIZ + 1;
  VF = 1;
  VL = VF + M;
  NWORK1 = VL + M;
  NWORK2 = NWORK1 + SMLSZP * SMLSZP;

  dlasdt(N, NLVL, ND, IWORK(INODE), IWORK(NDIML), IWORK(NDIMR), SMLSIZ);

  // for the nodes on bottom level of the tree, solve
  // their subproblems by DLASDQ.

  NDB1 = (ND.value + 1) ~/ 2;
  for (I = NDB1; I <= ND.value; I++) {
    // IC : center row of each node
    // NL : number of rows of left  subproblem
    // NR : number of rows of right subproblem
    // NLF: starting row of the left   subproblem
    // NRF: starting row of the right  subproblem

    I1 = I - 1;
    IC = IWORK[INODE + I1];
    NL = IWORK[NDIML + I1];
    NLP1 = NL + 1;
    NR = IWORK[NDIMR + I1];
    NLF = IC - NL;
    NRF = IC + 1;
    IDXQI = IDXQ + NLF - 2;
    VFI = VF + NLF - 1;
    VLI = VL + NLF - 1;
    SQREI = 1;
    if (ICOMPQ == 0) {
      dlaset('A', NLP1, NLP1, ZERO, ONE, WORK(NWORK1).asMatrix(SMLSZP), SMLSZP);
      dlasdq(
          'U',
          SQREI,
          NL,
          NLP1,
          NRU,
          NCC,
          D(NLF),
          E(NLF),
          WORK(NWORK1).asMatrix(SMLSZP),
          SMLSZP,
          WORK(NWORK2).asMatrix(NL),
          NL,
          WORK(NWORK2).asMatrix(NL),
          NL,
          WORK(NWORK2),
          INFO);
      ITEMP = NWORK1 + NL * SMLSZP;
      dcopy(NLP1, WORK(NWORK1), 1, WORK(VFI), 1);
      dcopy(NLP1, WORK(ITEMP), 1, WORK(VLI), 1);
    } else {
      dlaset('A', NL, NL, ZERO, ONE, U(NLF, 1), LDU);
      dlaset('A', NLP1, NLP1, ZERO, ONE, VT(NLF, 1), LDU);
      dlasdq('U', SQREI, NL, NLP1, NL, NCC, D(NLF), E(NLF), VT(NLF, 1), LDU,
          U(NLF, 1), LDU, U(NLF, 1), LDU, WORK(NWORK1), INFO);
      dcopy(NLP1, VT(NLF, 1).asArray(), 1, WORK(VFI), 1);
      dcopy(NLP1, VT(NLF, NLP1).asArray(), 1, WORK(VLI), 1);
    }
    if (INFO.value != 0) {
      return;
    }
    for (J = 1; J <= NL; J++) {
      IWORK[IDXQI + J] = J;
    }
    if ((I == ND.value) && (SQRE == 0)) {
      SQREI = 0;
    } else {
      SQREI = 1;
    }
    IDXQI = IDXQI + NLP1;
    VFI = VFI + NLP1;
    VLI = VLI + NLP1;
    NRP1 = NR + SQREI;
    if (ICOMPQ == 0) {
      dlaset('A', NRP1, NRP1, ZERO, ONE, WORK(NWORK1).asMatrix(SMLSZP), SMLSZP);
      dlasdq(
          'U',
          SQREI,
          NR,
          NRP1,
          NRU,
          NCC,
          D(NRF),
          E(NRF),
          WORK(NWORK1).asMatrix(SMLSZP),
          SMLSZP,
          WORK(NWORK2).asMatrix(NR),
          NR,
          WORK(NWORK2).asMatrix(NR),
          NR,
          WORK(NWORK2),
          INFO);
      ITEMP = NWORK1 + (NRP1 - 1) * SMLSZP;
      dcopy(NRP1, WORK(NWORK1), 1, WORK(VFI), 1);
      dcopy(NRP1, WORK(ITEMP), 1, WORK(VLI), 1);
    } else {
      dlaset('A', NR, NR, ZERO, ONE, U(NRF, 1), LDU);
      dlaset('A', NRP1, NRP1, ZERO, ONE, VT(NRF, 1), LDU);
      dlasdq('U', SQREI, NR, NRP1, NR, NCC, D(NRF), E(NRF), VT(NRF, 1), LDU,
          U(NRF, 1), LDU, U(NRF, 1), LDU, WORK(NWORK1), INFO);
      dcopy(NRP1, VT(NRF, 1).asArray(), 1, WORK(VFI), 1);
      dcopy(NRP1, VT(NRF, NRP1).asArray(), 1, WORK(VLI), 1);
    }
    if (INFO.value != 0) {
      return;
    }
    for (J = 1; J <= NR; J++) {
      IWORK[IDXQI + J] = J;
    }
  }

  // Now conquer each subproblem bottom-up.

  J = pow(2, NLVL.value).toInt();
  for (LVL = NLVL.value; LVL >= 1; LVL--) {
    LVL2 = LVL * 2 - 1;

    // Find the first node LF and last node LL on
    // the current level LVL.

    if (LVL == 1) {
      LF = 1;
      LL = 1;
    } else {
      LF = pow(2, LVL - 1).toInt();
      LL = 2 * LF - 1;
    }
    for (I = LF; I <= LL; I++) {
      IM1 = I - 1;
      IC = IWORK[INODE + IM1];
      NL = IWORK[NDIML + IM1];
      NR = IWORK[NDIMR + IM1];
      NLF = IC - NL;
      NRF = IC + 1;
      if (I == LL) {
        SQREI = SQRE;
      } else {
        SQREI = 1;
      }
      VFI = VF + NLF - 1;
      VLI = VL + NLF - 1;
      IDXQI = IDXQ + NLF - 1;
      ALPHA.value = D[IC];
      BETA.value = E[IC];
      if (ICOMPQ == 0) {
        dlasd6(
            ICOMPQ,
            NL,
            NR,
            SQREI,
            D(NLF),
            WORK(VFI),
            WORK(VLI),
            ALPHA,
            BETA,
            IWORK(IDXQI),
            PERM.asArray(),
            GIVPTR.box(1),
            GIVCOL,
            LDGCOL,
            GIVNUM,
            LDU,
            POLES,
            DIFL.asArray(),
            DIFR.asArray(),
            Z.asArray(),
            K.box(1),
            C.box(1),
            S.box(1),
            WORK(NWORK1),
            IWORK(IWK),
            INFO);
      } else {
        J = J - 1;
        dlasd6(
            ICOMPQ,
            NL,
            NR,
            SQREI,
            D(NLF),
            WORK(VFI),
            WORK(VLI),
            ALPHA,
            BETA,
            IWORK(IDXQI),
            PERM(NLF, LVL).asArray(),
            GIVPTR.box(J),
            GIVCOL(NLF, LVL2),
            LDGCOL,
            GIVNUM(NLF, LVL2),
            LDU,
            POLES(NLF, LVL2),
            DIFL(NLF, LVL).asArray(),
            DIFR(NLF, LVL2).asArray(),
            Z(NLF, LVL).asArray(),
            K.box(J),
            C.box(J),
            S.box(J),
            WORK(NWORK1),
            IWORK(IWK),
            INFO);
      }
      if (INFO.value != 0) {
        return;
      }
    }
  }
}
