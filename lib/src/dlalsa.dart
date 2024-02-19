import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlals0.dart';
import 'package:lapack/src/dlasdt.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlalsa(
  final int ICOMPQ,
  final int SMLSIZ,
  final int N,
  final int NRHS,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> BX_,
  final int LDBX,
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
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final B = B_.dim(LDB);
  final BX = BX_.dim(LDBX);
  final U = U_.dim(LDU);
  final K = K_.dim();
  final VT = VT_.dim(LDU);
  final DIFL = DIFL_.dim(LDU);
  final DIFR = DIFR_.dim(LDU);
  final Z = Z_.dim(LDU);
  final POLES = POLES_.dim(LDU);
  final GIVNUM = GIVNUM_.dim(LDU);
  final GIVPTR = GIVPTR_.dim();
  final GIVCOL = GIVCOL_.dim(LDGCOL);
  final PERM = PERM_.dim(LDGCOL);
  final C = C_.dim();
  final S = S_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();

  const ZERO = 0.0, ONE = 1.0;
  int I,
      I1,
      IC,
      IM1,
      INODE,
      J,
      LF,
      LL,
      LVL,
      LVL2,
      NDB1,
      NDIML,
      NDIMR,
      NL,
      NLF,
      NLP1,
      NR,
      NRF,
      NRP1,
      SQRE;
  final NLVL = Box(0), ND = Box(0);

  // Test the input parameters.

  INFO.value = 0;

  if ((ICOMPQ < 0) || (ICOMPQ > 1)) {
    INFO.value = -1;
  } else if (SMLSIZ < 3) {
    INFO.value = -2;
  } else if (N < SMLSIZ) {
    INFO.value = -3;
  } else if (NRHS < 1) {
    INFO.value = -4;
  } else if (LDB < N) {
    INFO.value = -6;
  } else if (LDBX < N) {
    INFO.value = -8;
  } else if (LDU < N) {
    INFO.value = -10;
  } else if (LDGCOL < N) {
    INFO.value = -19;
  }
  if (INFO.value != 0) {
    xerbla('DLALSA', -INFO.value);
    return;
  }

  // Book-keeping and  setting up the computation tree.

  INODE = 1;
  NDIML = INODE + N;
  NDIMR = NDIML + N;

  dlasdt(N, NLVL, ND, IWORK(INODE), IWORK(NDIML), IWORK(NDIMR), SMLSIZ);

  // The following code applies back the left singular vector factors.
  // For applying back the right singular vector factors, go to 50.

  if (ICOMPQ != 1) {
    // The nodes on the bottom level of the tree were solved
    // by DLASDQ. The corresponding left and right singular vector
    // matrices are in explicit form. First apply back the left
    // singular vector matrices.

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
      NR = IWORK[NDIMR + I1];
      NLF = IC - NL;
      NRF = IC + 1;
      dgemm('T', 'N', NL, NRHS, NL, ONE, U(NLF, 1), LDU, B(NLF, 1), LDB, ZERO,
          BX(NLF, 1), LDBX);
      dgemm('T', 'N', NR, NRHS, NR, ONE, U(NRF, 1), LDU, B(NRF, 1), LDB, ZERO,
          BX(NRF, 1), LDBX);
    }

    // Next copy the rows of B that correspond to unchanged rows
    // in the bidiagonal matrix to BX.

    for (I = 1; I <= ND.value; I++) {
      IC = IWORK[INODE + I - 1];
      dcopy(NRHS, B(IC, 1).asArray(), LDB, BX(IC, 1).asArray(), LDBX);
    }

    // Finally go through the left singular vector matrices of all
    // the other subproblems bottom-up on the tree.

    J = pow(2, NLVL.value).toInt();
    SQRE = 0;

    for (LVL = NLVL.value; LVL >= 1; LVL--) {
      LVL2 = 2 * LVL - 1;

      // find the first node LF and last node LL on
      // the current level LVL

      if (LVL == 1) {
        LF = 1;
        LL = 1;
      } else {
        LF = pow(2, (LVL - 1)).toInt();
        LL = 2 * LF - 1;
      }
      for (I = LF; I <= LL; I++) {
        IM1 = I - 1;
        IC = IWORK[INODE + IM1];
        NL = IWORK[NDIML + IM1];
        NR = IWORK[NDIMR + IM1];
        NLF = IC - NL;
        NRF = IC + 1;
        J = J - 1;
        dlals0(
            ICOMPQ,
            NL,
            NR,
            SQRE,
            NRHS,
            BX(NLF, 1),
            LDBX,
            B(NLF, 1),
            LDB,
            PERM(NLF, LVL).asArray(),
            GIVPTR[J],
            GIVCOL(NLF, LVL2),
            LDGCOL,
            GIVNUM(NLF, LVL2),
            LDU,
            POLES(NLF, LVL2),
            DIFL(NLF, LVL).asArray(),
            DIFR(NLF, LVL2),
            Z(NLF, LVL).asArray(),
            K[J],
            C[J],
            S[J],
            WORK,
            INFO);
      }
    }
    return;
  }

  // ICOMPQ = 1: applying back the right singular vector factors.

  // First now go through the right singular vector matrices of all
  // the tree nodes top-down.

  J = 0;
  for (LVL = 1; LVL <= NLVL.value; LVL++) {
    LVL2 = 2 * LVL - 1;

    // Find the first node LF and last node LL on
    // the current level LVL.

    if (LVL == 1) {
      LF = 1;
      LL = 1;
    } else {
      LF = pow(2, (LVL - 1)).toInt();
      LL = 2 * LF - 1;
    }
    for (I = LL; I >= LF; I--) {
      IM1 = I - 1;
      IC = IWORK[INODE + IM1];
      NL = IWORK[NDIML + IM1];
      NR = IWORK[NDIMR + IM1];
      NLF = IC - NL;
      NRF = IC + 1;
      if (I == LL) {
        SQRE = 0;
      } else {
        SQRE = 1;
      }
      J = J + 1;
      dlals0(
          ICOMPQ,
          NL,
          NR,
          SQRE,
          NRHS,
          B(NLF, 1),
          LDB,
          BX(NLF, 1),
          LDBX,
          PERM(NLF, LVL).asArray(),
          GIVPTR[J],
          GIVCOL(NLF, LVL2),
          LDGCOL,
          GIVNUM(NLF, LVL2),
          LDU,
          POLES(NLF, LVL2),
          DIFL(NLF, LVL).asArray(),
          DIFR(NLF, LVL2),
          Z(NLF, LVL).asArray(),
          K[J],
          C[J],
          S[J],
          WORK,
          INFO);
    }
  }

  // The nodes on the bottom level of the tree were solved
  // by DLASDQ. The corresponding right singular vector
  // matrices are in explicit form. Apply them back.

  NDB1 = (ND.value + 1) ~/ 2;
  for (I = NDB1; I <= ND.value; I++) {
    I1 = I - 1;
    IC = IWORK[INODE + I1];
    NL = IWORK[NDIML + I1];
    NR = IWORK[NDIMR + I1];
    NLP1 = NL + 1;
    if (I == ND.value) {
      NRP1 = NR;
    } else {
      NRP1 = NR + 1;
    }
    NLF = IC - NL;
    NRF = IC + 1;
    dgemm('T', 'N', NLP1, NRHS, NLP1, ONE, VT(NLF, 1), LDU, B(NLF, 1), LDB,
        ZERO, BX(NLF, 1), LDBX);
    dgemm('T', 'N', NRP1, NRHS, NRP1, ONE, VT(NRF, 1), LDU, B(NRF, 1), LDB,
        ZERO, BX(NRF, 1), LDBX);
  }
}
