// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlasdt.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlals0.dart';

void zlalsa(
  final int ICOMPQ,
  final int SMLSIZ,
  final int N,
  final int NRHS,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> BX_,
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
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final B = B_.having(ld: LDB);
  final BX = BX_.having(ld: LDBX);
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDU);
  final DIFL = DIFL_.having(ld: LDU);
  final DIFR = DIFR_.having(ld: LDU);
  final Z = Z_.having(ld: LDU);
  final POLES = POLES_.having(ld: LDU);
  final GIVNUM = GIVNUM_.having(ld: LDU);
  final C = C_.having();
  final S = S_.having();
  final GIVCOL = GIVCOL_.having(ld: LDGCOL);
  final K = K_.having();
  final GIVPTR = GIVPTR_.having();
  final PERM = PERM_.having(ld: LDGCOL);
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

  const ZERO = 0.0, ONE = 1.0;
  int I,
      I1,
      IC,
      IM1,
      INODE,
      J,
      JCOL,
      JIMAG,
      JREAL,
      JROW,
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
    xerbla('ZLALSA', -INFO.value);
    return;
  }

  // Book-keeping and  setting up the computation tree.

  INODE = 1;
  NDIML = INODE + N;
  NDIMR = NDIML + N;

  dlasdt(N, NLVL, ND, IWORK(INODE), IWORK(NDIML), IWORK(NDIMR), SMLSIZ);

  // The following code applies back the left singular vector factors.
  // For applying back the right singular vector factors, go to 170.

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

      // Since B and BX are complex, the following call to DGEMM
      // is performed in two steps (real and imaginary parts).

      // CALL DGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU,
      // $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )

      J = NL * NRHS * 2;
      for (JCOL = 1; JCOL <= NRHS; JCOL++) {
        for (JROW = NLF; JROW <= NLF + NL - 1; JROW++) {
          J++;
          RWORK[J] = B[JROW][JCOL].real;
        }
      }
      dgemm(
          'T',
          'N',
          NL,
          NRHS,
          NL,
          ONE,
          U(NLF, 1),
          LDU,
          RWORK(1 + NL * NRHS * 2).asMatrix(NL),
          NL,
          ZERO,
          RWORK(1).asMatrix(NL),
          NL);
      J = NL * NRHS * 2;
      for (JCOL = 1; JCOL <= NRHS; JCOL++) {
        for (JROW = NLF; JROW <= NLF + NL - 1; JROW++) {
          J++;
          RWORK[J] = B[JROW][JCOL].imaginary;
        }
      }
      dgemm(
          'T',
          'N',
          NL,
          NRHS,
          NL,
          ONE,
          U(NLF, 1),
          LDU,
          RWORK(1 + NL * NRHS * 2).asMatrix(NL),
          NL,
          ZERO,
          RWORK(1 + NL * NRHS).asMatrix(NL),
          NL);
      JREAL = 0;
      JIMAG = NL * NRHS;
      for (JCOL = 1; JCOL <= NRHS; JCOL++) {
        for (JROW = NLF; JROW <= NLF + NL - 1; JROW++) {
          JREAL++;
          JIMAG++;
          BX[JROW][JCOL] = Complex(RWORK[JREAL], RWORK[JIMAG]);
        }
      }

      // Since B and BX are complex, the following call to DGEMM
      // is performed in two steps (real and imaginary parts).

      // CALL DGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU,
// $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )

      J = NR * NRHS * 2;
      for (JCOL = 1; JCOL <= NRHS; JCOL++) {
        for (JROW = NRF; JROW <= NRF + NR - 1; JROW++) {
          J++;
          RWORK[J] = B[JROW][JCOL].real;
        }
      }
      dgemm(
          'T',
          'N',
          NR,
          NRHS,
          NR,
          ONE,
          U(NRF, 1),
          LDU,
          RWORK(1 + NR * NRHS * 2).asMatrix(NR),
          NR,
          ZERO,
          RWORK(1).asMatrix(NR),
          NR);
      J = NR * NRHS * 2;
      for (JCOL = 1; JCOL <= NRHS; JCOL++) {
        for (JROW = NRF; JROW <= NRF + NR - 1; JROW++) {
          J++;
          RWORK[J] = B[JROW][JCOL].imaginary;
        }
      }
      dgemm(
          'T',
          'N',
          NR,
          NRHS,
          NR,
          ONE,
          U(NRF, 1),
          LDU,
          RWORK(1 + NR * NRHS * 2).asMatrix(NR),
          NR,
          ZERO,
          RWORK(1 + NR * NRHS).asMatrix(NR),
          NR);
      JREAL = 0;
      JIMAG = NR * NRHS;
      for (JCOL = 1; JCOL <= NRHS; JCOL++) {
        for (JROW = NRF; JROW <= NRF + NR - 1; JROW++) {
          JREAL++;
          JIMAG++;
          BX[JROW][JCOL] = Complex(RWORK[JREAL], RWORK[JIMAG]);
        }
      }
    }

    // Next copy the rows of B that correspond to unchanged rows
    // in the bidiagonal matrix to BX.

    for (I = 1; I <= ND.value; I++) {
      IC = IWORK[INODE + I - 1];
      zcopy(NRHS, B(IC, 1).asArray(), LDB, BX(IC, 1).asArray(), LDBX);
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
        J--;
        zlals0(
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
            RWORK,
            INFO);
      }
    }
    return;

    // ICOMPQ = 1: applying back the right singular vector factors.
  }

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
      LF = pow(2, LVL - 1).toInt();
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
      J++;
      zlals0(
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
          RWORK,
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

    // Since B and BX are complex, the following call to DGEMM is
    // performed in two steps (real and imaginary parts).

    // CALL DGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU,
// $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )

    J = NLP1 * NRHS * 2;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      for (JROW = NLF; JROW <= NLF + NLP1 - 1; JROW++) {
        J++;
        RWORK[J] = B[JROW][JCOL].real;
      }
    }
    dgemm(
        'T',
        'N',
        NLP1,
        NRHS,
        NLP1,
        ONE,
        VT(NLF, 1),
        LDU,
        RWORK(1 + NLP1 * NRHS * 2).asMatrix(NLP1),
        NLP1,
        ZERO,
        RWORK(1).asMatrix(NLP1),
        NLP1);
    J = NLP1 * NRHS * 2;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      for (JROW = NLF; JROW <= NLF + NLP1 - 1; JROW++) {
        J++;
        RWORK[J] = B[JROW][JCOL].imaginary;
      }
    }
    dgemm(
        'T',
        'N',
        NLP1,
        NRHS,
        NLP1,
        ONE,
        VT(NLF, 1),
        LDU,
        RWORK(1 + NLP1 * NRHS * 2).asMatrix(NLP1),
        NLP1,
        ZERO,
        RWORK(1 + NLP1 * NRHS).asMatrix(NLP1),
        NLP1);
    JREAL = 0;
    JIMAG = NLP1 * NRHS;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      for (JROW = NLF; JROW <= NLF + NLP1 - 1; JROW++) {
        JREAL++;
        JIMAG++;
        BX[JROW][JCOL] = Complex(RWORK[JREAL], RWORK[JIMAG]);
      }
    }

    // Since B and BX are complex, the following call to DGEMM is
    // performed in two steps (real and imaginary parts).

    // CALL DGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU,
// $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )

    J = NRP1 * NRHS * 2;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      for (JROW = NRF; JROW <= NRF + NRP1 - 1; JROW++) {
        J++;
        RWORK[J] = B[JROW][JCOL].real;
      }
    }
    dgemm(
        'T',
        'N',
        NRP1,
        NRHS,
        NRP1,
        ONE,
        VT(NRF, 1),
        LDU,
        RWORK(1 + NRP1 * NRHS * 2).asMatrix(NRP1),
        NRP1,
        ZERO,
        RWORK(1).asMatrix(NRP1),
        NRP1);
    J = NRP1 * NRHS * 2;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      for (JROW = NRF; JROW <= NRF + NRP1 - 1; JROW++) {
        J++;
        RWORK[J] = B[JROW][JCOL].imaginary;
      }
    }
    dgemm(
        'T',
        'N',
        NRP1,
        NRHS,
        NRP1,
        ONE,
        VT(NRF, 1),
        LDU,
        RWORK(1 + NRP1 * NRHS * 2).asMatrix(NRP1),
        NRP1,
        ZERO,
        RWORK(1 + NRP1 * NRHS).asMatrix(NRP1),
        NRP1);
    JREAL = 0;
    JIMAG = NRP1 * NRHS;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      for (JROW = NRF; JROW <= NRF + NRP1 - 1; JROW++) {
        JREAL++;
        JIMAG++;
        BX[JROW][JCOL] = Complex(RWORK[JREAL], RWORK[JIMAG]);
      }
    }
  }
}
