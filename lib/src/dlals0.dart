import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlals0(
  final int ICOMPQ,
  final int NL,
  final int NR,
  final int SQRE,
  final int NRHS,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> BX_,
  final int LDBX,
  final Array<int> PERM_,
  final int GIVPTR,
  final Matrix<int> GIVCOL_,
  final int LDGCOL,
  final Matrix<double> GIVNUM_,
  final int LDGNUM,
  final Matrix<double> POLES_,
  final Array<double> DIFL_,
  final Matrix<double> DIFR_,
  final Array<double> Z_,
  final int K,
  final double C,
  final double S,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final B = B_.dim(LDB);
  final BX = BX_.dim(LDBX);
  final PERM = PERM_.dim();
  final GIVCOL = GIVCOL_.dim(LDGCOL);
  final GIVNUM = GIVNUM_.dim(LDGNUM);
  final POLES = POLES_.dim(LDGNUM);
  final DIFR = DIFR_.dim(LDGNUM);
  final DIFL = DIFL_.dim();
  final Z = Z_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0, NEGONE = -1.0;
  int I, J, M, N, NLP1;
  double DIFLJ, DIFRJ = 0, DJ, DSIGJ, DSIGJP = 0, TEMP;

  // Test the input parameters.

  INFO.value = 0;
  N = NL + NR + 1;

  if ((ICOMPQ < 0) || (ICOMPQ > 1)) {
    INFO.value = -1;
  } else if (NL < 1) {
    INFO.value = -2;
  } else if (NR < 1) {
    INFO.value = -3;
  } else if ((SQRE < 0) || (SQRE > 1)) {
    INFO.value = -4;
  } else if (NRHS < 1) {
    INFO.value = -5;
  } else if (LDB < N) {
    INFO.value = -7;
  } else if (LDBX < N) {
    INFO.value = -9;
  } else if (GIVPTR < 0) {
    INFO.value = -11;
  } else if (LDGCOL < N) {
    INFO.value = -13;
  } else if (LDGNUM < N) {
    INFO.value = -15;
  } else if (K < 1) {
    INFO.value = -20;
  }
  if (INFO.value != 0) {
    xerbla('DLALS0', -INFO.value);
    return;
  }

  M = N + SQRE;
  NLP1 = NL + 1;

  if (ICOMPQ == 0) {
    // Apply back orthogonal transformations from the left.

    // Step (1L): apply back the Givens rotations performed.

    for (I = 1; I <= GIVPTR; I++) {
      drot(NRHS, B(GIVCOL[I][2], 1).asArray(), LDB,
          B(GIVCOL[I][1], 1).asArray(), LDB, GIVNUM[I][2], GIVNUM[I][1]);
    }

    // Step (2L): permute rows of B.

    dcopy(NRHS, B(NLP1, 1).asArray(), LDB, BX(1, 1).asArray(), LDBX);
    for (I = 2; I <= N; I++) {
      dcopy(NRHS, B(PERM[I], 1).asArray(), LDB, BX(I, 1).asArray(), LDBX);
    }

    // Step (3L): apply the inverse of the left singular vector
    // matrix to BX.

    if (K == 1) {
      dcopy(NRHS, BX.asArray(), LDBX, B.asArray(), LDB);
      if (Z[1] < ZERO) {
        dscal(NRHS, NEGONE, B.asArray(), LDB);
      }
    } else {
      for (J = 1; J <= K; J++) {
        DIFLJ = DIFL[J];
        DJ = POLES[J][1];
        DSIGJ = -POLES[J][2];
        if (J < K) {
          DIFRJ = -DIFR[J][1];
          DSIGJP = -POLES[J + 1][2];
        }
        if ((Z[J] == ZERO) || (POLES[J][2] == ZERO)) {
          WORK[J] = ZERO;
        } else {
          WORK[J] = -POLES[J][2] * Z[J] / DIFLJ / (POLES[J][2] + DJ);
        }
        for (I = 1; I <= J - 1; I++) {
          if ((Z[I] == ZERO) || (POLES[I][2] == ZERO)) {
            WORK[I] = ZERO;
          } else {
            // Use calls to the subroutine DLAMC3 to enforce the
            // parentheses (x+y)+z. The goal is to prevent
            // optimizing compilers from doing x+(y+z).

            WORK[I] = POLES[I][2] *
                Z[I] /
                (dlamc3(POLES[I][2], DSIGJ) - DIFLJ) /
                (POLES[I][2] + DJ);
          }
        }
        for (I = J + 1; I <= K; I++) {
          if ((Z[I] == ZERO) || (POLES[I][2] == ZERO)) {
            WORK[I] = ZERO;
          } else {
            WORK[I] = POLES[I][2] *
                Z[I] /
                (dlamc3(POLES[I][2], DSIGJP) + DIFRJ) /
                (POLES[I][2] + DJ);
          }
        }
        WORK[1] = NEGONE;
        TEMP = dnrm2(K, WORK, 1);
        dgemv(
            'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO, B(J, 1).asArray(), LDB);
        dlascl('G', 0, 0, TEMP, ONE, 1, NRHS, B(J, 1), LDB, INFO);
      }
    }

    // Move the deflated rows of BX to B also.

    if (K < max(M, N)) {
      dlacpy('A', N - K, NRHS, BX(K + 1, 1), LDBX, B(K + 1, 1), LDB);
    }
  } else {
    // Apply back the right orthogonal transformations.

    // Step (1R): apply back the new right singular vector matrix
    // to B.

    if (K == 1) {
      dcopy(NRHS, B.asArray(), LDB, BX.asArray(), LDBX);
    } else {
      for (J = 1; J <= K; J++) {
        DSIGJ = POLES[J][2];
        if (Z[J] == ZERO) {
          WORK[J] = ZERO;
        } else {
          WORK[J] = -Z[J] / DIFL[J] / (DSIGJ + POLES[J][1]) / DIFR[J][2];
        }
        for (I = 1; I <= J - 1; I++) {
          if (Z[J] == ZERO) {
            WORK[I] = ZERO;
          } else {
            // Use calls to the subroutine DLAMC3 to enforce the
            // parentheses (x+y)+z. The goal is to prevent
            // optimizing compilers from doing x+(y+z).

            WORK[I] = Z[J] /
                (dlamc3(DSIGJ, -POLES[I + 1][2]) - DIFR[I][1]) /
                (DSIGJ + POLES[I][1]) /
                DIFR[I][2];
          }
        }
        for (I = J + 1; I <= K; I++) {
          if (Z[J] == ZERO) {
            WORK[I] = ZERO;
          } else {
            WORK[I] = Z[J] /
                (dlamc3(DSIGJ, -POLES[I][2]) - DIFL[I]) /
                (DSIGJ + POLES[I][1]) /
                DIFR[I][2];
          }
        }
        dgemv(
            'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO, BX(J, 1).asArray(), LDBX);
      }
    }

    // Step (2R): if SQRE = 1, apply back the rotation that is
    // related to the right null space of the subproblem.

    if (SQRE == 1) {
      dcopy(NRHS, B(M, 1).asArray(), LDB, BX(M, 1).asArray(), LDBX);
      drot(NRHS, BX(1, 1).asArray(), LDBX, BX(M, 1).asArray(), LDBX, C, S);
    }
    if (K < max(M, N)) {
      dlacpy('A', N - K, NRHS, B(K + 1, 1), LDB, BX(K + 1, 1), LDBX);
    }

    // Step (3R): permute rows of B.

    dcopy(NRHS, BX(1, 1).asArray(), LDBX, B(NLP1, 1).asArray(), LDB);
    if (SQRE == 1) {
      dcopy(NRHS, BX(M, 1).asArray(), LDBX, B(M, 1).asArray(), LDB);
    }
    for (I = 2; I <= N; I++) {
      dcopy(NRHS, BX(I, 1).asArray(), LDBX, B(PERM[I], 1).asArray(), LDB);
    }

    // Step (4R): apply back the Givens rotations performed.

    for (I = GIVPTR; I >= 1; I--) {
      drot(NRHS, B(GIVCOL[I][2], 1).asArray(), LDB,
          B(GIVCOL[I][1], 1).asArray(), LDB, GIVNUM[I][2], -GIVNUM[I][1]);
    }
  }
}
