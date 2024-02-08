import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlaeda(
  final int N,
  final int TLVLS,
  final int CURLVL,
  final int CURPBM,
  final Array<int> PRMPTR,
  final Array<int> PERM,
  final Array<int> GIVPTR,
  final Matrix<int> GIVCOL,
  final Matrix<double> GIVNUM,
  final Array<double> Q,
  final Array<int> QPTR,
  final Array<double> Z,
  final Array<double> ZTEMP,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  int BSIZ1, BSIZ2, CURR, I, K, MID, PSIZ1, PSIZ2, PTR, ZPTR1;

  // Test the input parameters.

  INFO.value = 0;

  if (N < 0) {
    INFO.value = -1;
  }
  if (INFO.value != 0) {
    xerbla('DLAEDA', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Determine location of first number in second half.

  MID = N ~/ 2 + 1;

  // Gather last/first rows of appropriate eigenblocks into center of Z

  PTR = 1;

  // Determine location of lowest level subproblem in the full storage
  // scheme

  CURR =
      PTR + CURPBM * pow(2, CURLVL).toInt() + pow(2, (CURLVL - 1)).toInt() - 1;

  // Determine size of these matrices.  We add HALF to the value of
  // the SQRT in case the machine underestimates one of these square
  // roots.

  BSIZ1 = (HALF + sqrt(QPTR[CURR + 1] - QPTR[CURR])).toInt();
  BSIZ2 = (HALF + sqrt(QPTR[CURR + 2] - QPTR[CURR + 1])).toInt();
  for (K = 1; K <= MID - BSIZ1 - 1; K++) {
    Z[K] = ZERO;
  }
  dcopy(BSIZ1, Q(QPTR[CURR] + BSIZ1 - 1), BSIZ1, Z(MID - BSIZ1), 1);
  dcopy(BSIZ2, Q(QPTR[CURR + 1]), BSIZ2, Z(MID), 1);
  for (K = MID + BSIZ2; K <= N; K++) {
    Z[K] = ZERO;
  }

  // Loop through remaining levels 1 -> CURLVL applying the Givens
  // rotations and permutation and then multiplying the center matrices
  // against the current Z.

  PTR = pow(2, TLVLS).toInt() + 1;
  for (K = 1; K <= CURLVL - 1; K++) {
    CURR = PTR +
        CURPBM * pow(2, (CURLVL - K)).toInt() +
        pow(2, (CURLVL - K - 1)).toInt() -
        1;
    PSIZ1 = PRMPTR[CURR + 1] - PRMPTR[CURR];
    PSIZ2 = PRMPTR[CURR + 2] - PRMPTR[CURR + 1];
    ZPTR1 = MID - PSIZ1;

    // Apply Givens at CURR and CURR+1

    for (I = GIVPTR[CURR]; I <= GIVPTR[CURR + 1] - 1; I++) {
      drot(
        1,
        Z(ZPTR1 + GIVCOL[1][I] - 1),
        1,
        Z(ZPTR1 + GIVCOL[2][I] - 1),
        1,
        GIVNUM[1][I],
        GIVNUM[2][I],
      );
    }
    for (I = GIVPTR[CURR + 1]; I <= GIVPTR[CURR + 2] - 1; I++) {
      drot(
        1,
        Z(MID - 1 + GIVCOL[1][I]),
        1,
        Z(MID - 1 + GIVCOL[2][I]),
        1,
        GIVNUM[1][I],
        GIVNUM[2][I],
      );
    }
    PSIZ1 = PRMPTR[CURR + 1] - PRMPTR[CURR];
    PSIZ2 = PRMPTR[CURR + 2] - PRMPTR[CURR + 1];
    for (I = 0; I <= PSIZ1 - 1; I++) {
      ZTEMP[I + 1] = Z[ZPTR1 + PERM[PRMPTR[CURR] + I] - 1];
    }
    for (I = 0; I <= PSIZ2 - 1; I++) {
      ZTEMP[PSIZ1 + I + 1] = Z[MID + PERM[PRMPTR[CURR + 1] + I] - 1];
    }

    // Multiply Blocks at CURR and CURR+1

    // Determine size of these matrices.  We add HALF to the value of
    // the SQRT in case the machine underestimates one of these
    // square roots.

    BSIZ1 = (HALF + sqrt(QPTR[CURR + 1] - QPTR[CURR])).toInt();
    BSIZ2 = (HALF + sqrt(QPTR[CURR + 2] - QPTR[CURR + 1])).toInt();
    if (BSIZ1 > 0) {
      dgemv(
        'T',
        BSIZ1,
        BSIZ1,
        ONE,
        Q(QPTR[CURR]).asMatrix(BSIZ1),
        BSIZ1,
        ZTEMP(1),
        1,
        ZERO,
        Z(ZPTR1),
        1,
      );
    }
    dcopy(PSIZ1 - BSIZ1, ZTEMP(BSIZ1 + 1), 1, Z(ZPTR1 + BSIZ1), 1);
    if (BSIZ2 > 0) {
      dgemv(
        'T',
        BSIZ2,
        BSIZ2,
        ONE,
        Q(QPTR[CURR + 1]).asMatrix(BSIZ2),
        BSIZ2,
        ZTEMP(PSIZ1 + 1),
        1,
        ZERO,
        Z(MID),
        1,
      );
    }
    dcopy(PSIZ2 - BSIZ2, ZTEMP(PSIZ1 + BSIZ2 + 1), 1, Z(MID + BSIZ2), 1);

    PTR = PTR + pow(2, (TLVLS - K)).toInt();
  }
}
