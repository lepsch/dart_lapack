import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlaswp.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgetrf(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine (version 3.X) --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NEGONE = -1.0;
  double SFMIN, TMP;
  int I, J, JP, NSTEP, NTOPIV, NPIVED, KAHEAD;
  int KSTART, IPIVSTART, JPIVSTART, KCOLS = 0;

  // Test the input parameters.

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DGETRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Compute machine safe minimum

  SFMIN = dlamch('S');

  NSTEP = min(M, N);
  for (J = 1; J <= NSTEP; J++) {
    KAHEAD = J & -J;
    KSTART = J + 1 - KAHEAD;
    KCOLS = min(KAHEAD, M - J);

    // Find pivot.

    JP = J - 1 + idamax(M - J + 1, A(J, J).asArray(), 1);
    IPIV[J] = JP;

    // Permute just this column.
    if (JP != J) {
      TMP = A[J][J];
      A[J][J] = A[JP][J];
      A[JP][J] = TMP;
    }

    // Apply pending permutations to L
    NTOPIV = 1;
    IPIVSTART = J;
    JPIVSTART = J - NTOPIV;
    while (NTOPIV < KAHEAD) {
      dlaswp(NTOPIV, A(1, JPIVSTART), LDA, IPIVSTART, J, IPIV, 1);
      IPIVSTART -= NTOPIV;
      NTOPIV *= 2;
      JPIVSTART -= NTOPIV;
    }

    // Permute U block to match L
    dlaswp(KCOLS, A(1, J + 1), LDA, KSTART, J, IPIV, 1);

    // Factor the current column
    if (A[J][J] != ZERO && !disnan(A[J][J])) {
      if (A[J][J].abs() >= SFMIN) {
        dscal(M - J, ONE / A[J][J], A(J + 1, J).asArray(), 1);
      } else {
        for (I = 1; I <= M - J; I++) {
          A[J + I][J] /= A[J][J];
        }
      }
    } else if (A[J][J] == ZERO && INFO.value == 0) {
      INFO.value = J;
    }

    // Solve for U block.
    dtrsm('Left', 'Lower', 'No transpose', 'Unit', KAHEAD, KCOLS, ONE,
        A(KSTART, KSTART), LDA, A(KSTART, J + 1), LDA);
    // Schur complement.
    dgemm(
        'No transpose',
        'No transpose',
        M - J,
        KCOLS,
        KAHEAD,
        NEGONE,
        A(J + 1, KSTART),
        LDA,
        A(KSTART, J + 1),
        LDA,
        ONE,
        A(J + 1, J + 1),
        LDA);
  }

  // Handle pivot permutations on the way out of the recursion
  NPIVED = NSTEP & -NSTEP;
  J = NSTEP - NPIVED;
  while (J > 0) {
    NTOPIV = J & -J;
    dlaswp(NTOPIV, A(1, J - NTOPIV + 1), LDA, J + 1, NSTEP, IPIV, 1);
    J -= NTOPIV;
  }

  // If short and wide, handle the rest of the columns.
  if (M < N) {
    dlaswp(N - M, A(1, M + KCOLS + 1), LDA, 1, M, IPIV, 1);
    dtrsm('Left', 'Lower', 'No transpose', 'Unit', M, N - M, ONE, A, LDA,
        A(1, M + KCOLS + 1), LDA);
  }
}
