import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaswp.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgetrf2(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Test the input parameters

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DGETRF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  if (M == 1) {
    // Use unblocked code for one row case
    // Just need to handle IPIV and INFO

    IPIV[1] = 1;
    if (A[1][1] == ZERO) INFO.value = 1;
  } else if (N == 1) {
    // Use unblocked code for one column case

    // Compute machine safe minimum

    final SFMIN = dlamch('S');

    // Find pivot and test for singularity

    final I = idamax(M, A(1, 1).asArray(), 1);
    IPIV[1] = I;
    if (A[I][1] != ZERO) {
      // Apply the interchange

      if (I != 1) {
        final TEMP = A[1][1];
        A[1][1] = A[I][1];
        A[I][1] = TEMP;
      }

      // Compute elements 2:M of the column

      if (A[1][1].abs() >= SFMIN) {
        dscal(M - 1, ONE / A[1][1], A(2, 1).asArray(), 1);
      } else {
        for (var I = 1; I <= M - 1; I++) {
          A[1 + I][1] /= A[1][1];
        }
      }
    } else {
      INFO.value = 1;
    }
  } else {
    // Use recursive code

    final N1 = min(M, N) ~/ 2;
    final N2 = N - N1;

    //        [ A11 ]
    // Factor [ --- ]
    //        [ A21 ]

    final IINFO = Box(0);
    dgetrf2(M, N1, A, LDA, IPIV, IINFO);
    if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value;

    //                       [ A12 ]
    // Apply interchanges to [ --- ]
    //                       [ A22 ]

    dlaswp(N2, A(1, N1 + 1), LDA, 1, N1, IPIV, 1);

    // Solve A12

    dtrsm('L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, A(1, N1 + 1), LDA);

    // Update A22

    dgemm('N', 'N', M - N1, N2, N1, -ONE, A(N1 + 1, 1), LDA, A(1, N1 + 1), LDA,
        ONE, A(N1 + 1, N1 + 1), LDA);

    // Factor A22

    dgetrf2(M - N1, N2, A(N1 + 1, N1 + 1), LDA, IPIV(N1 + 1), IINFO);

    // Adjust INFO and the pivot indices

    if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value + N1;
    for (var I = N1 + 1; I <= min(M, N); I++) {
      IPIV[I] += N1;
    }

    // Apply interchanges to A21

    dlaswp(N1, A(1, 1), LDA, N1 + 1, min(M, N), IPIV, 1);
  }
}
