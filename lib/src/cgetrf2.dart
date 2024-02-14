import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/slamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlaswp.dart';

void cgetrf2(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim();
  final IPIV = IPIV_.dim();
  double SFMIN;
  Complex TEMP;
  int I, N1, N2;
  final IINFO = Box(0);

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
    xerbla('CGETRF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  if (M == 1) {
    // Use unblocked code for one row case
    // Just need to handle IPIV and INFO.value

    IPIV[1] = 1;
    if (A[1][1] == Complex.zero) INFO.value = 1;
  } else if (N == 1) {
    // Use unblocked code for one column case

    // Compute machine safe minimum

    SFMIN = slamch('S');

    // Find pivot and test for singularity

    I = izamax(M, A[1], 1);
    IPIV[1] = I;
    if (A[I][1] != Complex.zero) {
      // Apply the interchange

      if (I != 1) {
        TEMP = A[1][1];
        A[1][1] = A[I][1];
        A[I][1] = TEMP;
      }

      // Compute elements 2:M of the column

      if ((A[1][1]).abs() >= SFMIN) {
        cscal(M - 1, Complex.one / A[1][1], A[2][1], 1);
      } else {
        for (I = 1; I <= M - 1; I++) {
          // 10
          A[1 + I][1] = A[1 + I][1] / A[1][1];
        } // 10
      }
    } else {
      INFO.value = 1;
    }
  } else {
    // Use recursive code

    N1 = min(M, N) ~/ 2;
    N2 = N - N1;

    // [ A11 ]
    // Factor [ --- ]
    // [ A21 ]

    cgetrf2(M, N1, A, LDA, IPIV, IINFO);
    if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value;

    // [ A12 ]
    // Apply interchanges to [ --- ]
    // [ A22 ]

    claswp(N2, A[1][N1 + 1], LDA, 1, N1, IPIV, 1);

    // Solve A12

    ctrsm('L', 'L', 'N', 'U', N1, N2, Complex.one, A, LDA, A[1][N1 + 1], LDA);

    // Update A22

    cgemm('N', 'N', M - N1, N2, N1, -Complex.one, A[N1 + 1][1], LDA,
        A[1][N1 + 1], LDA, Complex.one, A[N1 + 1][N1 + 1], LDA);

    // Factor A22

    cgetrf2(M - N1, N2, A(N1 + 1,N1 + 1), LDA, IPIV(N1 + 1), IINFO);

    // Adjust INFO.value and the pivot indices

    if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value + N1;
    for (I = N1 + 1; I <= min(M, N); I++) {
      // 20
      IPIV[I] = IPIV[I] + N1;
    } // 20

    // Apply interchanges to A21

    zlaswp(N1, A[1][1], LDA, N1 + 1, min(M, N), IPIV, 1);
  }
}
