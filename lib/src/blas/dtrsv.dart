import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dtrsv(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Array<double> X,
  final int INCX,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  double TEMP;
  int I, INFO, IX, J, JX, KX = 0;
  bool NOUNIT;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO = 2;
  } else if (!lsame(DIAG, 'U') && !lsame(DIAG, 'N')) {
    INFO = 3;
  } else if (N < 0) {
    INFO = 4;
  } else if (LDA < max(1, N)) {
    INFO = 6;
  } else if (INCX == 0) {
    INFO = 8;
  }
  if (INFO != 0) {
    xerbla('DTRSV ', INFO);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  NOUNIT = lsame(DIAG, 'N');

  // Set up the start point in X if the increment is not unity. This
  // will be  ( N - 1 )*INCX  too small for descending loops.

  if (INCX <= 0) {
    KX = 1 - (N - 1) * INCX;
  } else if (INCX != 1) {
    KX = 1;
  }

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through A.

  if (lsame(TRANS, 'N')) {
    // Form  x := inv( A )*x.

    if (lsame(UPLO, 'U')) {
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          if (X[J] != ZERO) {
            if (NOUNIT) X[J] = X[J] / A[J][J];
            TEMP = X[J];
            for (I = J - 1; I >= 1; I--) {
              X[I] = X[I] - TEMP * A[I][J];
            }
          }
        }
      } else {
        JX = KX + (N - 1) * INCX;
        for (J = N; J >= 1; J--) {
          if (X[JX] != ZERO) {
            if (NOUNIT) X[JX] = X[JX] / A[J][J];
            TEMP = X[JX];
            IX = JX;
            for (I = J - 1; I >= 1; I--) {
              IX = IX - INCX;
              X[IX] = X[IX] - TEMP * A[I][J];
            }
          }
          JX = JX - INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          if (X[J] != ZERO) {
            if (NOUNIT) X[J] = X[J] / A[J][J];
            TEMP = X[J];
            for (I = J + 1; I <= N; I++) {
              X[I] = X[I] - TEMP * A[I][J];
            }
          }
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          if (X[JX] != ZERO) {
            if (NOUNIT) X[JX] = X[JX] / A[J][J];
            TEMP = X[JX];
            IX = JX;
            for (I = J + 1; I <= N; I++) {
              IX = IX + INCX;
              X[IX] = X[IX] - TEMP * A[I][J];
            }
          }
          JX = JX + INCX;
        }
      }
    }
  } else {
    // Form  x := inv( A**T )*x.

    if (lsame(UPLO, 'U')) {
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          TEMP = X[J];
          for (I = 1; I <= J - 1; I++) {
            TEMP = TEMP - A[I][J] * X[I];
          }
          if (NOUNIT) TEMP = TEMP / A[J][J];
          X[J] = TEMP;
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          TEMP = X[JX];
          IX = KX;
          for (I = 1; I <= J - 1; I++) {
            TEMP = TEMP - A[I][J] * X[IX];
            IX = IX + INCX;
          }
          if (NOUNIT) TEMP = TEMP / A[J][J];
          X[JX] = TEMP;
          JX = JX + INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          TEMP = X[J];
          for (I = N; I >= J + 1; I--) {
            TEMP = TEMP - A[I][J] * X[I];
          }
          if (NOUNIT) TEMP = TEMP / A[J][J];
          X[J] = TEMP;
        }
      } else {
        KX = KX + (N - 1) * INCX;
        JX = KX;
        for (J = N; J >= 1; J--) {
          TEMP = X[JX];
          IX = KX;
          for (I = N; I >= J + 1; I--) {
            TEMP = TEMP - A[I][J] * X[IX];
            IX = IX - INCX;
          }
          if (NOUNIT) TEMP = TEMP / A[J][J];
          X[JX] = TEMP;
          JX = JX - INCX;
        }
      }
    }
  }
}
