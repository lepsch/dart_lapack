import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void ztrmv(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> X_,
  final int INCX,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having();

  // Test the input parameters.

  var INFO = 0;
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
    xerbla('ZTRMV ', INFO);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  final NOCONJ = lsame(TRANS, 'T');
  final NOUNIT = lsame(DIAG, 'N');

  // Set up the start point in X if the increment is not unity. This
  // will be  ( N - 1 )*INCX  too small for descending loops.

  var KX = switch (INCX) {
    <= 0 => 1 - (N - 1) * INCX,
    1 => 0,
    _ => 1,
  };

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through A.

  if (lsame(TRANS, 'N')) {
    // Form  x := A*x.

    if (lsame(UPLO, 'U')) {
      if (INCX == 1) {
        for (var J = 1; J <= N; J++) {
          if (X[J] != Complex.zero) {
            final TEMP = X[J];
            for (var I = 1; I <= J - 1; I++) {
              X[I] += TEMP * A[I][J];
            }
            if (NOUNIT) X[J] *= A[J][J];
          }
        }
      } else {
        var JX = KX;
        for (var J = 1; J <= N; J++) {
          if (X[JX] != Complex.zero) {
            final TEMP = X[JX];
            var IX = KX;
            for (var I = 1; I <= J - 1; I++) {
              X[IX] += TEMP * A[I][J];
              IX += INCX;
            }
            if (NOUNIT) X[JX] *= A[J][J];
          }
          JX += INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (var J = N; J >= 1; J--) {
          if (X[J] != Complex.zero) {
            final TEMP = X[J];
            for (var I = N; I >= J + 1; I--) {
              X[I] += TEMP * A[I][J];
            }
            if (NOUNIT) X[J] *= A[J][J];
          }
        }
      } else {
        KX += (N - 1) * INCX;
        var JX = KX;
        for (var J = N; J >= 1; J--) {
          if (X[JX] != Complex.zero) {
            final TEMP = X[JX];
            var IX = KX;
            for (var I = N; I >= J + 1; I--) {
              X[IX] += TEMP * A[I][J];
              IX -= INCX;
            }
            if (NOUNIT) X[JX] *= A[J][J];
          }
          JX -= INCX;
        }
      }
    }
  } else {
    // Form  x := A**T*x  or  x := A**H*x.

    if (lsame(UPLO, 'U')) {
      if (INCX == 1) {
        for (var J = N; J >= 1; J--) {
          var TEMP = X[J];
          if (NOCONJ) {
            if (NOUNIT) TEMP *= A[J][J];
            for (var I = J - 1; I >= 1; I--) {
              TEMP += A[I][J] * X[I];
            }
          } else {
            if (NOUNIT) TEMP *= A[J][J].conjugate();
            for (var I = J - 1; I >= 1; I--) {
              TEMP += A[I][J].conjugate() * X[I];
            }
          }
          X[J] = TEMP;
        }
      } else {
        var JX = KX + (N - 1) * INCX;
        for (var J = N; J >= 1; J--) {
          var TEMP = X[JX];
          var IX = JX;
          if (NOCONJ) {
            if (NOUNIT) TEMP *= A[J][J];
            for (var I = J - 1; I >= 1; I--) {
              IX -= INCX;
              TEMP += A[I][J] * X[IX];
            }
          } else {
            if (NOUNIT) TEMP *= A[J][J].conjugate();
            for (var I = J - 1; I >= 1; I--) {
              IX -= INCX;
              TEMP += A[I][J].conjugate() * X[IX];
            }
          }
          X[JX] = TEMP;
          JX -= INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (var J = 1; J <= N; J++) {
          var TEMP = X[J];
          if (NOCONJ) {
            if (NOUNIT) TEMP *= A[J][J];
            for (var I = J + 1; I <= N; I++) {
              TEMP += A[I][J] * X[I];
            }
          } else {
            if (NOUNIT) TEMP *= A[J][J].conjugate();
            for (var I = J + 1; I <= N; I++) {
              TEMP += A[I][J].conjugate() * X[I];
            }
          }
          X[J] = TEMP;
        }
      } else {
        var JX = KX;
        for (var J = 1; J <= N; J++) {
          var TEMP = X[JX];
          var IX = JX;
          if (NOCONJ) {
            if (NOUNIT) TEMP *= A[J][J];
            for (var I = J + 1; I <= N; I++) {
              IX += INCX;
              TEMP += A[I][J] * X[IX];
            }
          } else {
            if (NOUNIT) TEMP *= A[J][J].conjugate();
            for (var I = J + 1; I <= N; I++) {
              IX += INCX;
              TEMP += A[I][J].conjugate() * X[IX];
            }
          }
          X[JX] = TEMP;
          JX += INCX;
        }
      }
    }
  }
}
