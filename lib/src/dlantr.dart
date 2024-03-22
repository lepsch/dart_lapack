import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/matrix.dart';

double dlantr(
  final String NORM,
  final String UPLO,
  final String DIAG,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> WORK_,
) {
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  bool UDIAG;
  int I, J;
  double VALUE = 0;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (min(M, N) == 0) {
    VALUE = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    if (lsame(DIAG, 'U')) {
      VALUE = ONE;
      if (lsame(UPLO, 'U')) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= min(M, J - 1); I++) {
            SUM.value = A[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = J + 1; I <= M; I++) {
            SUM.value = A[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
        }
      }
    } else {
      VALUE = ZERO;
      if (lsame(UPLO, 'U')) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= min(M, J); I++) {
            SUM.value = A[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = J; I <= M; I++) {
            SUM.value = A[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
        }
      }
    }
  } else if ((lsame(NORM, 'O')) || (NORM == '1')) {
    // Find norm1(A).

    VALUE = ZERO;
    UDIAG = lsame(DIAG, 'U');
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        if (UDIAG && (J <= M)) {
          SUM.value = ONE;
          for (I = 1; I <= J - 1; I++) {
            SUM.value += A[I][J].abs();
          }
        } else {
          SUM.value = ZERO;
          for (I = 1; I <= min(M, J); I++) {
            SUM.value += A[I][J].abs();
          }
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    } else {
      for (J = 1; J <= N; J++) {
        if (UDIAG) {
          SUM.value = ONE;
          for (I = J + 1; I <= M; I++) {
            SUM.value += A[I][J].abs();
          }
        } else {
          SUM.value = ZERO;
          for (I = J; I <= M; I++) {
            SUM.value += A[I][J].abs();
          }
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    }
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= M; I++) {
          WORK[I] = ONE;
        }
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= min(M, J - 1); I++) {
            WORK[I] += A[I][J].abs();
          }
        }
      } else {
        for (I = 1; I <= M; I++) {
          WORK[I] = ZERO;
        }
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= min(M, J); I++) {
            WORK[I] += A[I][J].abs();
          }
        }
      }
    } else {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= min(M, N); I++) {
          WORK[I] = ONE;
        }
        for (I = N + 1; I <= M; I++) {
          WORK[I] = ZERO;
        }
        for (J = 1; J <= N; J++) {
          for (I = J + 1; I <= M; I++) {
            WORK[I] += A[I][J].abs();
          }
        }
      } else {
        for (I = 1; I <= M; I++) {
          WORK[I] = ZERO;
        }
        for (J = 1; J <= N; J++) {
          for (I = J; I <= M; I++) {
            WORK[I] += A[I][J].abs();
          }
        }
      }
    }
    VALUE = ZERO;
    for (I = 1; I <= M; I++) {
      SUM.value = WORK[I];
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    }
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = min(M, N).toDouble();
        for (J = 2; J <= N; J++) {
          dlassq(min(M, J - 1), A(1, J).asArray(), 1, SCALE, SUM);
        }
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        for (J = 1; J <= N; J++) {
          dlassq(min(M, J), A(1, J).asArray(), 1, SCALE, SUM);
        }
      }
    } else {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = min(M, N).toDouble();
        for (J = 1; J <= N; J++) {
          dlassq(M - J, A(min(M, J + 1), J).asArray(), 1, SCALE, SUM);
        }
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        for (J = 1; J <= N; J++) {
          dlassq(M - J + 1, A(J, J).asArray(), 1, SCALE, SUM);
        }
      }
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
