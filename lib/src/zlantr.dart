import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlantr(
  final String NORM,
  final String UPLO,
  final String DIAG,
  final int M,
  final int N,
  final Matrix<Complex> A_,
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
          // 20
          for (I = 1; I <= min(M, J - 1); I++) {
            // 10
            SUM.value = A[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 10
        } // 20
      } else {
        for (J = 1; J <= N; J++) {
          // 40
          for (I = J + 1; I <= M; I++) {
            // 30
            SUM.value = A[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 30
        } // 40
      }
    } else {
      VALUE = ZERO;
      if (lsame(UPLO, 'U')) {
        for (J = 1; J <= N; J++) {
          // 60
          for (I = 1; I <= min(M, J); I++) {
            // 50
            SUM.value = A[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 50
        } // 60
      } else {
        for (J = 1; J <= N; J++) {
          // 80
          for (I = J; I <= M; I++) {
            // 70
            SUM.value = A[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 70
        } // 80
      }
    }
  } else if ((lsame(NORM, 'O')) || (NORM == '1')) {
    // Find norm1(A).

    VALUE = ZERO;
    UDIAG = lsame(DIAG, 'U');
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        // 110
        if ((UDIAG) && (J <= M)) {
          SUM.value = ONE;
          for (I = 1; I <= J - 1; I++) {
            // 90
            SUM.value = SUM.value + A[I][J].abs();
          } // 90
        } else {
          SUM.value = ZERO;
          for (I = 1; I <= min(M, J); I++) {
            // 100
            SUM.value = SUM.value + A[I][J].abs();
          } // 100
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 110
    } else {
      for (J = 1; J <= N; J++) {
        // 140
        if (UDIAG) {
          SUM.value = ONE;
          for (I = J + 1; I <= M; I++) {
            // 120
            SUM.value = SUM.value + A[I][J].abs();
          } // 120
        } else {
          SUM.value = ZERO;
          for (I = J; I <= M; I++) {
            // 130
            SUM.value = SUM.value + A[I][J].abs();
          } // 130
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 140
    }
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= M; I++) {
          // 150
          WORK[I] = ONE;
        } // 150
        for (J = 1; J <= N; J++) {
          // 170
          for (I = 1; I <= min(M, J - 1); I++) {
            // 160
            WORK[I] = WORK[I] + A[I][J].abs();
          } // 160
        } // 170
      } else {
        for (I = 1; I <= M; I++) {
          // 180
          WORK[I] = ZERO;
        } // 180
        for (J = 1; J <= N; J++) {
          // 200
          for (I = 1; I <= min(M, J); I++) {
            // 190
            WORK[I] = WORK[I] + A[I][J].abs();
          } // 190
        } // 200
      }
    } else {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= min(M, N); I++) {
          // 210
          WORK[I] = ONE;
        } // 210
        for (I = N + 1; I <= M; I++) {
          // 220
          WORK[I] = ZERO;
        } // 220
        for (J = 1; J <= N; J++) {
          // 240
          for (I = J + 1; I <= M; I++) {
            // 230
            WORK[I] = WORK[I] + A[I][J].abs();
          } // 230
        } // 240
      } else {
        for (I = 1; I <= M; I++) {
          // 250
          WORK[I] = ZERO;
        } // 250
        for (J = 1; J <= N; J++) {
          // 270
          for (I = J; I <= M; I++) {
            // 260
            WORK[I] = WORK[I] + A[I][J].abs();
          } // 260
        } // 270
      }
    }
    VALUE = ZERO;
    for (I = 1; I <= M; I++) {
      // 280
      SUM.value = WORK[I];
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    } // 280
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = min(M, N).toDouble();
        for (J = 2; J <= N; J++) {
          // 290
          zlassq(min(M, J - 1), A(1, J).asArray(), 1, SCALE, SUM);
        } // 290
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        for (J = 1; J <= N; J++) {
          // 300
          zlassq(min(M, J), A(1, J).asArray(), 1, SCALE, SUM);
        } // 300
      }
    } else {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = min(M, N).toDouble();
        for (J = 1; J <= N; J++) {
          // 310
          zlassq(M - J, A(min(M, J + 1), J).asArray(), 1, SCALE, SUM);
        } // 310
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        for (J = 1; J <= N; J++) {
          // 320
          zlassq(M - J + 1, A(J, J).asArray(), 1, SCALE, SUM);
        } // 320
      }
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
