import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlantb(
  final String NORM,
  final String UPLO,
  final String DIAG,
  final int N,
  final int K,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<double> WORK_,
) {
  final AB = AB_.having(ld: LDAB);
  final WORK = WORK_.having();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  bool UDIAG;
  int I, J, L;
  double VALUE = 0;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N == 0) {
    VALUE = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    if (lsame(DIAG, 'U')) {
      VALUE = ONE;
      if (lsame(UPLO, 'U')) {
        for (J = 1; J <= N; J++) {
          // 20
          for (I = max(K + 2 - J, 1); I <= K; I++) {
            // 10
            SUM.value = AB[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 10
        } // 20
      } else {
        for (J = 1; J <= N; J++) {
          // 40
          for (I = 2; I <= min(N + 1 - J, K + 1); I++) {
            // 30
            SUM.value = AB[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 30
        } // 40
      }
    } else {
      VALUE = ZERO;
      if (lsame(UPLO, 'U')) {
        for (J = 1; J <= N; J++) {
          // 60
          for (I = max(K + 2 - J, 1); I <= K + 1; I++) {
            // 50
            SUM.value = AB[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 50
        } // 60
      } else {
        for (J = 1; J <= N; J++) {
          // 80
          for (I = 1; I <= min(N + 1 - J, K + 1); I++) {
            // 70
            SUM.value = AB[I][J].abs();
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
        if (UDIAG) {
          SUM.value = ONE;
          for (I = max(K + 2 - J, 1); I <= K; I++) {
            // 90
            SUM.value = SUM.value + AB[I][J].abs();
          } // 90
        } else {
          SUM.value = ZERO;
          for (I = max(K + 2 - J, 1); I <= K + 1; I++) {
            // 100
            SUM.value = SUM.value + AB[I][J].abs();
          } // 100
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 110
    } else {
      for (J = 1; J <= N; J++) {
        // 140
        if (UDIAG) {
          SUM.value = ONE;
          for (I = 2; I <= min(N + 1 - J, K + 1); I++) {
            // 120
            SUM.value = SUM.value + AB[I][J].abs();
          } // 120
        } else {
          SUM.value = ZERO;
          for (I = 1; I <= min(N + 1 - J, K + 1); I++) {
            // 130
            SUM.value = SUM.value + AB[I][J].abs();
          } // 130
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 140
    }
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    VALUE = ZERO;
    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= N; I++) {
          // 150
          WORK[I] = ONE;
        } // 150
        for (J = 1; J <= N; J++) {
          // 170
          L = K + 1 - J;
          for (I = max(1, J - K); I <= J - 1; I++) {
            // 160
            WORK[I] += (AB[L + I][J]).abs();
          } // 160
        } // 170
      } else {
        for (I = 1; I <= N; I++) {
          // 180
          WORK[I] = ZERO;
        } // 180
        for (J = 1; J <= N; J++) {
          // 200
          L = K + 1 - J;
          for (I = max(1, J - K); I <= J; I++) {
            // 190
            WORK[I] += (AB[L + I][J]).abs();
          } // 190
        } // 200
      }
    } else {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= N; I++) {
          // 210
          WORK[I] = ONE;
        } // 210
        for (J = 1; J <= N; J++) {
          // 230
          L = 1 - J;
          for (I = J + 1; I <= min(N, J + K); I++) {
            // 220
            WORK[I] += (AB[L + I][J]).abs();
          } // 220
        } // 230
      } else {
        for (I = 1; I <= N; I++) {
          // 240
          WORK[I] = ZERO;
        } // 240
        for (J = 1; J <= N; J++) {
          // 260
          L = 1 - J;
          for (I = J; I <= min(N, J + K); I++) {
            // 250
            WORK[I] += (AB[L + I][J]).abs();
          } // 250
        } // 260
      }
    }
    for (I = 1; I <= N; I++) {
      // 270
      SUM.value = WORK[I];
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    } // 270
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = N.toDouble();
        if (K > 0) {
          for (J = 2; J <= N; J++) {
            // 280
            zlassq(min(J - 1, K), AB(max(K + 2 - J, 1), J).asArray(), 1, SCALE,
                SUM);
          } // 280
        }
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        for (J = 1; J <= N; J++) {
          // 290
          zlassq(
              min(J, K + 1), AB(max(K + 2 - J, 1), J).asArray(), 1, SCALE, SUM);
        } // 290
      }
    } else {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = N.toDouble();
        if (K > 0) {
          for (J = 1; J <= N - 1; J++) {
            // 300
            zlassq(min(N - J, K), AB(2, J).asArray(), 1, SCALE, SUM);
          } // 300
        }
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        for (J = 1; J <= N; J++) {
          // 310
          zlassq(min(N - J + 1, K + 1), AB(1, J).asArray(), 1, SCALE, SUM);
        } // 310
      }
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
