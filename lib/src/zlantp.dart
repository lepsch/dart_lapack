import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlantp(
  final String NORM,
  final String UPLO,
  final String DIAG,
  final int N,
  final Array<Complex> AP_,
  final Array<double> WORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UDIAG;
  int I, J, K = 0;
  double VALUE = 0;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N == 0) {
    VALUE = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    K = 1;
    if (lsame(DIAG, 'U')) {
      VALUE = ONE;
      if (lsame(UPLO, 'U')) {
        for (J = 1; J <= N; J++) {
          // 20
          for (I = K; I <= K + J - 2; I++) {
            // 10
            SUM.value = AP[I].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 10
          K += J;
        } // 20
      } else {
        for (J = 1; J <= N; J++) {
          // 40
          for (I = K + 1; I <= K + N - J; I++) {
            // 30
            SUM.value = AP[I].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 30
          K += N - J + 1;
        } // 40
      }
    } else {
      VALUE = ZERO;
      if (lsame(UPLO, 'U')) {
        for (J = 1; J <= N; J++) {
          // 60
          for (I = K; I <= K + J - 1; I++) {
            // 50
            SUM.value = AP[I].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 50
          K += J;
        } // 60
      } else {
        for (J = 1; J <= N; J++) {
          // 80
          for (I = K; I <= K + N - J; I++) {
            // 70
            SUM.value = AP[I].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          } // 70
          K += N - J + 1;
        } // 80
      }
    }
  } else if ((lsame(NORM, 'O')) || (NORM == '1')) {
    // Find norm1(A).

    VALUE = ZERO;
    K = 1;
    UDIAG = lsame(DIAG, 'U');
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        // 110
        if (UDIAG) {
          SUM.value = ONE;
          for (I = K; I <= K + J - 2; I++) {
            // 90
            SUM.value = SUM.value + AP[I].abs();
          } // 90
        } else {
          SUM.value = ZERO;
          for (I = K; I <= K + J - 1; I++) {
            // 100
            SUM.value = SUM.value + AP[I].abs();
          } // 100
        }
        K += J;
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 110
    } else {
      for (J = 1; J <= N; J++) {
        // 140
        if (UDIAG) {
          SUM.value = ONE;
          for (I = K + 1; I <= K + N - J; I++) {
            // 120
            SUM.value = SUM.value + AP[I].abs();
          } // 120
        } else {
          SUM.value = ZERO;
          for (I = K; I <= K + N - J; I++) {
            // 130
            SUM.value = SUM.value + AP[I].abs();
          } // 130
        }
        K += N - J + 1;
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 140
    }
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    K = 1;
    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= N; I++) {
          // 150
          WORK[I] = ONE;
        } // 150
        for (J = 1; J <= N; J++) {
          // 170
          for (I = 1; I <= J - 1; I++) {
            // 160
            WORK[I] += AP[K].abs();
            K++;
          } // 160
          K++;
        } // 170
      } else {
        for (I = 1; I <= N; I++) {
          // 180
          WORK[I] = ZERO;
        } // 180
        for (J = 1; J <= N; J++) {
          // 200
          for (I = 1; I <= J; I++) {
            // 190
            WORK[I] += AP[K].abs();
            K++;
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
          K++;
          for (I = J + 1; I <= N; I++) {
            // 220
            WORK[I] += AP[K].abs();
            K++;
          } // 220
        } // 230
      } else {
        for (I = 1; I <= N; I++) {
          // 240
          WORK[I] = ZERO;
        } // 240
        for (J = 1; J <= N; J++) {
          // 260
          for (I = J; I <= N; I++) {
            // 250
            WORK[I] += AP[K].abs();
            K++;
          } // 250
        } // 260
      }
    }
    VALUE = ZERO;
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
        K = 2;
        for (J = 2; J <= N; J++) {
          // 280
          zlassq(J - 1, AP(K), 1, SCALE, SUM);
          K += J;
        } // 280
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        K = 1;
        for (J = 1; J <= N; J++) {
          // 290
          zlassq(J, AP(K), 1, SCALE, SUM);
          K += J;
        } // 290
      }
    } else {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = N.toDouble();
        K = 2;
        for (J = 1; J <= N - 1; J++) {
          // 300
          zlassq(N - J, AP(K), 1, SCALE, SUM);
          K += N - J + 1;
        } // 300
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        K = 1;
        for (J = 1; J <= N; J++) {
          // 310
          zlassq(N - J + 1, AP(K), 1, SCALE, SUM);
          K += N - J + 1;
        } // 310
      }
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
