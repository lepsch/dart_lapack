// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlassq.dart';

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
          for (I = max(K + 2 - J, 1); I <= K; I++) {
            SUM.value = AB[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = 2; I <= min(N + 1 - J, K + 1); I++) {
            SUM.value = AB[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
        }
      }
    } else {
      VALUE = ZERO;
      if (lsame(UPLO, 'U')) {
        for (J = 1; J <= N; J++) {
          for (I = max(K + 2 - J, 1); I <= K + 1; I++) {
            SUM.value = AB[I][J].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= min(N + 1 - J, K + 1); I++) {
            SUM.value = AB[I][J].abs();
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
        if (UDIAG) {
          SUM.value = ONE;
          for (I = max(K + 2 - J, 1); I <= K; I++) {
            SUM.value += AB[I][J].abs();
          }
        } else {
          SUM.value = ZERO;
          for (I = max(K + 2 - J, 1); I <= K + 1; I++) {
            SUM.value += AB[I][J].abs();
          }
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    } else {
      for (J = 1; J <= N; J++) {
        if (UDIAG) {
          SUM.value = ONE;
          for (I = 2; I <= min(N + 1 - J, K + 1); I++) {
            SUM.value += AB[I][J].abs();
          }
        } else {
          SUM.value = ZERO;
          for (I = 1; I <= min(N + 1 - J, K + 1); I++) {
            SUM.value += AB[I][J].abs();
          }
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    }
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    VALUE = ZERO;
    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= N; I++) {
          WORK[I] = ONE;
        }
        for (J = 1; J <= N; J++) {
          L = K + 1 - J;
          for (I = max(1, J - K); I <= J - 1; I++) {
            WORK[I] += AB[L + I][J].abs();
          }
        }
      } else {
        for (I = 1; I <= N; I++) {
          WORK[I] = ZERO;
        }
        for (J = 1; J <= N; J++) {
          L = K + 1 - J;
          for (I = max(1, J - K); I <= J; I++) {
            WORK[I] += AB[L + I][J].abs();
          }
        }
      }
    } else {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= N; I++) {
          WORK[I] = ONE;
        }
        for (J = 1; J <= N; J++) {
          L = 1 - J;
          for (I = J + 1; I <= min(N, J + K); I++) {
            WORK[I] += AB[L + I][J].abs();
          }
        }
      } else {
        for (I = 1; I <= N; I++) {
          WORK[I] = ZERO;
        }
        for (J = 1; J <= N; J++) {
          L = 1 - J;
          for (I = J; I <= min(N, J + K); I++) {
            WORK[I] += AB[L + I][J].abs();
          }
        }
      }
    }
    for (I = 1; I <= N; I++) {
      SUM.value = WORK[I];
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    }
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = N.toDouble();
        if (K > 0) {
          for (J = 2; J <= N; J++) {
            zlassq(min(J - 1, K), AB(max(K + 2 - J, 1), J).asArray(), 1, SCALE,
                SUM);
          }
        }
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        for (J = 1; J <= N; J++) {
          zlassq(
              min(J, K + 1), AB(max(K + 2 - J, 1), J).asArray(), 1, SCALE, SUM);
        }
      }
    } else {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = N.toDouble();
        if (K > 0) {
          for (J = 1; J <= N - 1; J++) {
            zlassq(min(N - J, K), AB(2, J).asArray(), 1, SCALE, SUM);
          }
        }
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        for (J = 1; J <= N; J++) {
          zlassq(min(N - J + 1, K + 1), AB(1, J).asArray(), 1, SCALE, SUM);
        }
      }
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
