// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/dlassq.dart';
import 'package:dart_lapack/src/matrix.dart';

double dlantp(
  final String NORM,
  final String UPLO,
  final String DIAG,
  final int N,
  final Array<double> AP_,
  final Array<double> WORK_,
) {
  final AP = AP_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UDIAG;
  int I, J, K;
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
          for (I = K; I <= K + J - 2; I++) {
            SUM.value = AP[I].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
          K += J;
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = K + 1; I <= K + N - J; I++) {
            SUM.value = AP[I].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
          K += N - J + 1;
        }
      }
    } else {
      VALUE = ZERO;
      if (lsame(UPLO, 'U')) {
        for (J = 1; J <= N; J++) {
          for (I = K; I <= K + J - 1; I++) {
            SUM.value = AP[I].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
          K += J;
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = K; I <= K + N - J; I++) {
            SUM.value = AP[I].abs();
            if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
          }
          K += N - J + 1;
        }
      }
    }
  } else if ((lsame(NORM, 'O')) || (NORM == '1')) {
    // Find norm1(A).

    VALUE = ZERO;
    K = 1;
    UDIAG = lsame(DIAG, 'U');
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        if (UDIAG) {
          SUM.value = ONE;
          for (I = K; I <= K + J - 2; I++) {
            SUM.value += AP[I].abs();
          }
        } else {
          SUM.value = ZERO;
          for (I = K; I <= K + J - 1; I++) {
            SUM.value += AP[I].abs();
          }
        }
        K += J;
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    } else {
      for (J = 1; J <= N; J++) {
        if (UDIAG) {
          SUM.value = ONE;
          for (I = K + 1; I <= K + N - J; I++) {
            SUM.value += AP[I].abs();
          }
        } else {
          SUM.value = ZERO;
          for (I = K; I <= K + N - J; I++) {
            SUM.value += AP[I].abs();
          }
        }
        K += N - J + 1;
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    }
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    K = 1;
    if (lsame(UPLO, 'U')) {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= N; I++) {
          WORK[I] = ONE;
        }
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J - 1; I++) {
            WORK[I] += AP[K].abs();
            K++;
          }
          K++;
        }
      } else {
        for (I = 1; I <= N; I++) {
          WORK[I] = ZERO;
        }
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            WORK[I] += AP[K].abs();
            K++;
          }
        }
      }
    } else {
      if (lsame(DIAG, 'U')) {
        for (I = 1; I <= N; I++) {
          WORK[I] = ONE;
        }
        for (J = 1; J <= N; J++) {
          K++;
          for (I = J + 1; I <= N; I++) {
            WORK[I] += AP[K].abs();
            K++;
          }
        }
      } else {
        for (I = 1; I <= N; I++) {
          WORK[I] = ZERO;
        }
        for (J = 1; J <= N; J++) {
          for (I = J; I <= N; I++) {
            WORK[I] += AP[K].abs();
            K++;
          }
        }
      }
    }
    VALUE = ZERO;
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
        K = 2;
        for (J = 2; J <= N; J++) {
          dlassq(J - 1, AP(K), 1, SCALE, SUM);
          K += J;
        }
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        K = 1;
        for (J = 1; J <= N; J++) {
          dlassq(J, AP(K), 1, SCALE, SUM);
          K += J;
        }
      }
    } else {
      if (lsame(DIAG, 'U')) {
        SCALE.value = ONE;
        SUM.value = N.toDouble();
        K = 2;
        for (J = 1; J <= N - 1; J++) {
          dlassq(N - J, AP(K), 1, SCALE, SUM);
          K += N - J + 1;
        }
      } else {
        SCALE.value = ZERO;
        SUM.value = ONE;
        K = 1;
        for (J = 1; J <= N; J++) {
          dlassq(N - J + 1, AP(K), 1, SCALE, SUM);
          K += N - J + 1;
        }
      }
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
