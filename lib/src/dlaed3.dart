// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlaed4.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlaed3(
  final int K,
  final int N,
  final int N1,
  final Array<double> D_,
  final Matrix<double> Q_,
  final int LDQ,
  final double RHO,
  final Array<double> DLAMBDA_,
  final Array<double> Q2_,
  final Array<int> INDX_,
  final Array<int> CTOT_,
  final Array<double> W_,
  final Array<double> S_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final Q = Q_.having(ld: LDQ);
  final DLAMBDA = DLAMBDA_.having();
  final Q2 = Q2_.having();
  final INDX = INDX_.having();
  final CTOT = CTOT_.having();
  final W = W_.having();
  final S = S_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, II, IQ2, J, N12, N2, N23;
  double TEMP;

  // Test the input parameters.

  INFO.value = 0;

  if (K < 0) {
    INFO.value = -1;
  } else if (N < K) {
    INFO.value = -2;
  } else if (LDQ < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DLAED3', -INFO.value);
    return;
  }

  // Quick return if possible

  if (K == 0) return;

  for (J = 1; J <= K; J++) {
    dlaed4(K, J, DLAMBDA, W, Q(1, J).asArray(), RHO, D.box(J), INFO);

    // If the zero finder fails, the computation is terminated.

    if (INFO.value != 0) return;
  }

  if (K == 2) {
    for (J = 1; J <= K; J++) {
      W[1] = Q[1][J];
      W[2] = Q[2][J];
      II = INDX[1];
      Q[1][J] = W[II];
      II = INDX[2];
      Q[2][J] = W[II];
    }
  } else if (K != 1) {
    // Compute updated W.

    dcopy(K, W, 1, S, 1);

    // Initialize W[I] = Q[I][I]

    dcopy(K, Q.asArray(), LDQ + 1, W, 1);
    for (J = 1; J <= K; J++) {
      for (I = 1; I <= J - 1; I++) {
        W[I] *= (Q[I][J] / (DLAMBDA[I] - DLAMBDA[J]));
      }
      for (I = J + 1; I <= K; I++) {
        W[I] *= (Q[I][J] / (DLAMBDA[I] - DLAMBDA[J]));
      }
    }
    for (I = 1; I <= K; I++) {
      W[I] = sign(sqrt(-W[I]), S[I]);
    }

    // Compute eigenvectors of the modified rank-1 modification.

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= K; I++) {
        S[I] = W[I] / Q[I][J];
      }
      TEMP = dnrm2(K, S, 1);
      for (I = 1; I <= K; I++) {
        II = INDX[I];
        Q[I][J] = S[II] / TEMP;
      }
    }

    // Compute the updated eigenvectors.
  }

  N2 = N - N1;
  N12 = CTOT[1] + CTOT[2];
  N23 = CTOT[2] + CTOT[3];

  dlacpy('A', N23, K, Q(CTOT[1] + 1, 1), LDQ, S.asMatrix(N23), N23);
  IQ2 = N1 * N12 + 1;
  if (N23 != 0) {
    dgemm('N', 'N', N2, K, N23, ONE, Q2(IQ2).asMatrix(N2), N2, S.asMatrix(N23),
        N23, ZERO, Q(N1 + 1, 1), LDQ);
  } else {
    dlaset('A', N2, K, ZERO, ZERO, Q(N1 + 1, 1), LDQ);
  }

  dlacpy('A', N12, K, Q, LDQ, S.asMatrix(N12), N12);
  if (N12 != 0) {
    dgemm('N', 'N', N1, K, N12, ONE, Q2.asMatrix(N1), N1, S.asMatrix(N12), N12,
        ZERO, Q, LDQ);
  } else {
    dlaset('A', N1, K, ZERO, ZERO, Q(1, 1), LDQ);
  }
}
