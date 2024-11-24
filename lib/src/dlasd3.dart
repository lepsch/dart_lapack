// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/dlasd4.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlasd3(
  final int NL,
  final int NR,
  final int SQRE,
  final int K,
  final Array<double> D_,
  final Matrix<double> Q_,
  final int LDQ,
  final Array<double> DSIGMA_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> U2_,
  final int LDU2,
  final Matrix<double> VT_,
  final int LDVT,
  final Matrix<double> VT2_,
  final int LDVT2,
  final Array<int> IDXC_,
  final Array<int> CTOT_,
  final Array<double> Z_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final Q = Q_.having(ld: LDQ);
  final DSIGMA = DSIGMA_.having();
  final U = U_.having(ld: LDU);
  final U2 = U2_.having(ld: LDU2);
  final VT = VT_.having(ld: LDVT);
  final VT2 = VT2_.having(ld: LDVT2);
  final IDXC = IDXC_.having();
  final CTOT = CTOT_.having();
  final Z = Z_.having();
  const ONE = 1.0, ZERO = 0.0, NEGONE = -1.0;
  int CTEMP, I, J, JC, KTEMP, M, N, NLP1, NLP2, NRP1;
  double RHO, TEMP;

  // Test the input parameters.

  INFO.value = 0;

  if (NL < 1) {
    INFO.value = -1;
  } else if (NR < 1) {
    INFO.value = -2;
  } else if ((SQRE != 1) && (SQRE != 0)) {
    INFO.value = -3;
  }

  N = NL + NR + 1;
  M = N + SQRE;
  NLP1 = NL + 1;
  NLP2 = NL + 2;

  if ((K < 1) || (K > N)) {
    INFO.value = -4;
  } else if (LDQ < K) {
    INFO.value = -7;
  } else if (LDU < N) {
    INFO.value = -10;
  } else if (LDU2 < N) {
    INFO.value = -12;
  } else if (LDVT < M) {
    INFO.value = -14;
  } else if (LDVT2 < M) {
    INFO.value = -16;
  }
  if (INFO.value != 0) {
    xerbla('DLASD3', -INFO.value);
    return;
  }

  // Quick return if possible

  if (K == 1) {
    D[1] = Z[1].abs();
    dcopy(M, VT2(1, 1).asArray(), LDVT2, VT(1, 1).asArray(), LDVT);
    if (Z[1] > ZERO) {
      dcopy(N, U2(1, 1).asArray(), 1, U(1, 1).asArray(), 1);
    } else {
      for (I = 1; I <= N; I++) {
        U[I][1] = -U2[I][1];
      }
    }
    return;
  }

  // Keep a copy of Z.

  dcopy(K, Z, 1, Q.asArray(), 1);

  // Normalize Z.

  RHO = dnrm2(K, Z, 1);
  dlascl('G', 0, 0, RHO, ONE, K, 1, Z.asMatrix(K), K, INFO);
  RHO *= RHO;

  // Find the new singular values.

  for (J = 1; J <= K; J++) {
    dlasd4(K, J, DSIGMA, Z, U(1, J).asArray(), RHO, D.box(J),
        VT(1, J).asArray(), INFO);

    // If the zero finder fails, report the convergence failure.

    if (INFO.value != 0) {
      return;
    }
  }

  // Compute updated Z.

  for (I = 1; I <= K; I++) {
    Z[I] = U[I][K] * VT[I][K];
    for (J = 1; J <= I - 1; J++) {
      Z[I] *= (U[I][J] *
          VT[I][J] /
          (DSIGMA[I] - DSIGMA[J]) /
          (DSIGMA[I] + DSIGMA[J]));
    }
    for (J = I; J <= K - 1; J++) {
      Z[I] *= (U[I][J] *
          VT[I][J] /
          (DSIGMA[I] - DSIGMA[J + 1]) /
          (DSIGMA[I] + DSIGMA[J + 1]));
    }
    Z[I] = sign(sqrt(Z[I].abs()), Q[I][1]);
  }

  // Compute left singular vectors of the modified diagonal matrix,
  // and store related information for the right singular vectors.

  for (I = 1; I <= K; I++) {
    VT[1][I] = Z[1] / U[1][I] / VT[1][I];
    U[1][I] = NEGONE;
    for (J = 2; J <= K; J++) {
      VT[J][I] = Z[J] / U[J][I] / VT[J][I];
      U[J][I] = DSIGMA[J] * VT[J][I];
    }
    TEMP = dnrm2(K, U(1, I).asArray(), 1);
    Q[1][I] = U[1][I] / TEMP;
    for (J = 2; J <= K; J++) {
      JC = IDXC[J];
      Q[J][I] = U[JC][I] / TEMP;
    }
  }

  // Update the left singular vector matrix.

  if (K == 2) {
    dgemm('N', 'N', N, K, K, ONE, U2, LDU2, Q, LDQ, ZERO, U, LDU);
  } else {
    if (CTOT[1] > 0) {
      dgemm('N', 'N', NL, K, CTOT[1], ONE, U2(1, 2), LDU2, Q(2, 1), LDQ, ZERO,
          U(1, 1), LDU);
      if (CTOT[3] > 0) {
        KTEMP = 2 + CTOT[1] + CTOT[2];
        dgemm('N', 'N', NL, K, CTOT[3], ONE, U2(1, KTEMP), LDU2, Q(KTEMP, 1),
            LDQ, ONE, U(1, 1), LDU);
      }
    } else if (CTOT[3] > 0) {
      KTEMP = 2 + CTOT[1] + CTOT[2];
      dgemm('N', 'N', NL, K, CTOT[3], ONE, U2(1, KTEMP), LDU2, Q(KTEMP, 1), LDQ,
          ZERO, U(1, 1), LDU);
    } else {
      dlacpy('F', NL, K, U2, LDU2, U, LDU);
    }
    dcopy(K, Q(1, 1).asArray(), LDQ, U(NLP1, 1).asArray(), LDU);
    KTEMP = 2 + CTOT[1];
    CTEMP = CTOT[2] + CTOT[3];
    dgemm('N', 'N', NR, K, CTEMP, ONE, U2(NLP2, KTEMP), LDU2, Q(KTEMP, 1), LDQ,
        ZERO, U(NLP2, 1), LDU);
  }

  // Generate the right singular vectors.

  for (I = 1; I <= K; I++) {
    TEMP = dnrm2(K, VT(1, I).asArray(), 1);
    Q[I][1] = VT[1][I] / TEMP;
    for (J = 2; J <= K; J++) {
      JC = IDXC[J];
      Q[I][J] = VT[JC][I] / TEMP;
    }
  }

  // Update the right singular vector matrix.

  if (K == 2) {
    dgemm('N', 'N', K, M, K, ONE, Q, LDQ, VT2, LDVT2, ZERO, VT, LDVT);
    return;
  }
  KTEMP = 1 + CTOT[1];
  dgemm('N', 'N', K, NLP1, KTEMP, ONE, Q(1, 1), LDQ, VT2(1, 1), LDVT2, ZERO,
      VT(1, 1), LDVT);
  KTEMP = 2 + CTOT[1] + CTOT[2];
  if (KTEMP <= LDVT2) {
    dgemm('N', 'N', K, NLP1, CTOT[3], ONE, Q(1, KTEMP), LDQ, VT2(KTEMP, 1),
        LDVT2, ONE, VT(1, 1), LDVT);
  }

  KTEMP = CTOT[1] + 1;
  NRP1 = NR + SQRE;
  if (KTEMP > 1) {
    for (I = 1; I <= K; I++) {
      Q[I][KTEMP] = Q[I][1];
    }
    for (I = NLP2; I <= M; I++) {
      VT2[KTEMP][I] = VT2[1][I];
    }
  }
  CTEMP = 1 + CTOT[2] + CTOT[3];
  dgemm('N', 'N', K, NRP1, CTEMP, ONE, Q(1, KTEMP), LDQ, VT2(KTEMP, NLP2),
      LDVT2, ZERO, VT(1, NLP2), LDVT);
}
