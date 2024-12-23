// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zdrot.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dlamrg.dart';
import 'package:dart_lapack/src/dlapy2.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacpy.dart';

void zlaed8(
  final Box<int> K,
  final int N,
  final int QSIZ,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<double> D_,
  final Box<double> RHO,
  final int CUTPNT,
  final Array<double> Z_,
  final Array<double> DLAMBDA_,
  final Matrix<Complex> Q2_,
  final int LDQ2,
  final Array<double> W_,
  final Array<int> INDXP_,
  final Array<int> INDX_,
  final Array<int> INDXQ_,
  final Array<int> PERM_,
  final Box<int> GIVPTR,
  final Matrix<int> GIVCOL_,
  final Matrix<double> GIVNUM_,
  final Box<int> INFO,
) {
  final Q = Q_.having(ld: LDQ);
  final Q2 = Q2_.having(ld: LDQ2);
  final D = D_.having();
  final Z = Z_.having();
  final DLAMBDA = DLAMBDA_.having();
  final W = W_.having();
  final GIVNUM = GIVNUM_.having(ld: 2);
  final INDXP = INDXP_.having();
  final INDX = INDX_.having();
  final INDXQ = INDXQ_.having();
  final PERM = PERM_.having();
  final GIVCOL = GIVCOL_.having(ld: 2);
  const MONE = -1.0, ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0;
  int I, IMAX, J, JLAM = 0, JMAX, JP, K2, N1, N1P1, N2;
  double C, EPS, S, T, TAU, TOL;

  // Test the input parameters.

  INFO.value = 0;

  if (N < 0) {
    INFO.value = -2;
  } else if (QSIZ < N) {
    INFO.value = -3;
  } else if (LDQ < max(1, N)) {
    INFO.value = -5;
  } else if (CUTPNT < min(1, N) || CUTPNT > N) {
    INFO.value = -8;
  } else if (LDQ2 < max(1, N)) {
    INFO.value = -12;
  }
  if (INFO.value != 0) {
    xerbla('ZLAED8', -INFO.value);
    return;
  }

  // Need to initialize GIVPTR to O here in case of quick exit
  // to prevent an unspecified code behavior (usually sigfault)
  // when IWORK array on entry to *stedc is not zeroed
  // (or at least some IWORK entries which used in *laed7 for GIVPTR).

  GIVPTR.value = 0;

  // Quick return if possible

  if (N == 0) return;

  N1 = CUTPNT;
  N2 = N - N1;
  N1P1 = N1 + 1;

  if (RHO.value < ZERO) {
    dscal(N2, MONE, Z(N1P1), 1);
  }

  // Normalize z so that norm(z) = 1

  T = ONE / sqrt(TWO);
  for (J = 1; J <= N; J++) {
    INDX[J] = J;
  }
  dscal(N, T, Z, 1);
  RHO.value = (TWO * RHO.value).abs();

  // Sort the eigenvalues into increasing order

  for (I = CUTPNT + 1; I <= N; I++) {
    INDXQ[I] += CUTPNT;
  }
  for (I = 1; I <= N; I++) {
    DLAMBDA[I] = D[INDXQ[I]];
    W[I] = Z[INDXQ[I]];
  }
  I = 1;
  J = CUTPNT + 1;
  dlamrg(N1, N2, DLAMBDA, 1, 1, INDX);
  for (I = 1; I <= N; I++) {
    D[I] = DLAMBDA[INDX[I]];
    Z[I] = W[INDX[I]];
  }

  // Calculate the allowable deflation tolerance

  IMAX = idamax(N, Z, 1);
  JMAX = idamax(N, D, 1);
  EPS = dlamch('Epsilon');
  TOL = EIGHT * EPS * D[JMAX].abs();

  // If the rank-1 modifier is small enough, no more needs to be done
  // -- except to reorganize Q so that its columns correspond with the
  // elements in D.

  if (RHO.value * Z[IMAX].abs() <= TOL) {
    K.value = 0;
    for (J = 1; J <= N; J++) {
      PERM[J] = INDXQ[INDX[J]];
      zcopy(QSIZ, Q(1, PERM[J]).asArray(), 1, Q2(1, J).asArray(), 1);
    }
    zlacpy('A', QSIZ, N, Q2(1, 1), LDQ2, Q(1, 1), LDQ);
    return;
  }

  // If there are multiple eigenvalues then the problem deflates.  Here
  // the number of equal eigenvalues are found.  As each equal
  // eigenvalue is found, an elementary reflector is computed to rotate
  // the corresponding eigensubspace so that the corresponding
  // components of Z are zero in this new basis.

  K.value = 0;
  K2 = N + 1;
  var deflate = false;
  for (J = 1; J <= N; J++) {
    if (RHO.value * Z[J].abs() <= TOL) {
      // Deflate due to small z component.

      K2--;
      INDXP[K2] = J;
      if (J == N) {
        deflate = true;
        break;
      }
    } else {
      JLAM = J;
      break;
    }
  }

  while (!deflate) {
    J++;
    if (J > N) break;
    if (RHO.value * Z[J].abs() <= TOL) {
      // Deflate due to small z component.

      K2--;
      INDXP[K2] = J;
    } else {
      // Check if eigenvalues are close enough to allow deflation.

      S = Z[JLAM];
      C = Z[J];

      // Find sqrt(a**2+b**2) without overflow or
      // destructive underflow.

      TAU = dlapy2(C, S);
      T = D[J] - D[JLAM];
      C /= TAU;
      S = -S / TAU;
      if ((T * C * S).abs() <= TOL) {
        // Deflation is possible.

        Z[J] = TAU;
        Z[JLAM] = ZERO;

        // Record the appropriate Givens rotation

        GIVPTR.value++;
        GIVCOL[1][GIVPTR.value] = INDXQ[INDX[JLAM]];
        GIVCOL[2][GIVPTR.value] = INDXQ[INDX[J]];
        GIVNUM[1][GIVPTR.value] = C;
        GIVNUM[2][GIVPTR.value] = S;
        zdrot(QSIZ, Q(1, INDXQ[INDX[JLAM]]).asArray(), 1,
            Q(1, INDXQ[INDX[J]]).asArray(), 1, C, S);
        T = D[JLAM] * C * C + D[J] * S * S;
        D[J] = D[JLAM] * S * S + D[J] * C * C;
        D[JLAM] = T;
        K2--;
        I = 1;
        while (true) {
          if (K2 + I <= N) {
            if (D[JLAM] < D[INDXP[K2 + I]]) {
              INDXP[K2 + I - 1] = INDXP[K2 + I];
              INDXP[K2 + I] = JLAM;
              I++;
              continue;
            } else {
              INDXP[K2 + I - 1] = JLAM;
            }
          } else {
            INDXP[K2 + I - 1] = JLAM;
          }
          break;
        }
        JLAM = J;
      } else {
        K.value++;
        W[K.value] = Z[JLAM];
        DLAMBDA[K.value] = D[JLAM];
        INDXP[K.value] = JLAM;
        JLAM = J;
      }
    }
  }

  if (!deflate) {
    // Record the last eigenvalue.

    K.value++;
    W[K.value] = Z[JLAM];
    DLAMBDA[K.value] = D[JLAM];
    INDXP[K.value] = JLAM;
  }

  // Sort the eigenvalues and corresponding eigenvectors into DLAMBDA
  // and Q2 respectively.  The eigenvalues/vectors which were not
  // deflated go into the first K slots of DLAMBDA and Q2 respectively,
  // while those which were deflated go into the last N - K slots.

  for (J = 1; J <= N; J++) {
    JP = INDXP[J];
    DLAMBDA[J] = D[JP];
    PERM[J] = INDXQ[INDX[JP]];
    zcopy(QSIZ, Q(1, PERM[J]).asArray(), 1, Q2(1, J).asArray(), 1);
  }

  // The deflated eigenvalues and their corresponding vectors go back
  // into the last N - K slots of D and Q respectively.

  if (K.value < N) {
    dcopy(N - K.value, DLAMBDA(K.value + 1), 1, D(K.value + 1), 1);
    zlacpy('A', QSIZ, N - K.value, Q2(1, K.value + 1), LDQ2, Q(1, K.value + 1),
        LDQ);
  }
}
