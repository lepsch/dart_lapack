// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacn2.dart';
import 'package:dart_lapack/src/zspmv.dart';
import 'package:dart_lapack/src/zsptrs.dart';

void zsprfs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<Complex> AP_,
  final Array<Complex> AFP_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> X_,
  final int LDX,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final AP = AP_.having();
  final AFP = AFP_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  const ITMAX = 5;
  const ZERO = 0.0;
  const TWO = 2.0;
  const THREE = 3.0;
  bool UPPER;
  int COUNT, I, IK, J, K, KK, NZ;
  double EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
  final ISAVE = Array<int>(3);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  } else if (LDX < max(1, N)) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('ZSPRFS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) {
    for (J = 1; J <= NRHS; J++) {
      FERR[J] = ZERO;
      BERR[J] = ZERO;
    }
    return;
  }

  // NZ = maximum number of nonzero elements in each row of A, plus 1

  NZ = N + 1;
  EPS = dlamch('Epsilon');
  SAFMIN = dlamch('Safe minimum');
  SAFE1 = NZ * SAFMIN;
  SAFE2 = SAFE1 / EPS;

  // Do for each right hand side

  for (J = 1; J <= NRHS; J++) {
    COUNT = 1;
    LSTRES = THREE;
    while (true) {
      // Loop until stopping criterion is satisfied.

      // Compute residual R = B - A * X

      zcopy(N, B(1, J).asArray(), 1, WORK, 1);
      zspmv(UPLO, N, -Complex.one, AP, X(1, J).asArray(), 1, Complex.one, WORK,
          1);

      // Compute componentwise relative backward error from formula

      // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

      // where abs(Z) is the componentwise absolute value of the matrix
      // or vector Z.  If the i-th component of the denominator is less
      // than SAFE2, then SAFE1 is added to the i-th components of the
      // numerator and denominator before dividing.

      for (I = 1; I <= N; I++) {
        RWORK[I] = B[I][J].cabs1();
      }

      // Compute abs(A)*abs(X) + abs(B).

      KK = 1;
      if (UPPER) {
        for (K = 1; K <= N; K++) {
          S = ZERO;
          XK = X[K][J].cabs1();
          IK = KK;
          for (I = 1; I <= K - 1; I++) {
            RWORK[I] += AP[IK].cabs1() * XK;
            S += AP[IK].cabs1() * X[I][J].cabs1();
            IK++;
          }
          RWORK[K] += AP[KK + K - 1].cabs1() * XK + S;
          KK += K;
        }
      } else {
        for (K = 1; K <= N; K++) {
          S = ZERO;
          XK = X[K][J].cabs1();
          RWORK[K] += AP[KK].cabs1() * XK;
          IK = KK + 1;
          for (I = K + 1; I <= N; I++) {
            RWORK[I] += AP[IK].cabs1() * XK;
            S += AP[IK].cabs1() * X[I][J].cabs1();
            IK++;
          }
          RWORK[K] += S;
          KK += (N - K + 1);
        }
      }
      S = ZERO;
      for (I = 1; I <= N; I++) {
        if (RWORK[I] > SAFE2) {
          S = max(S, WORK[I].cabs1() / RWORK[I]);
        } else {
          S = max(S, (WORK[I].cabs1() + SAFE1) / (RWORK[I] + SAFE1));
        }
      }
      BERR[J] = S;

      // Test stopping criterion. Continue iterating if
      //    1) The residual BERR(J) is larger than machine epsilon, and
      //    2) BERR(J) decreased by at least a factor of 2 during the
      //       last iteration, and
      //    3) At most ITMAX iterations tried.

      if (BERR[J] > EPS && TWO * BERR[J] <= LSTRES && COUNT <= ITMAX) {
        // Update solution and try again.

        zsptrs(UPLO, N, 1, AFP, IPIV, WORK.asMatrix(), N, INFO);
        zaxpy(N, Complex.one, WORK, 1, X(1, J).asArray(), 1);
        LSTRES = BERR[J];
        COUNT++;
        continue;
      }
      break;
    }

    // Bound error from formula

    // norm(X - XTRUE) / norm(X) <= FERR =
    // norm( abs(inv(A))*
    //    ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)

    // where
    //   norm(Z) is the magnitude of the largest component of Z
    //   inv(A) is the inverse of A
    //   abs(Z) is the componentwise absolute value of the matrix or
    //      vector Z
    //   NZ is the maximum number of nonzeros in any row of A, plus 1
    //   EPS is machine epsilon

    // The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
    // is incremented by SAFE1 if the i-th component of
    // abs(A)*abs(X) + abs(B) is less than SAFE2.

    // Use ZLACN2 to estimate the infinity-norm of the matrix
    //    inv(A) * diag(W),
    // where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))

    for (I = 1; I <= N; I++) {
      if (RWORK[I] > SAFE2) {
        RWORK[I] = WORK[I].cabs1() + NZ * EPS * RWORK[I];
      } else {
        RWORK[I] = WORK[I].cabs1() + NZ * EPS * RWORK[I] + SAFE1;
      }
    }

    KASE.value = 0;
    while (true) {
      zlacn2(N, WORK(N + 1), WORK, FERR(J), KASE, ISAVE);
      if (KASE.value == 0) break;
      if (KASE.value == 1) {
        // Multiply by diag(W)*inv(A**T).

        zsptrs(UPLO, N, 1, AFP, IPIV, WORK.asMatrix(), N, INFO);
        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
      } else if (KASE.value == 2) {
        // Multiply by inv(A)*diag(W).

        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
        zsptrs(UPLO, N, 1, AFP, IPIV, WORK.asMatrix(), N, INFO);
      }
    }

    // Normalize error.

    LSTRES = ZERO;
    for (I = 1; I <= N; I++) {
      LSTRES = max(LSTRES, X[I][J].cabs1());
    }
    if (LSTRES != ZERO) FERR[J] /= LSTRES;
  }
}
