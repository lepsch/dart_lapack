// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zhbmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacn2.dart';
import 'package:lapack/src/zpbtrs.dart';

void zpbrfs(
  final String UPLO,
  final int N,
  final int KD,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> AFB_,
  final int LDAFB,
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
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final AFB = AFB_.having(ld: LDAFB);
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final RWORK = RWORK_.having();

  const ITMAX = 5;
  const ZERO = 0.0;
  const TWO = 2.0;
  const THREE = 3.0;
  bool UPPER;
  int COUNT = 0, I, J, K, L, NZ;
  double EPS, LSTRES = 0, S, SAFE1, SAFE2, SAFMIN, XK;
  final ISAVE = Array<int>(3);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KD < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDAB < KD + 1) {
    INFO.value = -6;
  } else if (LDAFB < KD + 1) {
    INFO.value = -8;
  } else if (LDB < max(1, N)) {
    INFO.value = -10;
  } else if (LDX < max(1, N)) {
    INFO.value = -12;
  }
  if (INFO.value != 0) {
    xerbla('ZPBRFS', -INFO.value);
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

  NZ = min(N + 1, 2 * KD + 2);
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
      zhbmv(UPLO, N, KD, -Complex.one, AB, LDAB, X(1, J).asArray(), 1,
          Complex.one, WORK, 1);

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

      if (UPPER) {
        for (K = 1; K <= N; K++) {
          S = ZERO;
          XK = X[K][J].cabs1();
          L = KD + 1 - K;
          for (I = max(1, K - KD); I <= K - 1; I++) {
            RWORK[I] += AB[L + I][K].cabs1() * XK;
            S += AB[L + I][K].cabs1() * X[I][J].cabs1();
          }
          RWORK[K] += AB[KD + 1][K].real.abs() * XK + S;
        }
      } else {
        for (K = 1; K <= N; K++) {
          S = ZERO;
          XK = X[K][J].cabs1();
          RWORK[K] += AB[1][K].real.abs() * XK;
          L = 1 - K;
          for (I = K + 1; I <= min(N, K + KD); I++) {
            RWORK[I] += AB[L + I][K].cabs1() * XK;
            S += AB[L + I][K].cabs1() * X[I][J].cabs1();
          }
          RWORK[K] += S;
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

        zpbtrs(UPLO, N, KD, 1, AFB, LDAFB, WORK.asMatrix(), N, INFO);
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
        // Multiply by diag(W)*inv(A**H).

        zpbtrs(UPLO, N, KD, 1, AFB, LDAFB, WORK.asMatrix(), N, INFO);
        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
      } else if (KASE.value == 2) {
        // Multiply by inv(A)*diag(W).

        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
        zpbtrs(UPLO, N, KD, 1, AFB, LDAFB, WORK.asMatrix(), N, INFO);
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
