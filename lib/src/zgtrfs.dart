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
import 'package:dart_lapack/src/zgttrs.dart';
import 'package:dart_lapack/src/zlacn2.dart';
import 'package:dart_lapack/src/zlagtm.dart';

void zgtrfs(
  final String TRANS,
  final int N,
  final int NRHS,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Array<Complex> DLF_,
  final Array<Complex> DF_,
  final Array<Complex> DUF_,
  final Array<Complex> DU2_,
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
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DLF = DLF_.having();
  final DF = DF_.having();
  final DUF = DUF_.having();
  final DU2 = DU2_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  const ITMAX = 5;
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0;
  const THREE = 3.0;
  bool NOTRAN;
  String TRANSN, TRANST;
  int COUNT = 0, I, J, NZ;
  double EPS, LSTRES = 0, S, SAFE1, SAFE2, SAFMIN;
  final ISAVE = Array<int>(3);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  NOTRAN = lsame(TRANS, 'N');
  if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -13;
  } else if (LDX < max(1, N)) {
    INFO.value = -15;
  }
  if (INFO.value != 0) {
    xerbla('ZGTRFS', -INFO.value);
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

  if (NOTRAN) {
    TRANSN = 'N';
    TRANST = 'C';
  } else {
    TRANSN = 'C';
    TRANST = 'N';
  }

  // NZ = maximum number of nonzero elements in each row of A, plus 1

  NZ = 4;
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

      // Compute residual R = B - op(A) * X,
      // where op(A) = A, A**T, or A**H, depending on TRANS.

      zcopy(N, B(1, J).asArray(), 1, WORK, 1);
      zlagtm(
          TRANS, N, 1, -ONE, DL, D, DU, X(1, J), LDX, ONE, WORK.asMatrix(N), N);

      // Compute abs(op(A))*abs(x) + abs(b) for use in the backward
      // error bound.

      if (NOTRAN) {
        if (N == 1) {
          RWORK[1] = B[1][J].cabs1() + D[1].cabs1() * X[1][J].cabs1();
        } else {
          RWORK[1] = B[1][J].cabs1() +
              D[1].cabs1() * X[1][J].cabs1() +
              DU[1].cabs1() * X[2][J].cabs1();
          for (I = 2; I <= N - 1; I++) {
            RWORK[I] = B[I][J].cabs1() +
                DL[I - 1].cabs1() * X[I - 1][J].cabs1() +
                D[I].cabs1() * X[I][J].cabs1() +
                DU[I].cabs1() * X[I + 1][J].cabs1();
          }
          RWORK[N] = B[N][J].cabs1() +
              DL[N - 1].cabs1() * X[N - 1][J].cabs1() +
              D[N].cabs1() * X[N][J].cabs1();
        }
      } else {
        if (N == 1) {
          RWORK[1] = B[1][J].cabs1() + D[1].cabs1() * X[1][J].cabs1();
        } else {
          RWORK[1] = B[1][J].cabs1() +
              D[1].cabs1() * X[1][J].cabs1() +
              DL[1].cabs1() * X[2][J].cabs1();
          for (I = 2; I <= N - 1; I++) {
            RWORK[I] = B[I][J].cabs1() +
                DU[I - 1].cabs1() * X[I - 1][J].cabs1() +
                D[I].cabs1() * X[I][J].cabs1() +
                DL[I].cabs1() * X[I + 1][J].cabs1();
          }
          RWORK[N] = B[N][J].cabs1() +
              DU[N - 1].cabs1() * X[N - 1][J].cabs1() +
              D[N].cabs1() * X[N][J].cabs1();
        }
      }

      // Compute componentwise relative backward error from formula

      // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

      // where abs(Z) is the componentwise absolute value of the matrix
      // or vector Z.  If the i-th component of the denominator is less
      // than SAFE2, then SAFE1 is added to the i-th components of the
      // numerator and denominator before dividing.

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

        zgttrs(TRANS, N, 1, DLF, DF, DUF, DU2, IPIV, WORK.asMatrix(N), N, INFO);
        zaxpy(N, Complex.one, WORK, 1, X(1, J).asArray(), 1);
        LSTRES = BERR[J];
        COUNT++;
        continue;
      }
      break;
    }

    // Bound error from formula

    // norm(X - XTRUE) / norm(X) <= FERR =
    // norm( abs(inv(op(A)))*
    //    ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)

    // where
    //   norm(Z) is the magnitude of the largest component of Z
    //   inv(op(A)) is the inverse of op(A)
    //   abs(Z) is the componentwise absolute value of the matrix or
    //      vector Z
    //   NZ is the maximum number of nonzeros in any row of A, plus 1
    //   EPS is machine epsilon

    // The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
    // is incremented by SAFE1 if the i-th component of
    // abs(op(A))*abs(X) + abs(B) is less than SAFE2.

    // Use ZLACN2 to estimate the infinity-norm of the matrix
    //    inv(op(A)) * diag(W),
    // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

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
        // Multiply by diag(W)*inv(op(A)**H).

        zgttrs(
            TRANST, N, 1, DLF, DF, DUF, DU2, IPIV, WORK.asMatrix(N), N, INFO);
        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
      } else {
        // Multiply by inv(op(A))*diag(W).

        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
        zgttrs(
            TRANSN, N, 1, DLF, DF, DUF, DU2, IPIV, WORK.asMatrix(N), N, INFO);
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
