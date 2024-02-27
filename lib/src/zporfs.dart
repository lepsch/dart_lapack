import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zhemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacn2.dart';
import 'package:lapack/src/zpotrs.dart';

void zporfs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
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
  final A = A_.dim(LDA);
  final AF = AF_.dim(LDAF);
  final B = B_.dim(LDB);
  final X = X_.dim(LDX);
  final WORK = WORK_.dim();
  final FERR = FERR_.dim();
  final BERR = BERR_.dim();
  final RWORK = RWORK_.dim();
  const ITMAX = 5;
  const ZERO = 0.0;
  const TWO = 2.0;
  const THREE = 3.0;
  bool UPPER;
  int COUNT, I, J, K, NZ;
  double EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
  final ISAVE = Array<int>(3);
  final KASE = Box(0);

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDAF < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDX < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('ZPORFS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) {
    for (J = 1; J <= NRHS; J++) {
      // 10
      FERR[J] = ZERO;
      BERR[J] = ZERO;
    } // 10
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
    // 140

    COUNT = 1;
    LSTRES = THREE;
    while (true) {
      // Loop until stopping criterion is satisfied.

      // Compute residual R = B - A * X

      zcopy(N, B(1, J).asArray(), 1, WORK, 1);
      zhemv(UPLO, N, -Complex.one, A, LDA, X(1, J).asArray(), 1, Complex.one,
          WORK, 1);

      // Compute componentwise relative backward error from formula

      // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

      // where abs(Z) is the componentwise absolute value of the matrix
      // or vector Z.  If the i-th component of the denominator is less
      // than SAFE2, then SAFE1 is added to the i-th components of the
      // numerator and denominator before dividing.

      for (I = 1; I <= N; I++) {
        // 30
        RWORK[I] = CABS1(B[I][J]);
      } // 30

      // Compute abs(A)*abs(X) + abs(B).

      if (UPPER) {
        for (K = 1; K <= N; K++) {
          // 50
          S = ZERO;
          XK = CABS1(X[K][J]);
          for (I = 1; I <= K - 1; I++) {
            // 40
            RWORK[I] = RWORK[I] + CABS1(A[I][K]) * XK;
            S = S + CABS1(A[I][K]) * CABS1(X[I][J]);
          } // 40
          RWORK[K] = RWORK[K] + A[K][K].toDouble().abs() * XK + S;
        } // 50
      } else {
        for (K = 1; K <= N; K++) {
          // 70
          S = ZERO;
          XK = CABS1(X[K][J]);
          RWORK[K] = RWORK[K] + A[K][K].toDouble().abs() * XK;
          for (I = K + 1; I <= N; I++) {
            // 60
            RWORK[I] = RWORK[I] + CABS1(A[I][K]) * XK;
            S = S + CABS1(A[I][K]) * CABS1(X[I][J]);
          } // 60
          RWORK[K] = RWORK[K] + S;
        } // 70
      }
      S = ZERO;
      for (I = 1; I <= N; I++) {
        // 80
        if (RWORK[I] > SAFE2) {
          S = max(S, CABS1(WORK[I]) / RWORK[I]);
        } else {
          S = max(S, (CABS1(WORK[I]) + SAFE1) / (RWORK[I] + SAFE1));
        }
      } // 80
      BERR[J] = S;

      // Test stopping criterion. Continue iterating if
      //    1) The residual BERR(J) is larger than machine epsilon, and
      //    2) BERR(J) decreased by at least a factor of 2 during the
      //       last iteration, and
      //    3) At most ITMAX iterations tried.

      if (BERR[J] > EPS && TWO * BERR[J] <= LSTRES && COUNT <= ITMAX) {
        // Update solution and try again.

        zpotrs(UPLO, N, 1, AF, LDAF, WORK.asMatrix(), N, INFO);
        zaxpy(N, Complex.one, WORK, 1, X(1, J).asArray(), 1);
        LSTRES = BERR[J];
        COUNT = COUNT + 1;
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
      // 90
      if (RWORK[I] > SAFE2) {
        RWORK[I] = CABS1(WORK[I]) + NZ * EPS * RWORK[I];
      } else {
        RWORK[I] = CABS1(WORK[I]) + NZ * EPS * RWORK[I] + SAFE1;
      }
    } // 90

    KASE.value = 0;
    while (true) {
      zlacn2(N, WORK(N + 1), WORK, FERR(J), KASE, ISAVE);
      if (KASE.value == 0) break;
      if (KASE.value == 1) {
        // Multiply by diag(W)*inv(A**H).

        zpotrs(UPLO, N, 1, AF, LDAF, WORK.asMatrix(), N, INFO);
        for (I = 1; I <= N; I++) {
          // 110
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        } // 110
      } else if (KASE.value == 2) {
        // Multiply by inv(A)*diag(W).

        for (I = 1; I <= N; I++) {
          // 120
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        } // 120
        zpotrs(UPLO, N, 1, AF, LDAF, WORK.asMatrix(), N, INFO);
      }
    }

    // Normalize error.

    LSTRES = ZERO;
    for (I = 1; I <= N; I++) {
      // 130
      LSTRES = max(LSTRES, CABS1(X[I][J]));
    } // 130
    if (LSTRES != ZERO) FERR[J] = FERR[J] / LSTRES;
  } // 140
}
