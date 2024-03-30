import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/ztrmv.dart';
import 'package:lapack/src/blas/ztrsv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacn2.dart';

void ztrrfs(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
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
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  bool NOTRAN, NOUNIT, UPPER;
  String TRANSN, TRANST;
  int I, J, K, NZ;
  double EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
  final ISAVE = Array<int>(3);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  NOTRAN = lsame(TRANS, 'N');
  NOUNIT = lsame(DIAG, 'N');

  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (NRHS < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDX < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('ZTRRFS', -INFO.value);
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

  NZ = N + 1;
  EPS = dlamch('Epsilon');
  SAFMIN = dlamch('Safe minimum');
  SAFE1 = NZ * SAFMIN;
  SAFE2 = SAFE1 / EPS;

  // Do for each right hand side

  for (J = 1; J <= NRHS; J++) {
    // Compute residual R = B - op(A) * X,
    // where op(A) = A, A**T, or A**H, depending on TRANS.

    zcopy(N, X(1, J).asArray(), 1, WORK, 1);
    ztrmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1);
    zaxpy(N, -Complex.one, B(1, J).asArray(), 1, WORK, 1);

    // Compute componentwise relative backward error from formula

    // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

    // where abs(Z) is the componentwise absolute value of the matrix
    // or vector Z.  If the i-th component of the denominator is less
    // than SAFE2, then SAFE1 is added to the i-th components of the
    // numerator and denominator before dividing.

    for (I = 1; I <= N; I++) {
      RWORK[I] = B[I][J].cabs1();
    }

    if (NOTRAN) {
      // Compute abs(A)*abs(X) + abs(B).

      if (UPPER) {
        if (NOUNIT) {
          for (K = 1; K <= N; K++) {
            XK = X[K][J].cabs1();
            for (I = 1; I <= K; I++) {
              RWORK[I] += A[I][K].cabs1() * XK;
            }
          }
        } else {
          for (K = 1; K <= N; K++) {
            XK = X[K][J].cabs1();
            for (I = 1; I <= K - 1; I++) {
              RWORK[I] += A[I][K].cabs1() * XK;
            }
            RWORK[K] += XK;
          }
        }
      } else {
        if (NOUNIT) {
          for (K = 1; K <= N; K++) {
            XK = X[K][J].cabs1();
            for (I = K; I <= N; I++) {
              RWORK[I] += A[I][K].cabs1() * XK;
            }
          }
        } else {
          for (K = 1; K <= N; K++) {
            XK = X[K][J].cabs1();
            for (I = K + 1; I <= N; I++) {
              RWORK[I] += A[I][K].cabs1() * XK;
            }
            RWORK[K] += XK;
          }
        }
      }
    } else {
      // Compute abs(A**H)*abs(X) + abs(B).

      if (UPPER) {
        if (NOUNIT) {
          for (K = 1; K <= N; K++) {
            S = ZERO;
            for (I = 1; I <= K; I++) {
              S += A[I][K].cabs1() * X[I][J].cabs1();
            }
            RWORK[K] += S;
          }
        } else {
          for (K = 1; K <= N; K++) {
            S = X[K][J].cabs1();
            for (I = 1; I <= K - 1; I++) {
              S += A[I][K].cabs1() * X[I][J].cabs1();
            }
            RWORK[K] += S;
          }
        }
      } else {
        if (NOUNIT) {
          for (K = 1; K <= N; K++) {
            S = ZERO;
            for (I = K; I <= N; I++) {
              S += A[I][K].cabs1() * X[I][J].cabs1();
            }
            RWORK[K] += S;
          }
        } else {
          for (K = 1; K <= N; K++) {
            S = X[K][J].cabs1();
            for (I = K + 1; I <= N; I++) {
              S += A[I][K].cabs1() * X[I][J].cabs1();
            }
            RWORK[K] += S;
          }
        }
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

        ztrsv(UPLO, TRANST, DIAG, N, A, LDA, WORK, 1);
        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
      } else {
        // Multiply by inv(op(A))*diag(W).

        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
        ztrsv(UPLO, TRANSN, DIAG, N, A, LDA, WORK, 1);
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
