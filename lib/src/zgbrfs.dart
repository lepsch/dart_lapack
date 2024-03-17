import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgbmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgbtrs.dart';
import 'package:lapack/src/zlacn2.dart';

void zgbrfs(
  final String TRANS,
  final int N,
  final int KL,
  final int KU,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> AFB_,
  final int LDAFB,
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
  final AB = AB_.having(ld: LDAB);
  final AFB = AFB_.having(ld: LDAFB);
  final IPIV = IPIV_.having();
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
  bool NOTRAN;
  String TRANSN, TRANST;
  int COUNT, I, J, K, KK, NZ;
  double EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
  final ISAVE = Array<int>(3);
  final KASE = Box(0);

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  // Test the input parameters.

  INFO.value = 0;
  NOTRAN = lsame(TRANS, 'N');
  if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KL < 0) {
    INFO.value = -3;
  } else if (KU < 0) {
    INFO.value = -4;
  } else if (NRHS < 0) {
    INFO.value = -5;
  } else if (LDAB < KL + KU + 1) {
    INFO.value = -7;
  } else if (LDAFB < 2 * KL + KU + 1) {
    INFO.value = -9;
  } else if (LDB < max(1, N)) {
    INFO.value = -12;
  } else if (LDX < max(1, N)) {
    INFO.value = -14;
  }
  if (INFO.value != 0) {
    xerbla('ZGBRFS', -INFO.value);
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

  NZ = min(KL + KU + 2, N + 1);
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
      zgbmv(TRANS, N, N, KL, KU, -Complex.one, AB, LDAB, X(1, J).asArray(), 1,
          Complex.one, WORK, 1);

      // Compute componentwise relative backward error from formula

      // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

      // where abs(Z) is the componentwise absolute value of the matrix
      // or vector Z.  If the i-th component of the denominator is less
      // than SAFE2, then SAFE1 is added to the i-th components of the
      // numerator and denominator before dividing.

      for (I = 1; I <= N; I++) {
        RWORK[I] = CABS1(B[I][J]);
      }

      // Compute abs(op(A))*abs(X) + abs(B).

      if (NOTRAN) {
        for (K = 1; K <= N; K++) {
          KK = KU + 1 - K;
          XK = CABS1(X[K][J]);
          for (I = max(1, K - KU); I <= min(N, K + KL); I++) {
            RWORK[I] += CABS1(AB[KK + I][K]) * XK;
          }
        }
      } else {
        for (K = 1; K <= N; K++) {
          S = ZERO;
          KK = KU + 1 - K;
          for (I = max(1, K - KU); I <= min(N, K + KL); I++) {
            S += CABS1(AB[KK + I][K]) * CABS1(X[I][J]);
          }
          RWORK[K] += S;
        }
      }
      S = ZERO;
      for (I = 1; I <= N; I++) {
        if (RWORK[I] > SAFE2) {
          S = max(S, CABS1(WORK[I]) / RWORK[I]);
        } else {
          S = max(S, (CABS1(WORK[I]) + SAFE1) / (RWORK[I] + SAFE1));
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

        zgbtrs(
            TRANS, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK.asMatrix(N), N, INFO);
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
        RWORK[I] = CABS1(WORK[I]) + NZ * EPS * RWORK[I];
      } else {
        RWORK[I] = CABS1(WORK[I]) + NZ * EPS * RWORK[I] + SAFE1;
      }
    }

    KASE.value = 0;
    while (true) {
      zlacn2(N, WORK(N + 1), WORK, FERR.box(J), KASE, ISAVE);
      if (KASE.value == 0) break;
      if (KASE.value == 1) {
        // Multiply by diag(W)*inv(op(A)**H).

        zgbtrs(
            TRANST, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK.asMatrix(N), N, INFO);
        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
      } else {
        // Multiply by inv(op(A))*diag(W).

        for (I = 1; I <= N; I++) {
          WORK[I] = RWORK[I].toComplex() * WORK[I];
        }
        zgbtrs(
            TRANSN, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK.asMatrix(N), N, INFO);
      }
    }

    // Normalize error.

    LSTRES = ZERO;
    for (I = 1; I <= N; I++) {
      LSTRES = max(LSTRES, CABS1(X[I][J]));
    }
    if (LSTRES != ZERO) FERR[J] /= LSTRES;
  }
}
