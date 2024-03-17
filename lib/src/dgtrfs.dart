import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgttrs.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/dlagtm.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgtrfs(
  final String TRANS,
  final int N,
  final int NRHS,
  final Array<double> DL_,
  final Array<double> D_,
  final Array<double> DU_,
  final Array<double> DLF_,
  final Array<double> DF_,
  final Array<double> DUF_,
  final Array<double> DU2_,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> X_,
  final int LDX,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DLF = DLF_.having();
  final DF = DF_.having();
  final DUF = DUF_.having();
  final DU2 = DU2_.having();
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ITMAX = 5;
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0;
  const THREE = 3.0;
  bool NOTRAN;
  String TRANSN, TRANST;
  int COUNT, I, J, NZ;
  double EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN;
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
    xerbla('DGTRFS', -INFO.value);
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
    TRANST = 'T';
  } else {
    TRANSN = 'T';
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

      dcopy(N, B(1, J).asArray(), 1, WORK(N + 1), 1);
      dlagtm(TRANS, N, 1, -ONE, DL, D, DU, X(1, J), LDX, ONE,
          WORK(N + 1).asMatrix(N), N);

      // Compute abs(op(A))*abs(x) + abs(b) for use in the backward
      // error bound.

      if (NOTRAN) {
        if (N == 1) {
          WORK[1] = (B[1][J]).abs() + (D[1] * X[1][J]).abs();
        } else {
          WORK[1] = (B[1][J]).abs() +
              (D[1] * X[1][J]).abs() +
              (DU[1] * X[2][J]).abs();
          for (I = 2; I <= N - 1; I++) {
            WORK[I] = (B[I][J]).abs() +
                (DL[I - 1] * X[I - 1][J]).abs() +
                (D[I] * X[I][J]).abs() +
                (DU[I] * X[I + 1][J]).abs();
          }
          WORK[N] = (B[N][J]).abs() +
              (DL[N - 1] * X[N - 1][J]).abs() +
              (D[N] * X[N][J]).abs();
        }
      } else {
        if (N == 1) {
          WORK[1] = (B[1][J]).abs() + (D[1] * X[1][J]).abs();
        } else {
          WORK[1] = (B[1][J]).abs() +
              (D[1] * X[1][J]).abs() +
              (DL[1] * X[2][J]).abs();
          for (I = 2; I <= N - 1; I++) {
            WORK[I] = (B[I][J]).abs() +
                (DU[I - 1] * X[I - 1][J]).abs() +
                (D[I] * X[I][J]).abs() +
                (DL[I] * X[I + 1][J]).abs();
          }
          WORK[N] = (B[N][J]).abs() +
              (DU[N - 1] * X[N - 1][J]).abs() +
              (D[N] * X[N][J]).abs();
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
        if (WORK[I] > SAFE2) {
          S = max(S, (WORK[N + I]).abs() / WORK[I]);
        } else {
          S = max(S, ((WORK[N + I]).abs() + SAFE1) / (WORK[I] + SAFE1));
        }
      }
      BERR[J] = S;

      // Test stopping criterion. Continue iterating if
      //    1) The residual BERR[J] is larger than machine epsilon, and
      //    2) BERR[J] decreased by at least a factor of 2 during the
      //       last iteration, and
      //    3) At most ITMAX iterations tried.

      if (BERR[J] > EPS && TWO * BERR[J] <= LSTRES && COUNT <= ITMAX) {
        // Update solution and try again.

        dgttrs(TRANS, N, 1, DLF, DF, DUF, DU2, IPIV, WORK(N + 1).asMatrix(N), N,
            INFO);
        daxpy(N, ONE, WORK(N + 1), 1, X(1, J).asArray(), 1);
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

    // Use DLACN2 to estimate the infinity-norm of the matrix
    //    inv(op(A)) * diag(W),
    // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

    for (I = 1; I <= N; I++) {
      if (WORK[I] > SAFE2) {
        WORK[I] = (WORK[N + I]).abs() + NZ * EPS * WORK[I];
      } else {
        WORK[I] = (WORK[N + I]).abs() + NZ * EPS * WORK[I] + SAFE1;
      }
    }

    KASE.value = 0;
    while (true) {
      dlacn2(N, WORK(2 * N + 1), WORK(N + 1), IWORK, FERR.box(J), KASE, ISAVE);
      if (KASE.value == 0) break;

      if (KASE.value == 1) {
        // Multiply by diag(W)*inv(op(A)**T).

        dgttrs(TRANST, N, 1, DLF, DF, DUF, DU2, IPIV, WORK(N + 1).asMatrix(N),
            N, INFO);
        for (I = 1; I <= N; I++) {
          WORK[N + I] = WORK[I] * WORK[N + I];
        }
      } else {
        // Multiply by inv(op(A))*diag(W).

        for (I = 1; I <= N; I++) {
          WORK[N + I] = WORK[I] * WORK[N + I];
        }
        dgttrs(TRANSN, N, 1, DLF, DF, DUF, DU2, IPIV, WORK(N + 1).asMatrix(N),
            N, INFO);
      }
    }

    // Normalize error.

    LSTRES = ZERO;
    for (I = 1; I <= N; I++) {
      LSTRES = max(LSTRES, (X[I][J]).abs());
    }
    if (LSTRES != ZERO) FERR[J] /= LSTRES;
  }
}
