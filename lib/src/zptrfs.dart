import 'dart:math';

import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zpttrs.dart';

void zptrfs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<Complex> E_,
  final Array<double> DF_,
  final Array<Complex> EF_,
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
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final E = E_.having();
  final EF = EF_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final D = D_.having();
  final DF = DF_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  const ITMAX = 5;
  const ZERO = 0.0;
  const ONE = 1.0;
  const TWO = 2.0;
  const THREE = 3.0;
  bool UPPER;
  int COUNT = 0, I, IX, J, NZ;
  double EPS, LSTRES = 0, S, SAFE1, SAFE2, SAFMIN;
  Complex BI, CX, DX, EX;

  double CABS1(Complex ZDUM) => ZDUM.real.abs() + ZDUM.imaginary.abs();

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
    INFO.value = -9;
  } else if (LDX < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('ZPTRFS', -INFO.value);
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

      // Compute residual R = B - A * X.  Also compute
      // abs(A)*abs(x) + abs(b) for use in the backward error bound.

      if (UPPER) {
        if (N == 1) {
          BI = B[1][J];
          DX = D[1].toComplex() * X[1][J];
          WORK[1] = BI - DX;
          RWORK[1] = CABS1(BI) + CABS1(DX);
        } else {
          BI = B[1][J];
          DX = D[1].toComplex() * X[1][J];
          EX = E[1] * X[2][J];
          WORK[1] = BI - DX - EX;
          RWORK[1] = CABS1(BI) + CABS1(DX) + CABS1(E[1]) * CABS1(X[2][J]);
          for (I = 2; I <= N - 1; I++) {
            BI = B[I][J];
            CX = E[I - 1].conjugate() * X[I - 1][J];
            DX = D[I].toComplex() * X[I][J];
            EX = E[I] * X[I + 1][J];
            WORK[I] = BI - CX - DX - EX;
            RWORK[I] = CABS1(BI) +
                CABS1(E[I - 1]) * CABS1(X[I - 1][J]) +
                CABS1(DX) +
                CABS1(E[I]) * CABS1(X[I + 1][J]);
          }
          BI = B[N][J];
          CX = E[N - 1].conjugate() * X[N - 1][J];
          DX = D[N].toComplex() * X[N][J];
          WORK[N] = BI - CX - DX;
          RWORK[N] =
              CABS1(BI) + CABS1(E[N - 1]) * CABS1(X[N - 1][J]) + CABS1(DX);
        }
      } else {
        if (N == 1) {
          BI = B[1][J];
          DX = D[1].toComplex() * X[1][J];
          WORK[1] = BI - DX;
          RWORK[1] = CABS1(BI) + CABS1(DX);
        } else {
          BI = B[1][J];
          DX = D[1].toComplex() * X[1][J];
          EX = E[1].conjugate() * X[2][J];
          WORK[1] = BI - DX - EX;
          RWORK[1] = CABS1(BI) + CABS1(DX) + CABS1(E[1]) * CABS1(X[2][J]);
          for (I = 2; I <= N - 1; I++) {
            BI = B[I][J];
            CX = E[I - 1] * X[I - 1][J];
            DX = D[I].toComplex() * X[I][J];
            EX = E[I].conjugate() * X[I + 1][J];
            WORK[I] = BI - CX - DX - EX;
            RWORK[I] = CABS1(BI) +
                CABS1(E[I - 1]) * CABS1(X[I - 1][J]) +
                CABS1(DX) +
                CABS1(E[I]) * CABS1(X[I + 1][J]);
          }
          BI = B[N][J];
          CX = E[N - 1] * X[N - 1][J];
          DX = D[N].toComplex() * X[N][J];
          WORK[N] = BI - CX - DX;
          RWORK[N] =
              CABS1(BI) + CABS1(E[N - 1]) * CABS1(X[N - 1][J]) + CABS1(DX);
        }
      }

      // Compute componentwise relative backward error from formula

      // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

      // where abs(Z) is the componentwise absolute value of the matrix
      // or vector Z.  If the i-th component of the denominator is less
      // than SAFE2, then SAFE1 is added to the i-th components of the
      // numerator and denominator before dividing.

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

        zpttrs(UPLO, N, 1, DF, EF, WORK.asMatrix(), N, INFO);
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

    for (I = 1; I <= N; I++) {
      if (RWORK[I] > SAFE2) {
        RWORK[I] = CABS1(WORK[I]) + NZ * EPS * RWORK[I];
      } else {
        RWORK[I] = CABS1(WORK[I]) + NZ * EPS * RWORK[I] + SAFE1;
      }
    }
    IX = idamax(N, RWORK, 1);
    FERR[J] = RWORK[IX];

    // Estimate the norm of inv(A).

    // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

    // m(i,j) =  abs(A(i,j)), i = j,
    // m(i,j) = -abs(A(i,j)), i != j,

    // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.

    // Solve M(L) * x = e.

    RWORK[1] = ONE;
    for (I = 2; I <= N; I++) {
      RWORK[I] = ONE + RWORK[I - 1] * EF[I - 1].abs();
    }

    // Solve D * M(L)**H * x = b.

    RWORK[N] /= DF[N];
    for (I = N - 1; I >= 1; I--) {
      RWORK[I] /= DF[I] + RWORK[I + 1] * EF[I].abs();
    }

    // Compute norm(inv(A)) = max(x(i)), 1<=i<=n.

    IX = idamax(N, RWORK, 1);
    FERR[J] *= RWORK[IX].abs();

    // Normalize error.

    LSTRES = ZERO;
    for (I = 1; I <= N; I++) {
      LSTRES = max(LSTRES, X[I][J].abs());
    }
    if (LSTRES != ZERO) FERR[J] /= LSTRES;
  }
}
