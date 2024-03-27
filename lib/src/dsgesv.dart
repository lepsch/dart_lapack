import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgetrf.dart';
import 'package:lapack/src/dgetrs.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlag2s.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/sgetrf.dart';
import 'package:lapack/src/sgetrs.dart';
import 'package:lapack/src/slag2d.dart';
import 'package:lapack/src/xerbla.dart';

void dsgesv(
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> WORK_,
  final Array<double> SWORK_,
  final Box<int> ITER,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having(ld: N);
  final SWORK = SWORK_.having();
  const DOITREF = true;
  const ITERMAX = 30;
  const BWDMAX = 1.0e+00;
  const NEGONE = -1.0, ONE = 1.0;
  int I, IITER, PTSA, PTSX;
  double ANRM, CTE, EPS, RNRM, XNRM;

  INFO.value = 0;
  ITER.value = 0;

  // Test the input parameters.

  if (N < 0) {
    INFO.value = -1;
  } else if (NRHS < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  } else if (LDX < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DSGESV', -INFO.value);
    return;
  }

  // Quick return if (N == 0).

  if (N == 0) return;

  skipSinglePrecision:
  while (true) {
    // Skip single precision iterative refinement if a priori slower
    // than double precision factorization.

    // ignore: dead_code
    if (!DOITREF) {
      ITER.value = -1;
      break;
    }

    // Compute some constants.

    ANRM = dlange('I', N, N, A, LDA, WORK.asArray());
    EPS = dlamch('Epsilon');
    CTE = ANRM * EPS * sqrt(N) * BWDMAX;

    // Set the indices PTSA, PTSX for referencing SA and SX in SWORK.

    PTSA = 1;
    PTSX = PTSA + N * N;

    // Convert B from double precision to single precision and store the
    // result in SX.

    dlag2s(N, NRHS, B, LDB, SWORK(PTSX).asMatrix(N), N, INFO);

    if (INFO.value != 0) {
      ITER.value = -2;
      break;
    }

    // Convert A from double precision to single precision and store the
    // result in SA.

    dlag2s(N, N, A, LDA, SWORK(PTSA).asMatrix(N), N, INFO);

    if (INFO.value != 0) {
      ITER.value = -2;
      break;
    }

    // Compute the LU factorization of SA.

    sgetrf(N, N, SWORK(PTSA).asMatrix(N), N, IPIV, INFO);

    if (INFO.value != 0) {
      ITER.value = -3;
      break;
    }

    // Solve the system SA*SX = SB.

    sgetrs('No transpose', N, NRHS, SWORK(PTSA).asMatrix(N), N, IPIV,
        SWORK(PTSX).asMatrix(N), N, INFO);

    // Convert SX back to double precision

    slag2d(N, NRHS, SWORK(PTSX).asMatrix(N), N, X, LDX, INFO);

    // Compute R = B - AX (R is WORK).

    dlacpy('All', N, NRHS, B, LDB, WORK, N);

    dgemm('No Transpose', 'No Transpose', N, NRHS, N, NEGONE, A, LDA, X, LDX,
        ONE, WORK, N);

    // Check whether the NRHS normwise backward errors satisfy the
    // stopping criterion. If yes, set ITER=0 and return.
    var stopped = false;
    for (I = 1; I <= NRHS; I++) {
      XNRM = X[idamax(N, X(1, I).asArray(), 1)][I].abs();
      RNRM = WORK[idamax(N, WORK(1, I).asArray(), 1)][I].abs();
      if (RNRM > XNRM * CTE) {
        stopped = true;
        break;
      }
    }

    if (!stopped) {
      // If we are here, the NRHS normwise backward errors satisfy the
      // stopping criterion. We are good to exit.

      ITER.value = 0;
      return;
    }

    for (IITER = 1; IITER <= ITERMAX; IITER++) {
      // Convert R (in WORK) from double precision to single precision
      // and store the result in SX.

      dlag2s(N, NRHS, WORK, N, SWORK(PTSX).asMatrix(N), N, INFO);

      if (INFO.value != 0) {
        ITER.value = -2;
        break skipSinglePrecision;
      }

      // Solve the system SA*SX = SR.

      sgetrs('No transpose', N, NRHS, SWORK(PTSA).asMatrix(N), N, IPIV,
          SWORK(PTSX).asMatrix(N), N, INFO);

      // Convert SX back to double precision and update the current
      // iterate.

      slag2d(N, NRHS, SWORK(PTSX).asMatrix(N), N, WORK, N, INFO);

      for (I = 1; I <= NRHS; I++) {
        daxpy(N, ONE, WORK(1, I).asArray(), 1, X(1, I).asArray(), 1);
      }

      // Compute R = B - AX (R is WORK).

      dlacpy('All', N, NRHS, B, LDB, WORK, N);

      dgemm('No Transpose', 'No Transpose', N, NRHS, N, NEGONE, A, LDA, X, LDX,
          ONE, WORK, N);

      // Check whether the NRHS normwise backward errors satisfy the
      // stopping criterion. If yes, set ITER.value=IITER>0 and return.

      var stopped = false;
      for (I = 1; I <= NRHS; I++) {
        XNRM = X[idamax(N, X(1, I).asArray(), 1)][I];
        RNRM = WORK[idamax(N, WORK(1, I).asArray(), 1)][I];
        if (RNRM > XNRM * CTE) {
          stopped = true;
          break;
        }
      }

      if (!stopped) {
        // If we are here, the NRHS normwise backward errors satisfy the
        // stopping criterion, we are good to exit.

        ITER.value = IITER;

        return;
      }
    }

    // If we are at this place of the code, this is because we have
    // performed ITER.value=ITERMAX iterations and never satisfied the
    // stopping criterion, set up the ITER.value flag accordingly and follow up
    // on double precision routine.

    ITER.value = -ITERMAX - 1;
  }

  // Single-precision iterative refinement failed to converge to a
  // satisfactory solution, so we resort to double precision.

  dgetrf(N, N, A, LDA, IPIV, INFO);

  if (INFO.value != 0) return;

  dlacpy('All', N, NRHS, B, LDB, X, LDX);
  dgetrs('No transpose', N, NRHS, A, LDA, IPIV, X, LDX, INFO);
}
