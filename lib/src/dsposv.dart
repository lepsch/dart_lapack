import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dsymm.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlag2s.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlat2s.dart';
import 'package:lapack/src/dpotrs.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/slag2d.dart';
import 'package:lapack/src/spotrf.dart';
import 'package:lapack/src/spotrs.dart';
import 'package:lapack/src/dpotrf.dart';
import 'package:lapack/src/xerbla.dart';

void dsposv(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
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
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final SWORK = SWORK_.having();
  final WORK = WORK_.having(ld: N);

  const DOITREF = true;
  const ITERMAX = 30;
  const BWDMAX = 1.0e+00;
  const NEGONE = -1.0, ONE = 1.0;
  int I, IITER, PTSA, PTSX;
  double ANRM, CTE, EPS, RNRM, XNRM;

  INFO.value = 0;
  ITER.value = 0;

  // Test the input parameters.

  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  } else if (LDX < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DSPOSV', -INFO.value);
    return;
  }

  // Quick return if (N == 0).

  if (N == 0) return;

  fallbackToFullPrecision:
  while (true) {
    // Skip single precision iterative refinement if a priori slower
    // than double precision factorization.

    // ignore: dead_code
    if (!DOITREF) {
      ITER.value = -1;
      break fallbackToFullPrecision;
    }

    // Compute some constants.

    ANRM = dlansy('I', UPLO, N, A, LDA, WORK.asArray());
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
      break fallbackToFullPrecision;
    }

    // Convert A from double precision to single precision and store the
    // result in SA.

    dlat2s(UPLO, N, A, LDA, SWORK(PTSA).asMatrix(N), N, INFO);

    if (INFO.value != 0) {
      ITER.value = -2;
      break fallbackToFullPrecision;
    }

    // Compute the Cholesky factorization of SA.

    spotrf(UPLO, N, SWORK(PTSA).asMatrix(N), N, INFO);

    if (INFO.value != 0) {
      ITER.value = -3;
      break fallbackToFullPrecision;
    }

    // Solve the system SA*SX = SB.

    spotrs(UPLO, N, NRHS, SWORK(PTSA).asMatrix(N), N, SWORK(PTSX).asMatrix(N),
        N, INFO);

    // Convert SX back to double precision

    slag2d(N, NRHS, SWORK(PTSX).asMatrix(N), N, X, LDX, INFO);

    // Compute R = B - AX (R is WORK).

    dlacpy('All', N, NRHS, B, LDB, WORK, N);

    dsymm('Left', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, WORK, N);

    // Check whether the NRHS normwise backward errors satisfy the
    // stopping criterion. If yes, set ITER=0 and return.
    var satisfy = true;
    for (I = 1; I <= NRHS; I++) {
      XNRM = X[idamax(N, X(1, I).asArray(), 1)][I].abs();
      RNRM = WORK[idamax(N, WORK(1, I).asArray(), 1)][I].abs();
      if (RNRM > XNRM * CTE) {
        satisfy = false;
        break;
      }
    }

    if (satisfy) {
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
        break fallbackToFullPrecision;
      }

      // Solve the system SA*SX = SR.

      spotrs(UPLO, N, NRHS, SWORK(PTSA).asMatrix(N), N, SWORK(PTSX).asMatrix(N),
          N, INFO);

      // Convert SX back to double precision and update the current
      // iterate.

      slag2d(N, NRHS, SWORK(PTSX).asMatrix(N), N, WORK, N, INFO);

      for (I = 1; I <= NRHS; I++) {
        daxpy(N, ONE, WORK(1, I).asArray(), 1, X(1, I).asArray(), 1);
      }

      // Compute R = B - AX (R is WORK).

      dlacpy('All', N, NRHS, B, LDB, WORK, N);

      dsymm('L', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, WORK, N);

      // Check whether the NRHS normwise backward errors satisfy the
      // stopping criterion. If yes, set ITER=IITER>0 and return.

      var satisfy = true;
      for (I = 1; I <= NRHS; I++) {
        XNRM = X[idamax(N, X(1, I).asArray(), 1)][I].abs();
        RNRM = WORK[idamax(N, WORK(1, I).asArray(), 1)][I].abs();
        if (RNRM > XNRM * CTE) {
          satisfy = false;
          break;
        }
      }

      if (satisfy) {
        // If we are here, the NRHS normwise backward errors satisfy the
        // stopping criterion, we are good to exit.

        ITER.value = IITER;
        return;
      }
    }

    // If we are at this place of the code, this is because we have
    // performed ITER=ITERMAX iterations and never satisfied the
    // stopping criterion, set up the ITER flag accordingly and follow
    // up on double precision routine.

    ITER.value = -ITERMAX - 1;

    break;
  }

  // Single-precision iterative refinement failed to converge to a
  // satisfactory solution, so we resort to double precision.

  dpotrf(UPLO, N, A, LDA, INFO);

  if (INFO.value != 0) return;

  dlacpy('All', N, NRHS, B, LDB, X, LDX);
  dpotrs(UPLO, N, NRHS, A, LDA, X, LDX, INFO);
}
