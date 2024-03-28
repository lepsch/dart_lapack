import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zhemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/clag2z.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/cpotrf.dart';
import 'package:lapack/src/cpotrs.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlag2c.dart';
import 'package:lapack/src/zlanhe.dart';
import 'package:lapack/src/zlat2c.dart';
import 'package:lapack/src/zpotrf.dart';
import 'package:lapack/src/zpotrs.dart';

void zcposv(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> WORK_,
  final Array<Complex> SWORK_,
  final Array<double> RWORK_,
  final Box<int> ITER,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having(ld: N);
  final SWORK = SWORK_.having();
  final RWORK = RWORK_.having();
  const DOITREF = false;
  const ITERMAX = 30;
  const BWDMAX = 1.0e+00;
  const NEGONE = Complex(-1.0e+00, 0.0e+00);
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
    xerbla('ZCPOSV', -INFO.value);
    return;
  }

  // Quick return if (N == 0).

  if (N == 0) return;

  skipSinglePrecision:
  while (true) {
    // Skip single precision iterative refinement if a priori slower
    // than double precision factorization.
    if (!DOITREF) {
      ITER.value = -1;
      break skipSinglePrecision;
    }

    // Compute some constants.

    // ignore: dead_code
    ANRM = zlanhe('I', UPLO, N, A, LDA, RWORK);
    EPS = dlamch('Epsilon');
    CTE = ANRM * EPS * sqrt(N) * BWDMAX;

    // Set the indices PTSA, PTSX for referencing SA and SX in SWORK.

    PTSA = 1;
    PTSX = PTSA + N * N;

    // Convert B from double precision to single precision and store the
    // result in SX.

    zlag2c(N, NRHS, B, LDB, SWORK(PTSX).asMatrix(N), N, INFO);

    if (INFO.value != 0) {
      ITER.value = -2;
      break skipSinglePrecision;
    }

    // Convert A from double precision to single precision and store the
    // result in SA.

    zlat2c(UPLO, N, A, LDA, SWORK(PTSA).asMatrix(N), N, INFO);

    if (INFO.value != 0) {
      ITER.value = -2;
      break skipSinglePrecision;
    }

    // Compute the Cholesky factorization of SA.

    cpotrf(UPLO, N, SWORK(PTSA).asMatrix(N), N, INFO);

    if (INFO.value != 0) {
      ITER.value = -3;
      break skipSinglePrecision;
    }

    // Solve the system SA*SX = SB.

    cpotrs(UPLO, N, NRHS, SWORK(PTSA).asMatrix(N), N, SWORK(PTSX).asMatrix(N),
        N, INFO);

    // Convert SX back to Complex

    clag2z(N, NRHS, SWORK(PTSX).asMatrix(N), N, X, LDX, INFO);

    // Compute R = B - AX (R is WORK).

    zlacpy('All', N, NRHS, B, LDB, WORK, N);

    zhemm('Left', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, Complex.one, WORK, N);

    // Check whether the NRHS normwise backward errors satisfy the
    // stopping criterion. If yes, set ITER.value=0 and return.

    var satisfy = true;
    for (I = 1; I <= NRHS; I++) {
      XNRM = X[izamax(N, X(1, I).asArray(), 1)][I].cabs1();
      RNRM = WORK[izamax(N, WORK(1, I).asArray(), 1)][I].cabs1();
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

      zlag2c(N, NRHS, WORK, N, SWORK(PTSX).asMatrix(N), N, INFO);

      if (INFO.value != 0) {
        ITER.value = -2;
        break skipSinglePrecision;
      }

      // Solve the system SA*SX = SR.

      cpotrs(UPLO, N, NRHS, SWORK(PTSA).asMatrix(N), N, SWORK(PTSX).asMatrix(N),
          N, INFO);

      // Convert SX back to double precision and update the current
      // iterate.

      clag2z(N, NRHS, SWORK(PTSX).asMatrix(N), N, WORK, N, INFO);

      for (I = 1; I <= NRHS; I++) {
        zaxpy(N, Complex.one, WORK(1, I).asArray(), 1, X(1, I).asArray(), 1);
      }

      // Compute R = B - AX (R is WORK).

      zlacpy('All', N, NRHS, B, LDB, WORK, N);

      zhemm('L', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, Complex.one, WORK, N);

      // Check whether the NRHS normwise backward errors satisfy the
      // stopping criterion. If yes, set ITER.value=IITER>0 and return.

      var satisfy = true;
      for (I = 1; I <= NRHS; I++) {
        XNRM = X[izamax(N, X(1, I).asArray(), 1)][I].cabs1();
        RNRM = WORK[izamax(N, WORK(1, I).asArray(), 1)][I].cabs1();
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

  zpotrf(UPLO, N, A, LDA, INFO);

  if (INFO.value != 0) return;

  zlacpy('All', N, NRHS, B, LDB, X, LDX);
  zpotrs(UPLO, N, NRHS, A, LDA, X, LDX, INFO);
}
