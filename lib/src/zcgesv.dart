import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/cgetrf.dart';
import 'package:lapack/src/cgetrs.dart';
import 'package:lapack/src/clag2z.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgetrf.dart';
import 'package:lapack/src/zgetrs.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlag2c.dart';
import 'package:lapack/src/zlange.dart';

void zcgesv(
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
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
  final A = A_.dim(LDA);
  final IPIV = IPIV_.dim();
  final B = B_.dim(LDB);
  final X = X_.dim(LDX);
  final SWORK = SWORK_.dim(N);
  final WORK = WORK_.dim(N);
  final RWORK = RWORK_.dim();
  const DOITREF = true;
  const ITERMAX = 30;
  const BWDMAX = 1.0e+00;
  const NEGONE = Complex(-1.0e+00, 0.0e+00), ONE = Complex.one;
  int I, IITER, PTSA, PTSX;
  double ANRM, CTE, EPS, RNRM, XNRM;
  // Complex ZDUM;
  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

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
    xerbla('ZCGESV', -INFO.value);
    return;
  }

  // Quick return if (N == 0).

  if (N == 0) return;

  // Skip single precision iterative refinement if a priori slower
  // than double precision factorization.

  skipSinglePrecision:
  while (true) {
    // ignore: dead_code
    if (!DOITREF) {
      ITER.value = -1;
      break skipSinglePrecision;
    }

    // Compute some constants.

    ANRM = zlange('I', N, N, A, LDA, RWORK);
    EPS = dlamch('Epsilon');
    CTE = ANRM * EPS * sqrt(N.toDouble()) * BWDMAX;

    // Set the indices PTSA, PTSX for referencing SA and SX in SWORK.

    PTSA = 1;
    PTSX = PTSA + N * N;

    // Convert B from double precision to single precision and store the
    // result in SX.

    zlag2c(N, NRHS, B, LDB, SWORK(PTSX), N, INFO);

    if (INFO.value != 0) {
      ITER.value = -2;
      break skipSinglePrecision;
    }

    // Convert A from double precision to single precision and store the
    // result in SA.

    zlag2c(N, N, A, LDA, SWORK(PTSA), N, INFO);

    if (INFO.value != 0) {
      ITER.value = -2;
      break skipSinglePrecision;
    }

    // Compute the LU factorization of SA.

    cgetrf(N, N, SWORK(PTSA).asMatrix(N), N, IPIV, INFO);

    if (INFO.value != 0) {
      ITER.value = -3;
      break skipSinglePrecision;
    }

    // Solve the system SA*SX = SB.

    cgetrs('No transpose', N, NRHS, SWORK(PTSA).asMatrix(N), N, IPIV, SWORK(PTSX).asMatrix(N), N, INFO);

    // Convert SX back to double precision

    clag2z(N, NRHS, SWORK(PTSX).asMatrix(N), N, X, LDX, INFO);

    // Compute R = B - AX (R is WORK).

    zlacpy('All', N, NRHS, B, LDB, WORK, N);

    zgemm('No Transpose', 'No Transpose', N, NRHS, N, NEGONE, A, LDA, X, LDX,
        ONE, WORK, N);

    // Check whether the NRHS normwise backward errors satisfy the
    // stopping criterion. If yes, set ITER.value=0 and return.

    var satisfy = true;
    for (I = 1; I <= NRHS; I++) {
      XNRM = CABS1(X[izamax(N, X(1, I).asArray(), 1)][I]);
      RNRM = CABS1(WORK[izamax(N, WORK(1, I).asArray(), 1)][I]);
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
      // 30

      // Convert R (in WORK) from double precision to single precision
      // and store the result in SX.

      zlag2c(N, NRHS, WORK, N, SWORK(PTSX), N, INFO);

      if (INFO.value != 0) {
        ITER.value = -2;
        break skipSinglePrecision;
      }

      // Solve the system SA*SX = SR.

      cgetrs('No transpose', N, NRHS, SWORK(PTSA).asMatrix(N), N, IPIV,
          SWORK(PTSX).asMatrix(N), N, INFO);

      // Convert SX back to double precision and update the current
      // iterate.

      clag2z(N, NRHS, SWORK(PTSX).asMatrix(N), N, WORK, N, INFO);

      for (I = 1; I <= NRHS; I++) {
        zaxpy(N, ONE, WORK(1, I).asArray(), 1, X(1, I).asArray(), 1);
      }

      // Compute R = B - AX (R is WORK).

      zlacpy('All', N, NRHS, B, LDB, WORK, N);

      zgemm('No Transpose', 'No Transpose', N, NRHS, N, NEGONE, A, LDA, X, LDX,
          ONE, WORK, N);

      // Check whether the NRHS normwise backward errors satisfy the
      // stopping criterion. If yes, set ITER.value=IITER>0 and return.
      var satisfy = true;
      for (I = 1; I <= NRHS; I++) {
        XNRM = CABS1(X[izamax(N, X(1, I).asArray(), 1)][I]);
        RNRM = CABS1(WORK[izamax(N, WORK(1, I).asArray(), 1)][I]);
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
    } // 30

    // If we are at this place of the code, this is because we have
    // performed ITER.value=ITERMAX iterations and never satisfied the stopping
    // criterion, set up the ITER.value flag accordingly and follow up on double
    // precision routine.

    ITER.value = -ITERMAX - 1;
  } // 40

  // Single-precision iterative refinement failed to converge to a
  // satisfactory solution, so we resort to double precision.

  zgetrf(N, N, A, LDA, IPIV, INFO);

  if (INFO.value != 0) return;

  zlacpy('All', N, NRHS, B, LDB, X, LDX);
  zgetrs('No transpose', N, NRHS, A, LDA, IPIV, X, LDX, INFO);
}
