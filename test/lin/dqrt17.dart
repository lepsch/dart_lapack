import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

double dqrt17(
  final String TRANS,
  final int IRESID,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> C_,
  final Array<double> WORK_,
  final int LWORK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDB);
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, ONE = 1.0;
  final RWORK = Array<double>(1);

  final int NROWS, NCOLS;
  if (lsame(TRANS, 'N')) {
    NROWS = M;
    NCOLS = N;
  } else if (lsame(TRANS, 'T')) {
    NROWS = N;
    NCOLS = M;
  } else {
    xerbla('DQRT17', 1);
    return ZERO;
  }

  if (LWORK < NCOLS * NRHS) {
    xerbla('DQRT17', 13);
    return ZERO;
  }

  if (M <= 0 || N <= 0 || NRHS <= 0) {
    return ZERO;
  }

  final NORMA = dlange('One-norm', M, N, A, LDA, RWORK);
  final SMLNUM = dlamch('Safe minimum') / dlamch('Precision');

  // compute residual and scale it

  dlacpy('All', NROWS, NRHS, B, LDB, C, LDB);
  dgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, -ONE, A, LDA, X, LDX, ONE, C,
      LDB);
  final NORMRS = dlange('Max', NROWS, NRHS, C, LDB, RWORK);
  final int ISCL;
  if (NORMRS > SMLNUM) {
    ISCL = 1;
    final INFO = Box(0);
    dlascl('General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO);
  } else {
    ISCL = 0;
  }

  // compute R**T * op(A)

  dgemm('Transpose', TRANS, NRHS, NCOLS, NROWS, ONE, C, LDB, A, LDA, ZERO,
      WORK.asMatrix(), NRHS);

  // compute and properly scale error

  var ERR = dlange('One-norm', NRHS, NCOLS, WORK.asMatrix(), NRHS, RWORK);
  if (NORMA != ZERO) ERR /= NORMA;

  if (ISCL == 1) ERR *= NORMRS;

  if (IRESID == 1) {
    final NORMB = dlange('One-norm', NROWS, NRHS, B, LDB, RWORK);
    if (NORMB != ZERO) ERR /= NORMB;
  } else {
    if (NORMRS != ZERO) ERR /= NORMRS;
  }

  return ERR / (dlamch('Epsilon') * max(M, max(N, NRHS)));
}
