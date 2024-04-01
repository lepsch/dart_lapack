import 'dart:math';

import 'package:lapack/lapack.dart';

double zqrt17(
  final String TRANS,
  final int IRESID,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> C_,
  final Array<Complex> WORK_,
  final int LWORK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDB);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  final RWORK = Array<double>(1);

  final int NROWS, NCOLS;
  if (lsame(TRANS, 'N')) {
    NROWS = M;
    NCOLS = N;
  } else if (lsame(TRANS, 'C')) {
    NROWS = N;
    NCOLS = M;
  } else {
    xerbla('ZQRT17', 1);
    return ZERO;
  }

  if (LWORK < NCOLS * NRHS) {
    xerbla('ZQRT17', 13);
    return ZERO;
  }

  if (M <= 0 || N <= 0 || NRHS <= 0) return ZERO;

  final NORMA = zlange('One-norm', M, N, A, LDA, RWORK);
  final SMLNUM = dlamch('Safe minimum') / dlamch('Precision');

  // compute residual and scale it

  zlacpy('All', NROWS, NRHS, B, LDB, C, LDB);
  zgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, Complex(-ONE), A, LDA, X,
      LDX, Complex.one, C, LDB);
  final int ISCL;
  final NORMRS = zlange('Max', NROWS, NRHS, C, LDB, RWORK);
  if (NORMRS > SMLNUM) {
    ISCL = 1;
    final INFO = Box(0);
    zlascl('General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO);
  } else {
    ISCL = 0;
  }

  // compute R**H * op(A)

  zgemm('Conjugate transpose', TRANS, NRHS, NCOLS, NROWS, Complex.one, C, LDB,
      A, LDA, Complex.zero, WORK.asMatrix(), NRHS);

  // compute and properly scale error

  var ERR = zlange('One-norm', NRHS, NCOLS, WORK.asMatrix(), NRHS, RWORK);
  if (NORMA != ZERO) ERR /= NORMA;

  if (ISCL == 1) ERR *= NORMRS;

  if (IRESID == 1) {
    final NORMB = zlange('One-norm', NROWS, NRHS, B, LDB, RWORK);
    if (NORMB != ZERO) ERR /= NORMB;
  } else {
    if (NORMRS != ZERO) ERR /= NORMRS;
  }

  return ERR / (dlamch('Epsilon') * (max(M, max(N, NRHS))));
}
