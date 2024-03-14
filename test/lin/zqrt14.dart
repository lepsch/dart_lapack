import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zgelq2.dart';
import 'package:lapack/src/zgeqr2.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlascl.dart';

import 'xerbla.dart';

double zqrt14(
  final String TRANS,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> X_,
  final int LDX,
  final Array<Complex> WORK_,
  final int LWORK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, ONE = 1.0;
  final RWORK = Array<double>(1);
  final INFO = Box(0);

  final int LDWORK;
  final bool TPSD;
  if (lsame(TRANS, 'N')) {
    LDWORK = M + NRHS;
    TPSD = false;
    if (LWORK < (M + NRHS) * (N + 2)) {
      xerbla('ZQRT14', 10);
      return ZERO;
    } else if (N <= 0 || NRHS <= 0) {
      return ZERO;
    }
  } else if (lsame(TRANS, 'C')) {
    LDWORK = M;
    TPSD = true;
    if (LWORK < (N + NRHS) * (M + 2)) {
      xerbla('ZQRT14', 10);
      return ZERO;
    } else if (M <= 0 || NRHS <= 0) {
      return ZERO;
    }
  } else {
    xerbla('ZQRT14', 1);
    return ZERO;
  }

  // Copy and scale A

  zlacpy('All', M, N, A, LDA, WORK.asMatrix(), LDWORK);
  final ANRM = zlange('M', M, N, WORK.asMatrix(), LDWORK, RWORK);
  if (ANRM != ZERO) {
    zlascl('G', 0, 0, ANRM, ONE, M, N, WORK.asMatrix(), LDWORK, INFO);
  }

  // Copy X or X' into the right place and scale it

  var ERR = ZERO;
  if (TPSD) {
    // Copy X into columns n+1:n+nrhs of work

    zlacpy('All', M, NRHS, X, LDX, WORK(N * LDWORK + 1).asMatrix(), LDWORK);
    final XNRM =
        zlange('M', M, NRHS, WORK(N * LDWORK + 1).asMatrix(), LDWORK, RWORK);
    if (XNRM != ZERO) {
      zlascl('G', 0, 0, XNRM, ONE, M, NRHS, WORK(N * LDWORK + 1).asMatrix(),
          LDWORK, INFO);
    }

    // Compute QR factorization of X

    zgeqr2(M, N + NRHS, WORK.asMatrix(), LDWORK, WORK(LDWORK * (N + NRHS) + 1),
        WORK(LDWORK * (N + NRHS) + min(M, N + NRHS).toInt() + 1), INFO);

    // Compute largest entry in upper triangle of
    // work(n+1:m,n+1:n+nrhs)

    for (var J = N + 1; J <= N + NRHS; J++) {
      for (var I = N + 1; I <= min(M, J); I++) {
        ERR = max(ERR, WORK[I + (J - 1) * M].abs());
      }
    }
  } else {
    // Copy X' into rows m+1:m+nrhs of work

    for (var I = 1; I <= N; I++) {
      for (var J = 1; J <= NRHS; J++) {
        WORK[M + J + (I - 1) * LDWORK] = X[I][J].conjugate();
      }
    }

    final XNRM = zlange('M', NRHS, N, WORK(M + 1).asMatrix(), LDWORK, RWORK);
    if (XNRM != ZERO) {
      zlascl(
          'G', 0, 0, XNRM, ONE, NRHS, N, WORK(M + 1).asMatrix(), LDWORK, INFO);
    }

    // Compute LQ factorization of work

    zgelq2(LDWORK, N, WORK.asMatrix(), LDWORK, WORK(LDWORK * N + 1),
        WORK(LDWORK * (N + 1) + 1), INFO);

    // Compute largest entry in lower triangle in
    // work(m+1:m+nrhs,m+1:n)

    for (var J = M + 1; J <= N; J++) {
      for (var I = J; I <= LDWORK; I++) {
        ERR = max(ERR, WORK[I + (J - 1) * LDWORK].abs());
      }
    }
  }

  return ERR / (max(M, max(N, NRHS)) * dlamch('Epsilon'));
}
