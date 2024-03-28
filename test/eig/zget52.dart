import 'dart:math';

import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';

void zget52(
  final bool LEFT,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> E_,
  final int LDE,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final E = E_.having(ld: LDE);
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having(length: 2);

  const ZERO = 0.0, ONE = 1.0;
  String NORMAB, TRANS;
  int J, JVEC;
  double ABMAX,
      ALFMAX,
      ANORM,
      BETMAX,
      BNORM,
      ENORM,
      ENRMER,
      ERRNRM,
      SAFMAX,
      SAFMIN,
      SCALE,
      TEMP1,
      ULP;
  Complex ACOEFF, ALPHAI, BCOEFF, BETAI;

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0) return;

  SAFMIN = dlamch('Safe minimum');
  SAFMAX = ONE / SAFMIN;
  ULP = dlamch('Epsilon') * dlamch('Base');

  if (LEFT) {
    TRANS = 'C';
    NORMAB = 'I';
  } else {
    TRANS = 'N';
    NORMAB = 'O';
  }

  // Norm of A, B, and E:

  ANORM = max(zlange(NORMAB, N, N, A, LDA, RWORK), SAFMIN);
  BNORM = max(zlange(NORMAB, N, N, B, LDB, RWORK), SAFMIN);
  ENORM = max(zlange('O', N, N, E, LDE, RWORK), ULP);
  ALFMAX = SAFMAX / max(ONE, BNORM);
  BETMAX = SAFMAX / max(ONE, ANORM);

  // Compute error matrix.
  // Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )

  for (JVEC = 1; JVEC <= N; JVEC++) {
    ALPHAI = ALPHA[JVEC];
    BETAI = BETA[JVEC];
    ABMAX = max(ALPHAI.cabs1(), BETAI.cabs1());
    if (ALPHAI.cabs1() > ALFMAX || BETAI.cabs1() > BETMAX || ABMAX < ONE) {
      SCALE = ONE / max(ABMAX, SAFMIN);
      ALPHAI = SCALE.toComplex() * ALPHAI;
      BETAI = SCALE.toComplex() * BETAI;
    }
    SCALE =
        ONE / max(ALPHAI.cabs1() * BNORM, max(BETAI.cabs1() * ANORM, SAFMIN));
    ACOEFF = SCALE.toComplex() * BETAI;
    BCOEFF = SCALE.toComplex() * ALPHAI;
    if (LEFT) {
      ACOEFF = ACOEFF.conjugate();
      BCOEFF = BCOEFF.conjugate();
    }
    zgemv(TRANS, N, N, ACOEFF, A, LDA, E(1, JVEC).asArray(), 1, Complex.zero,
        WORK(N * (JVEC - 1) + 1), 1);
    zgemv(TRANS, N, N, -BCOEFF, B, LDA, E(1, JVEC).asArray(), 1, Complex.one,
        WORK(N * (JVEC - 1) + 1), 1);
  }

  ERRNRM = zlange('One', N, N, WORK.asMatrix(), N, RWORK) / ENORM;

  // Compute RESULT(1)

  RESULT[1] = ERRNRM / ULP;

  // Normalization of E:

  ENRMER = ZERO;
  for (JVEC = 1; JVEC <= N; JVEC++) {
    TEMP1 = ZERO;
    for (J = 1; J <= N; J++) {
      TEMP1 = max(TEMP1, E[J][JVEC].cabs1());
    }
    ENRMER = max(ENRMER, (TEMP1 - ONE).abs());
  }

  // Compute RESULT(2) : the normalization error in E.

  RESULT[2] = ENRMER / (N * ULP);
}
