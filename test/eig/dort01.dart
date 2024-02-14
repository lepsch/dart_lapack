import 'dart:math';

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dort01(
  final String ROWCOL,
  final int M,
  final int N,
  final Matrix<double> U_,
  final int LDU,
  final Array<double> WORK_,
  final int LWORK,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final U = U_.dim(LDU);
  final WORK = WORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  String TRANSU;
  int I, J, K, LDWORK, MNMIN;
  double EPS, TMP;

  RESID.value = ZERO;

  // Quick return if possible

  if (M <= 0 || N <= 0) return;

  EPS = dlamch('Precision');
  if (M < N || (M == N && lsame(ROWCOL, 'R'))) {
    TRANSU = 'N';
    K = N;
  } else {
    TRANSU = 'T';
    K = M;
  }
  MNMIN = min(M, N);

  if ((MNMIN + 1) * MNMIN <= LWORK) {
    LDWORK = MNMIN;
  } else {
    LDWORK = 0;
  }
  if (LDWORK > 0) {
    // Compute I - U*U' or I - U'*U.

    dlaset('Upper', MNMIN, MNMIN, ZERO, ONE, WORK.asMatrix(LDWORK), LDWORK);
    dsyrk('Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK.asMatrix(LDWORK),
        LDWORK);

    // Compute norm( I - U*U' ) / ( K * EPS ) .

    RESID.value = dlansy('1', 'Upper', MNMIN, WORK.asMatrix(LDWORK), LDWORK,
        WORK(LDWORK * MNMIN + 1));
    RESID.value = (RESID.value / K.toDouble()) / EPS;
  } else if (TRANSU == 'T') {
    // Find the maximum element in abs( I - U'*U ) / ( m * EPS )

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= J; I++) {
        if (I != J) {
          TMP = ZERO;
        } else {
          TMP = ONE;
        }
        TMP = TMP - ddot(M, U(1, I).asArray(), 1, U(1, J).asArray(), 1);
        RESID.value = max(RESID.value, (TMP).abs());
      }
    }
    RESID.value = (RESID.value / M.toDouble()) / EPS;
  } else {
    // Find the maximum element in abs( I - U*U' ) / ( n * EPS )

    for (J = 1; J <= M; J++) {
      for (I = 1; I <= J; I++) {
        if (I != J) {
          TMP = ZERO;
        } else {
          TMP = ONE;
        }
        TMP = TMP - ddot(N, U(J, 1).asArray(), LDU, U(I, 1).asArray(), LDU);
        RESID.value = max(RESID.value, (TMP).abs());
      }
    }
    RESID.value = (RESID.value / N.toDouble()) / EPS;
  }
}
