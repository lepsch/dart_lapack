import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorcsd.dart';
import 'package:lapack/src/dorcsd2by1.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dcsdts(
  final int M,
  final int P,
  final int Q,
  final Matrix<double> X,
  final Matrix<double> XF,
  final int LDX,
  final Matrix<double> U1,
  final int LDU1,
  final Matrix<double> U2,
  final int LDU2,
  final Matrix<double> V1T,
  final int LDV1T,
  final Matrix<double> V2T,
  final int LDV2T,
  final Array<double> THETA,
  final Array<int> IWORK,
  final Array<double> WORK,
  final int LWORK,
  final Array<double> RWORK,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const REALONE = 1.0, REALZERO = 0.0;
  const ZERO = 0.0, ONE = 1.0;
  const PIOVER2 = 1.57079632679489661923132169163975144210;
  int I, R;
  double EPS2, RESID, ULP, ULPINV;
  final INFO = Box(0);

  ULP = dlamch('Precision');
  ULPINV = REALONE / ULP;

  // The first half of the routine checks the 2-by-2 CSD

  dlaset('Full', M, M, ZERO, ONE, WORK.asMatrix(LDX), LDX);
  dsyrk('Upper', 'Conjugate transpose', M, M, -ONE, X, LDX, ONE,
      WORK.asMatrix(LDX), LDX);
  if (M > 0) {
    EPS2 = max(
        ULP, dlange('1', M, M, WORK.asMatrix(LDX), LDX, RWORK) / M.toDouble());
  } else {
    EPS2 = ULP;
  }
  R = min(min(P, M - P), min(Q, M - Q));

  // Copy the matrix X to the array XF.

  dlacpy('Full', M, M, X, LDX, XF, LDX);

  // Compute the CSD

  dorcsd(
      'Y',
      'Y',
      'Y',
      'Y',
      'N',
      'D',
      M,
      P,
      Q,
      XF(1,1),
      LDX,
      XF(1,Q + 1),
      LDX,
      XF(P + 1,1),
      LDX,
      XF(P + 1,Q + 1),
      LDX,
      THETA,
      U1,
      LDU1,
      U2,
      LDU2,
      V1T,
      LDV1T,
      V2T,
      LDV2T,
      WORK,
      LWORK,
      IWORK,
      INFO);

  // Compute XF := diag(U1,U2)'*X*diag(V1,V2) - [D11 D12; D21 D22]

  dlacpy('Full', M, M, X, LDX, XF, LDX);

  dgemm('No transpose', 'Conjugate transpose', P, Q, Q, ONE, XF, LDX, V1T,
      LDV1T, ZERO, WORK.asMatrix(LDX), LDX);

  dgemm('Conjugate transpose', 'No transpose', P, Q, P, ONE, U1, LDU1,
      WORK.asMatrix(LDX), LDX, ZERO, XF, LDX);

  for (I = 1; I <= min(P, Q) - R; I++) {
    XF[I][I] = XF[I][I] - ONE;
  }
  for (I = 1; I <= R; I++) {
    XF[min(P, Q) - R + I][min(P, Q) - R + I] =
        XF[min(P, Q) - R + I][min(P, Q) - R + I] - cos(THETA[I]);
  }

  dgemm('No transpose', 'Conjugate transpose', P, M - Q, M - Q, ONE,
      XF(1, Q + 1), LDX, V2T, LDV2T, ZERO, WORK.asMatrix(LDX), LDX);

  dgemm('Conjugate transpose', 'No transpose', P, M - Q, P, ONE, U1, LDU1,
      WORK.asMatrix(LDX), LDX, ZERO, XF(1, Q + 1), LDX);

  for (I = 1; I <= min(P, M - Q) - R; I++) {
    XF[P - I + 1][M - I + 1] = XF[P - I + 1][M - I + 1] + ONE;
  }
  for (I = 1; I <= R; I++) {
    XF[P - (min(P, M - Q) - R) + 1 - I][M - (min(P, M - Q) - R) + 1 - I] =
        XF[P - (min(P, M - Q) - R) + 1 - I][M - (min(P, M - Q) - R) + 1 - I] +
            sin(THETA[R - I + 1]);
  }

  dgemm('No transpose', 'Conjugate transpose', M - P, Q, Q, ONE, XF(P + 1, 1),
      LDX, V1T, LDV1T, ZERO, WORK.asMatrix(LDX), LDX);

  dgemm('Conjugate transpose', 'No transpose', M - P, Q, M - P, ONE, U2, LDU2,
      WORK.asMatrix(LDX), LDX, ZERO, XF(P + 1, 1), LDX);

  for (I = 1; I <= min(M - P, Q) - R; I++) {
    XF[M - I + 1][Q - I + 1] = XF[M - I + 1][Q - I + 1] - ONE;
  }
  for (I = 1; I <= R; I++) {
    XF[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] =
        XF[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] -
            sin(THETA[R - I + 1]);
  }

  dgemm('No transpose', 'Conjugate transpose', M - P, M - Q, M - Q, ONE,
      XF(P + 1, Q + 1), LDX, V2T, LDV2T, ZERO, WORK.asMatrix(LDX), LDX);

  dgemm('Conjugate transpose', 'No transpose', M - P, M - Q, M - P, ONE, U2,
      LDU2, WORK.asMatrix(LDX), LDX, ZERO, XF(P + 1, Q + 1), LDX);

  for (I = 1; I <= min(M - P, M - Q) - R; I++) {
    XF[P + I][Q + I] = XF[P + I][Q + I] - ONE;
  }
  for (I = 1; I <= R; I++) {
    XF[P + (min(M - P, M - Q) - R) + I][Q + (min(M - P, M - Q) - R) + I] =
        XF[P + (min(M - P, M - Q) - R) + I][Q + (min(M - P, M - Q) - R) + I] -
            cos(THETA[I]);
  }

  // Compute norm( U1'*X11*V1 - D11 ) / ( max(1,P,Q)*EPS2 ) .

  RESID = dlange('1', P, Q, XF, LDX, RWORK);
  RESULT[1] = (RESID / (max(1, max(P, Q)))).toDouble() / EPS2;

  // Compute norm( U1'*X12*V2 - D12 ) / ( max(1,P,M-Q)*EPS2 ) .

  RESID = dlange('1', P, M - Q, XF(1, Q + 1), LDX, RWORK);
  RESULT[2] = (RESID / (max(1, max(P, M - Q)))).toDouble() / EPS2;

  // Compute norm( U2'*X21*V1 - D21 ) / ( max(1,M-P,Q)*EPS2 ) .

  RESID = dlange('1', M - P, Q, XF(P + 1, 1), LDX, RWORK);
  RESULT[3] = (RESID / (max(1, max(M - P, Q)))).toDouble() / EPS2;

  // Compute norm( U2'*X22*V2 - D22 ) / ( max(1,M-P,M-Q)*EPS2 ) .

  RESID = dlange('1', M - P, M - Q, XF(P + 1, Q + 1), LDX, RWORK);
  RESULT[4] = (RESID / (max(1, max(M - P, M - Q)))).toDouble() / EPS2;

  // Compute I - U1'*U1

  dlaset('Full', P, P, ZERO, ONE, WORK.asMatrix(LDU1), LDU1);
  dsyrk('Upper', 'Conjugate transpose', P, P, -ONE, U1, LDU1, ONE,
      WORK.asMatrix(LDU1), LDU1);

  // Compute norm( I - U'*U ) / ( max(1,P) * ULP ) .

  RESID = dlansy('1', 'Upper', P, WORK, LDU1, RWORK);
  RESULT[5] = (RESID / (max(1, P))).toDouble() / ULP;

  // Compute I - U2'*U2

  dlaset('Full', M - P, M - P, ZERO, ONE, WORK.asMatrix(LDU2), LDU2);
  dsyrk('Upper', 'Conjugate transpose', M - P, M - P, -ONE, U2, LDU2, ONE,
      WORK.asMatrix(LDU2), LDU2);

  // Compute norm( I - U2'*U2 ) / ( max(1,M-P) * ULP ) .

  RESID = dlansy('1', 'Upper', M - P, WORK, LDU2, RWORK);
  RESULT[6] = (RESID / (max(1, M - P))).toDouble() / ULP;

  // Compute I - V1T*V1T'

  dlaset('Full', Q, Q, ZERO, ONE, WORK.asMatrix(LDV1T), LDV1T);
  dsyrk('Upper', 'No transpose', Q, Q, -ONE, V1T, LDV1T, ONE,
      WORK.asMatrix(LDV1T), LDV1T);

  // Compute norm( I - V1T*V1T' ) / ( max(1,Q) * ULP ) .

  RESID = dlansy('1', 'Upper', Q, WORK.asMatrix(LDV1T), LDV1T, RWORK);
  RESULT[7] = (RESID / (max(1, Q))).toDouble() / ULP;

  // Compute I - V2T*V2T'

  dlaset('Full', M - Q, M - Q, ZERO, ONE, WORK.asMatrix(LDV2T), LDV2T);
  dsyrk('Upper', 'No transpose', M - Q, M - Q, -ONE, V2T, LDV2T, ONE,
      WORK.asMatrix(LDV2T), LDV2T);

  // Compute norm( I - V2T*V2T' ) / ( max(1,M-Q) * ULP ) .

  RESID = dlansy('1', 'Upper', M - Q, WORK.asMatrix(LDV2T), LDV2T, RWORK);
  RESULT[8] = (RESID / (max(1, M - Q))).toDouble() / ULP;

  // Check sorting

  RESULT[9] = REALZERO;
  for (I = 1; I <= R; I++) {
    if (THETA[I] < REALZERO || THETA[I] > PIOVER2) {
      RESULT[9] = ULPINV;
    }
    if (I > 1) {
      if (THETA[I] < THETA[I - 1]) {
        RESULT[9] = ULPINV;
      }
    }
  }

  // The second half of the routine checks the 2-by-1 CSD

  dlaset('Full', Q, Q, ZERO, ONE, WORK.asMatrix(LDX), LDX);
  dsyrk('Upper', 'Conjugate transpose', Q, M, -ONE, X, LDX, ONE,
      WORK.asMatrix(LDX), LDX);
  if (M > 0) {
    EPS2 = max(
        ULP, dlange('1', Q, Q, WORK.asMatrix(LDX), LDX, RWORK) / M.toDouble());
  } else {
    EPS2 = ULP;
  }
  R = min(min(P, M - P), min(Q, M - Q));

  // Copy the matrix [ X11; X21 ] to the array XF.

  dlacpy('Full', M, Q, X, LDX, XF, LDX);

  // Compute the CSD

  dorcsd2by1('Y', 'Y', 'Y', M, P, Q, XF[1][1], LDX, XF[P + 1][1], LDX, THETA,
      U1, LDU1, U2, LDU2, V1T, LDV1T, WORK, LWORK, IWORK, INFO);

  // Compute [X11;X21] := diag(U1,U2)'*[X11;X21]*V1 - [D11;D21]

  dgemm('No transpose', 'Conjugate transpose', P, Q, Q, ONE, X, LDX, V1T, LDV1T,
      ZERO, WORK.asMatrix(LDX), LDX);

  dgemm('Conjugate transpose', 'No transpose', P, Q, P, ONE, U1, LDU1,
      WORK.asMatrix(LDX), LDX, ZERO, X, LDX);

  for (I = 1; I <= min(P, Q) - R; I++) {
    X[I][I] = X[I][I] - ONE;
  }
  for (I = 1; I <= R; I++) {
    X[min(P, Q) - R + I][min(P, Q) - R + I] =
        X[min(P, Q) - R + I][min(P, Q) - R + I] - cos(THETA[I]);
  }

  dgemm('No transpose', 'Conjugate transpose', M - P, Q, Q, ONE, X(P + 1, 1),
      LDX, V1T, LDV1T, ZERO, WORK.asMatrix(LDX), LDX);

  dgemm('Conjugate transpose', 'No transpose', M - P, Q, M - P, ONE, U2, LDU2,
      WORK.asMatrix(LDX), LDX, ZERO, X(P + 1, 1), LDX);

  for (I = 1; I <= min(M - P, Q) - R; I++) {
    X[M - I + 1][Q - I + 1] = X[M - I + 1][Q - I + 1] - ONE;
  }
  for (I = 1; I <= R; I++) {
    X[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] =
        X[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] -
            sin(THETA[R - I + 1]);
  }

  // Compute norm( U1'*X11*V1 - D11 ) / ( max(1,P,Q)*EPS2 ) .

  RESID = dlange('1', P, Q, X, LDX, RWORK);
  RESULT[10] = (RESID / (max(1, max(P, Q)))).toDouble() / EPS2;

  // Compute norm( U2'*X21*V1 - D21 ) / ( max(1,M-P,Q)*EPS2 ) .

  RESID = dlange('1', M - P, Q, X(P + 1, 1), LDX, RWORK);
  RESULT[11] = (RESID / (max(1, max(M - P, Q)))).toDouble() / EPS2;

  // Compute I - U1'*U1

  dlaset('Full', P, P, ZERO, ONE, WORK.asMatrix(LDU1), LDU1);
  dsyrk('Upper', 'Conjugate transpose', P, P, -ONE, U1, LDU1, ONE,
      WORK.asMatrix(LDU1), LDU1);

  // Compute norm( I - U1'*U1 ) / ( max(1,P) * ULP ) .

  RESID = dlansy('1', 'Upper', P, WORK.asMatrix(LDU1), LDU1, RWORK);
  RESULT[12] = (RESID / (max(1, P))).toDouble() / ULP;

  // Compute I - U2'*U2

  dlaset('Full', M - P, M - P, ZERO, ONE, WORK.asMatrix(LDU2), LDU2);
  dsyrk('Upper', 'Conjugate transpose', M - P, M - P, -ONE, U2, LDU2, ONE,
      WORK.asMatrix(LDU2), LDU2);

  // Compute norm( I - U2'*U2 ) / ( max(1,M-P) * ULP ) .

  RESID = dlansy('1', 'Upper', M - P, WORK.asMatrix(LDU2), LDU2, RWORK);
  RESULT[13] = (RESID / (max(1, M - P))).toDouble() / ULP;

  // Compute I - V1T*V1T'

  dlaset('Full', Q, Q, ZERO, ONE, WORK.asMatrix(LDV1T), LDV1T);
  dsyrk('Upper', 'No transpose', Q, Q, -ONE, V1T, LDV1T, ONE,
      WORK.asMatrix(LDV1T), LDV1T);

  // Compute norm( I - V1T*V1T' ) / ( max(1,Q) * ULP ) .

  RESID = dlansy('1', 'Upper', Q, WORK.asMatrix(LDV1T), LDV1T, RWORK);
  RESULT[14] = (RESID / (max(1, Q))).toDouble() / ULP;

  // Check sorting

  RESULT[15] = REALZERO;
  for (I = 1; I <= R; I++) {
    if (THETA[I] < REALZERO || THETA[I] > PIOVER2) {
      RESULT[15] = ULPINV;
    }
    if (I > 1) {
      if (THETA[I] < THETA[I - 1]) {
        RESULT[15] = ULPINV;
      }
    }
  }
}
