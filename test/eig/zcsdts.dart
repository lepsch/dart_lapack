import 'dart:math';

import 'package:collection/collection.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlanhe.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zuncsd.dart';
import 'package:lapack/src/zuncsd2by1.dart';

void zcsdts(
  final int M,
  final int P,
  final int Q,
  final Matrix<Complex> X_,
  final Matrix<Complex> XF_,
  final int LDX,
  final Matrix<Complex> U1_,
  final int LDU1,
  final Matrix<Complex> U2_,
  final int LDU2,
  final Matrix<Complex> V1T_,
  final int LDV1T,
  final Matrix<Complex> V2T_,
  final int LDV2T,
  final Array<double> THETA_,
  final Array<int> IWORK_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.dim(LDU1);
  final XF = XF_.dim(LDU1);
  final U1 = U1_.dim(LDU1);
  final U2 = U2_.dim(LDU2);
  final V1T = V1T_.dim(LDV1T);
  final V2T = V2T_.dim(LDV2T);
  final THETA = THETA_.dim();
  final IWORK = IWORK_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final RESULT = RESULT_.dim(15);

  const REALONE = 1.0, REALZERO = 0.0;
  const PIOVER2 = 1.57079632679489661923132169163975144210;
  int I, R;
  double EPS2, RESID, ULP, ULPINV;
  final INFO = Box(0);

  ULP = dlamch('Precision');
  ULPINV = REALONE / ULP;

  // The first half of the routine checks the 2-by-2 CSD

  zlaset('Full', M, M, Complex.zero, Complex.one, WORK.asMatrix(), LDX);
  zherk('Upper', 'Conjugate transpose', M, M, -REALONE, X, LDX, REALONE,
      WORK.asMatrix(), LDX);
  if (M > 0) {
    EPS2 =
        max(ULP, zlange('1', M, M, WORK.asMatrix(), LDX, RWORK) / M.toDouble());
  } else {
    EPS2 = ULP;
  }
  R = [P, M - P, Q, M - Q].min;

  // Copy the matrix X to the array XF.

  zlacpy('Full', M, M, X, LDX, XF, LDX);

  // Compute the CSD

  zuncsd(
      'Y',
      'Y',
      'Y',
      'Y',
      'N',
      'D',
      M,
      P,
      Q,
      XF(1, 1),
      LDX,
      XF(1, Q + 1),
      LDX,
      XF(P + 1, 1),
      LDX,
      XF(P + 1, Q + 1),
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
      RWORK,
      17 * (R + 2),
      IWORK,
      INFO);

  // Compute XF := diag(U1,U2)'*X*diag(V1,V2) - [D11 D12; D21 D22]

  zlacpy('Full', M, M, X, LDX, XF, LDX);

  zgemm('No transpose', 'Conjugate transpose', P, Q, Q, Complex.one, XF, LDX,
      V1T, LDV1T, Complex.zero, WORK.asMatrix(), LDX);

  zgemm('Conjugate transpose', 'No transpose', P, Q, P, Complex.one, U1, LDU1,
      WORK.asMatrix(), LDX, Complex.zero, XF, LDX);

  for (I = 1; I <= min(P, Q) - R; I++) {
    XF[I][I] = XF[I][I] - Complex.one;
  }
  for (I = 1; I <= R; I++) {
    XF[min(P, Q) - R + I][min(P, Q) - R + I] =
        XF[min(P, Q) - R + I][min(P, Q) - R + I] - cos(THETA[I]).toComplex();
  }

  zgemm('No transpose', 'Conjugate transpose', P, M - Q, M - Q, Complex.one,
      XF(1, Q + 1), LDX, V2T, LDV2T, Complex.zero, WORK.asMatrix(), LDX);

  zgemm('Conjugate transpose', 'No transpose', P, M - Q, P, Complex.one, U1,
      LDU1, WORK.asMatrix(), LDX, Complex.zero, XF(1, Q + 1), LDX);

  for (I = 1; I <= min(P, M - Q) - R; I++) {
    XF[P - I + 1][M - I + 1] = XF[P - I + 1][M - I + 1] + Complex.one;
  }
  for (I = 1; I <= R; I++) {
    XF[P - (min(P, M - Q) - R) + 1 - I][M - (min(P, M - Q) - R) + 1 - I] =
        XF[P - (min(P, M - Q) - R) + 1 - I][M - (min(P, M - Q) - R) + 1 - I] +
            sin(THETA[R - I + 1]).toComplex();
  }

  zgemm('No transpose', 'Conjugate transpose', M - P, Q, Q, Complex.one,
      XF(P + 1, 1), LDX, V1T, LDV1T, Complex.zero, WORK.asMatrix(), LDX);

  zgemm('Conjugate transpose', 'No transpose', M - P, Q, M - P, Complex.one, U2,
      LDU2, WORK.asMatrix(), LDX, Complex.zero, XF(P + 1, 1), LDX);

  for (I = 1; I <= min(M - P, Q) - R; I++) {
    XF[M - I + 1][Q - I + 1] = XF[M - I + 1][Q - I + 1] - Complex.one;
  }
  for (I = 1; I <= R; I++) {
    XF[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] =
        XF[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] -
            sin(THETA[R - I + 1]).toComplex();
  }

  zgemm('No transpose', 'Conjugate transpose', M - P, M - Q, M - Q, Complex.one,
      XF(P + 1, Q + 1), LDX, V2T, LDV2T, Complex.zero, WORK.asMatrix(), LDX);

  zgemm('Conjugate transpose', 'No transpose', M - P, M - Q, M - P, Complex.one,
      U2, LDU2, WORK.asMatrix(), LDX, Complex.zero, XF(P + 1, Q + 1), LDX);

  for (I = 1; I <= min(M - P, M - Q) - R; I++) {
    XF[P + I][Q + I] = XF[P + I][Q + I] - Complex.one;
  }
  for (I = 1; I <= R; I++) {
    XF[P + (min(M - P, M - Q) - R) + I][Q + (min(M - P, M - Q) - R) + I] =
        XF[P + (min(M - P, M - Q) - R) + I][Q + (min(M - P, M - Q) - R) + I] -
            cos(THETA[I]).toComplex();
  }

  // Compute norm( U1'*X11*V1 - D11 ) / ( max(1,P,Q)*EPS2 ) .

  RESID = zlange('1', P, Q, XF, LDX, RWORK);
  RESULT[1] = (RESID / (max(1, max(P, Q)))) / EPS2;

  // Compute norm( U1'*X12*V2 - D12 ) / ( max(1,P,M-Q)*EPS2 ) .

  RESID = zlange('1', P, M - Q, XF(1, Q + 1), LDX, RWORK);
  RESULT[2] = (RESID / (max(1, max(P, M - Q)))) / EPS2;

  // Compute norm( U2'*X21*V1 - D21 ) / ( max(1,M-P,Q)*EPS2 ) .

  RESID = zlange('1', M - P, Q, XF(P + 1, 1), LDX, RWORK);
  RESULT[3] = (RESID / (max(1, max(M - P, Q)))) / EPS2;

  // Compute norm( U2'*X22*V2 - D22 ) / ( max(1,M-P,M-Q)*EPS2 ) .

  RESID = zlange('1', M - P, M - Q, XF(P + 1, Q + 1), LDX, RWORK);
  RESULT[4] = (RESID / (max(1, max(M - P, M - Q)))) / EPS2;

  // Compute I - U1'*U1

  zlaset('Full', P, P, Complex.zero, Complex.one, WORK.asMatrix(), LDU1);
  zherk('Upper', 'Conjugate transpose', P, P, -REALONE, U1, LDU1, REALONE,
      WORK.asMatrix(), LDU1);

  // Compute norm( I - U'*U ) / ( max(1,P) * ULP ) .

  RESID = zlanhe('1', 'Upper', P, WORK.asMatrix(), LDU1, RWORK);
  RESULT[5] = (RESID / (max(1, P))) / ULP;

  // Compute I - U2'*U2

  zlaset(
      'Full', M - P, M - P, Complex.zero, Complex.one, WORK.asMatrix(), LDU2);
  zherk('Upper', 'Conjugate transpose', M - P, M - P, -REALONE, U2, LDU2,
      REALONE, WORK.asMatrix(), LDU2);

  // Compute norm( I - U2'*U2 ) / ( max(1,M-P) * ULP ) .

  RESID = zlanhe('1', 'Upper', M - P, WORK.asMatrix(), LDU2, RWORK);
  RESULT[6] = (RESID / (max(1, M - P))) / ULP;

  // Compute I - V1T*V1T'

  zlaset('Full', Q, Q, Complex.zero, Complex.one, WORK.asMatrix(), LDV1T);
  zherk('Upper', 'No transpose', Q, Q, -REALONE, V1T, LDV1T, REALONE,
      WORK.asMatrix(), LDV1T);

  // Compute norm( I - V1T*V1T' ) / ( max(1,Q) * ULP ) .

  RESID = zlanhe('1', 'Upper', Q, WORK.asMatrix(), LDV1T, RWORK);
  RESULT[7] = (RESID / (max(1, Q))) / ULP;

  // Compute I - V2T*V2T'

  zlaset(
      'Full', M - Q, M - Q, Complex.zero, Complex.one, WORK.asMatrix(), LDV2T);
  zherk('Upper', 'No transpose', M - Q, M - Q, -REALONE, V2T, LDV2T, REALONE,
      WORK.asMatrix(), LDV2T);

  // Compute norm( I - V2T*V2T' ) / ( max(1,M-Q) * ULP ) .

  RESID = zlanhe('1', 'Upper', M - Q, WORK.asMatrix(), LDV2T, RWORK);
  RESULT[8] = (RESID / (max(1, M - Q))) / ULP;

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

  zlaset('Full', Q, Q, Complex.zero, Complex.one, WORK.asMatrix(), LDX);
  zherk('Upper', 'Conjugate transpose', Q, M, -REALONE, X, LDX, REALONE,
      WORK.asMatrix(), LDX);
  if (M > 0) {
    EPS2 =
        max(ULP, zlange('1', Q, Q, WORK.asMatrix(), LDX, RWORK) / M.toDouble());
  } else {
    EPS2 = ULP;
  }
  R = min(min(P, M - P), min(Q, M - Q));

  // Copy the matrix X to the array XF.

  zlacpy('Full', M, M, X, LDX, XF, LDX);

  // Compute the CSD

  zuncsd2by1(
      'Y',
      'Y',
      'Y',
      M,
      P,
      Q,
      XF(1, 1),
      LDX,
      XF(P + 1, 1),
      LDX,
      THETA,
      U1,
      LDU1,
      U2,
      LDU2,
      V1T,
      LDV1T,
      WORK,
      LWORK,
      RWORK,
      17 * (R + 2),
      IWORK,
      INFO);

  // Compute [X11;X21] := diag(U1,U2)'*[X11;X21]*V1 - [D11;D21]

  zgemm('No transpose', 'Conjugate transpose', P, Q, Q, Complex.one, X, LDX,
      V1T, LDV1T, Complex.zero, WORK.asMatrix(), LDX);

  zgemm('Conjugate transpose', 'No transpose', P, Q, P, Complex.one, U1, LDU1,
      WORK.asMatrix(), LDX, Complex.zero, X, LDX);

  for (I = 1; I <= min(P, Q) - R; I++) {
    X[I][I] = X[I][I] - Complex.one;
  }
  for (I = 1; I <= R; I++) {
    X[min(P, Q) - R + I][min(P, Q) - R + I] =
        X[min(P, Q) - R + I][min(P, Q) - R + I] - cos(THETA[I]).toComplex();
  }

  zgemm('No transpose', 'Conjugate transpose', M - P, Q, Q, Complex.one,
      X(P + 1, 1), LDX, V1T, LDV1T, Complex.zero, WORK.asMatrix(), LDX);

  zgemm('Conjugate transpose', 'No transpose', M - P, Q, M - P, Complex.one, U2,
      LDU2, WORK.asMatrix(), LDX, Complex.zero, X(P + 1, 1), LDX);

  for (I = 1; I <= min(M - P, Q) - R; I++) {
    X[M - I + 1][Q - I + 1] = X[M - I + 1][Q - I + 1] - Complex.one;
  }
  for (I = 1; I <= R; I++) {
    X[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] =
        X[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] -
            sin(THETA[R - I + 1]).toComplex();
  }

  // Compute norm( U1'*X11*V1 - D11 ) / ( max(1,P,Q)*EPS2 ) .

  RESID = zlange('1', P, Q, X, LDX, RWORK);
  RESULT[10] = (RESID / (max(1, max(P, Q)))) / EPS2;

  // Compute norm( U2'*X21*V1 - D21 ) / ( max(1,M-P,Q)*EPS2 ) .

  RESID = zlange('1', M - P, Q, X(P + 1, 1), LDX, RWORK);
  RESULT[11] = (RESID / (max(1, max(M - P, Q)))) / EPS2;

  // Compute I - U1'*U1

  zlaset('Full', P, P, Complex.zero, Complex.one, WORK.asMatrix(), LDU1);
  zherk('Upper', 'Conjugate transpose', P, P, -REALONE, U1, LDU1, REALONE,
      WORK.asMatrix(), LDU1);

  // Compute norm( I - U'*U ) / ( max(1,P) * ULP ) .

  RESID = zlanhe('1', 'Upper', P, WORK.asMatrix(), LDU1, RWORK);
  RESULT[12] = (RESID / (max(1, P))) / ULP;

  // Compute I - U2'*U2

  zlaset(
      'Full', M - P, M - P, Complex.zero, Complex.one, WORK.asMatrix(), LDU2);
  zherk('Upper', 'Conjugate transpose', M - P, M - P, -REALONE, U2, LDU2,
      REALONE, WORK.asMatrix(), LDU2);

  // Compute norm( I - U2'*U2 ) / ( max(1,M-P) * ULP ) .

  RESID = zlanhe('1', 'Upper', M - P, WORK.asMatrix(), LDU2, RWORK);
  RESULT[13] = (RESID / (max(1, M - P))) / ULP;

  // Compute I - V1T*V1T'

  zlaset('Full', Q, Q, Complex.zero, Complex.one, WORK.asMatrix(), LDV1T);
  zherk('Upper', 'No transpose', Q, Q, -REALONE, V1T, LDV1T, REALONE,
      WORK.asMatrix(), LDV1T);

  // Compute norm( I - V1T*V1T' ) / ( max(1,Q) * ULP ) .

  RESID = zlanhe('1', 'Upper', Q, WORK.asMatrix(), LDV1T, RWORK);
  RESULT[14] = (RESID / (max(1, Q))) / ULP;

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
