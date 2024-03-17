import 'dart:math';

import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zbdt03(
  final String UPLO,
  final int N,
  final int KD,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> U_,
  final int LDU,
  final Array<double> S_,
  final Matrix<Complex> VT_,
  final int LDVT,
  final Array<Complex> WORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final S = S_.having();
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double BNORM, EPS;

  // Quick return if possible

  RESID.value = ZERO;
  if (N <= 0) return;

  // Compute B - U * S * V' one column at a time.

  BNORM = ZERO;
  if (KD >= 1) {
    // B is bidiagonal.

    if (lsame(UPLO, 'U')) {
      // B is upper bidiagonal.

      for (J = 1; J <= N; J++) {
        // 20
        for (I = 1; I <= N; I++) {
          // 10
          WORK[N + I] = S[I].toComplex() * VT[I][J];
        } // 10
        zgemv('No transpose', N, N, -Complex.one, U, LDU, WORK(N + 1), 1,
            Complex.zero, WORK, 1);
        WORK[J] += D[J].toComplex();
        if (J > 1) {
          WORK[J - 1] += E[J - 1].toComplex();
          BNORM = max(BNORM, (D[J]).abs() + (E[J - 1]).abs());
        } else {
          BNORM = max(BNORM, (D[J]).abs());
        }
        RESID.value = max(RESID.value, dzasum(N, WORK, 1));
      } // 20
    } else {
      // B is lower bidiagonal.

      for (J = 1; J <= N; J++) {
        // 40
        for (I = 1; I <= N; I++) {
          // 30
          WORK[N + I] = S[I].toComplex() * VT[I][J];
        } // 30
        zgemv('No transpose', N, N, -Complex.one, U, LDU, WORK(N + 1), 1,
            Complex.zero, WORK, 1);
        WORK[J] += D[J].toComplex();
        if (J < N) {
          WORK[J + 1] += E[J].toComplex();
          BNORM = max(BNORM, (D[J]).abs() + (E[J]).abs());
        } else {
          BNORM = max(BNORM, (D[J]).abs());
        }
        RESID.value = max(RESID.value, dzasum(N, WORK, 1));
      } // 40
    }
  } else {
    // B is diagonal.

    for (J = 1; J <= N; J++) {
      // 60
      for (I = 1; I <= N; I++) {
        // 50
        WORK[N + I] = S[I].toComplex() * VT[I][J];
      } // 50
      zgemv('No transpose', N, N, -Complex.one, U, LDU, WORK(N + 1), 1,
          Complex.zero, WORK, 1);
      WORK[J] += D[J].toComplex();
      RESID.value = max(RESID.value, dzasum(N, WORK, 1));
    } // 60
    J = idamax(N, D, 1);
    BNORM = (D[J]).abs();
  }

  // Compute norm(B - U * S * V') / ( n * norm(B) * EPS )

  EPS = dlamch('Precision');

  if (BNORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    if (BNORM >= RESID.value) {
      RESID.value = (RESID.value / BNORM) / (N.toDouble() * EPS);
    } else {
      if (BNORM < ONE) {
        RESID.value = (min(RESID.value, (N).toDouble() * BNORM) / BNORM) /
            (N.toDouble() * EPS);
      } else {
        RESID.value =
            min(RESID.value / BNORM, (N).toDouble()) / (N.toDouble() * EPS);
      }
    }
  }
}
