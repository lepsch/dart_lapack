import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgbtrs.dart';
import 'package:lapack/src/zlacn2.dart';

double zla_gbrcond_x(
  final String TRANS,
  final int N,
  final int KL,
  final int KU,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> AFB_,
  final int LDAFB,
  final Array<int> IPIV_,
  final Array<Complex> X_,
  final Box<int> INFO,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.dim(LDAB);
  final AFB = AFB_.dim(LDAFB);
  final IPIV = IPIV_.dim();
  final X = X_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  bool NOTRANS;
  int I, J, KD, KE;
  double ANORM, TMP;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0);
  final KASE = Box(0);

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  INFO.value = 0;
  NOTRANS = lsame(TRANS, 'N');
  if (!NOTRANS && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KL < 0 || KL > N - 1) {
    INFO.value = -3;
  } else if (KU < 0 || KU > N - 1) {
    INFO.value = -4;
  } else if (LDAB < KL + KU + 1) {
    INFO.value = -6;
  } else if (LDAFB < 2 * KL + KU + 1) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZLA_GBRCOND_X', -INFO.value);
    return 0;
  }

  // Compute norm of op(A)*op2(C).

  KD = KU + 1;
  KE = KL + 1;
  ANORM = 0.0;
  if (NOTRANS) {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      for (J = max(I - KL, 1); J <= min(I + KU, N); J++) {
        TMP = TMP + CABS1(AB[KD + I - J][J] * X[J]);
      }
      RWORK[I] = TMP;
      ANORM = max(ANORM, TMP);
    }
  } else {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      for (J = max(I - KL, 1); J <= min(I + KU, N); J++) {
        TMP = TMP + CABS1(AB[KE - I + J][I] * X[J]);
      }
      RWORK[I] = TMP;
      ANORM = max(ANORM, TMP);
    }
  }

  // Quick return if possible.

  if (N == 0) {
    return 1.0;
  } else if (ANORM == 0.0) {
    return 0;
  }

  // Estimate the norm of inv(op(A)).

  AINVNM.value = 0.0;

  KASE.value = 0;
  while (true) {
    zlacn2(N, WORK(N + 1), WORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;
    if (KASE.value == 2) {
      // Multiply by R.

      for (I = 1; I <= N; I++) {
        WORK[I] = WORK[I] * RWORK[I].toComplex();
      }

      if (NOTRANS) {
        zgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK.asMatrix(N),
            N, INFO);
      } else {
        zgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV,
            WORK.asMatrix(N), N, INFO);
      }

      // Multiply by inv(X).

      for (I = 1; I <= N; I++) {
        WORK[I] = WORK[I] / X[I];
      }
    } else {
      // Multiply by inv(X**H).

      for (I = 1; I <= N; I++) {
        WORK[I] = WORK[I] / X[I];
      }

      if (NOTRANS) {
        zgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV,
            WORK.asMatrix(N), N, INFO);
      } else {
        zgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK.asMatrix(N),
            N, INFO);
      }

      // Multiply by R.

      for (I = 1; I <= N; I++) {
        WORK[I] = WORK[I] * RWORK[I].toComplex();
      }
    }
  }

  // Compute the estimate of the reciprocal condition number.

  return AINVNM.value != 0.0 ? 1.0 / AINVNM.value : 0;
}
