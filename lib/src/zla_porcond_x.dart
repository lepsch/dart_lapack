import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacn2.dart';
import 'package:lapack/src/zpotrs.dart';

double zla_porcond_x(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
  final Array<Complex> X_,
  final Box<int> INFO,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final AF = AF_.dim(LDAF);
  final X = X_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  int I, J;
  double ANORM, TMP;
  bool UP, UPPER;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0);
  final KASE = Box(0);

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LDAF < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZLA_PORCOND_X', -INFO.value);
    return 0;
  }
  UP = false;
  if (lsame(UPLO, 'U')) UP = true;

  // Compute norm of op(A)*op2(C).

  ANORM = 0.0;
  if (UP) {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      for (J = 1; J <= I; J++) {
        TMP = TMP + CABS1(A[J][I] * X[J]);
      }
      for (J = I + 1; J <= N; J++) {
        TMP = TMP + CABS1(A[I][J] * X[J]);
      }
      RWORK[I] = TMP;
      ANORM = max(ANORM, TMP);
    }
  } else {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      for (J = 1; J <= I; J++) {
        TMP = TMP + CABS1(A[I][J] * X[J]);
      }
      for (J = I + 1; J <= N; J++) {
        TMP = TMP + CABS1(A[J][I] * X[J]);
      }
      RWORK[I] = TMP;
      ANORM = max(ANORM, TMP);
    }
  }

  // Quick return if possible.

  if (N == 0) {
    return 1;
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

      if (UP) {
        zpotrs('U', N, 1, AF, LDAF, WORK.asMatrix(N), N, INFO);
      } else {
        zpotrs('L', N, 1, AF, LDAF, WORK.asMatrix(N), N, INFO);
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

      if (UP) {
        zpotrs('U', N, 1, AF, LDAF, WORK.asMatrix(N), N, INFO);
      } else {
        zpotrs('L', N, 1, AF, LDAF, WORK.asMatrix(N), N, INFO);
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
