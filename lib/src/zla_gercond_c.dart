import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgetrs.dart';
import 'package:lapack/src/zlacn2.dart';

double zla_gercond_c(
  final String TRANS,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
  final Array<int> IPIV_,
  final Array<double> C_,
  final bool CAPPLY,
  final Box<int> INFO,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final AF = AF_.dim(LDAF);
  final IPIV = IPIV_.dim();
  final WORK = WORK_.dim();
  final C = C_.dim();
  final RWORK = RWORK_.dim();
  bool NOTRANS;
  int I, J;
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
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LDAF < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZLA_GERCOND_C', -INFO.value);
    return 0;
  }

  // Compute norm of op(A)*op2(C).

  ANORM = 0.0;
  if (NOTRANS) {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CAPPLY) {
        for (J = 1; J <= N; J++) {
          TMP = TMP + CABS1(A[I][J]) / C[J];
        }
      } else {
        for (J = 1; J <= N; J++) {
          TMP = TMP + CABS1(A[I][J]);
        }
      }
      RWORK[I] = TMP;
      ANORM = max(ANORM, TMP);
    }
  } else {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CAPPLY) {
        for (J = 1; J <= N; J++) {
          TMP = TMP + CABS1(A[J][I]) / C[J];
        }
      } else {
        for (J = 1; J <= N; J++) {
          TMP = TMP + CABS1(A[J][I]);
        }
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

      if (NOTRANS) {
        zgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
      } else {
        zgetrs('Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N,
            INFO);
      }

      // Multiply by inv(C).

      if (CAPPLY) {
        for (I = 1; I <= N; I++) {
          WORK[I] = WORK[I] * C[I].toComplex();
        }
      }
    } else {
      // Multiply by inv(C**H).

      if (CAPPLY) {
        for (I = 1; I <= N; I++) {
          WORK[I] = WORK[I] * C[I].toComplex();
        }
      }

      if (NOTRANS) {
        zgetrs('Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N,
            INFO);
      } else {
        zgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
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
