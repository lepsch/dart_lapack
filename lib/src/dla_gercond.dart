import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgetrs.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

double dla_gercond(
  final String TRANS,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AF_,
  final int LDAF,
  final Array<int> IPIV_,
  final int CMODE,
  final Array<double> C_,
  final Box<int> INFO,
  final Array<double> WORK_,
  final Array<int> IWORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final IPIV = IPIV_.having();
  final C = C_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();

  bool NOTRANS;
  int I, J;
  double TMP;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0);
  final KASE = Box(0);

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
    xerbla('DLA_GERCOND', -INFO.value);
    return 0;
  }
  if (N == 0) {
    return 1.0;
  }

  // Compute the equilibration matrix R such that
  // inv(R)*A*C has unit 1-norm.

  if (NOTRANS) {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CMODE == 1) {
        for (J = 1; J <= N; J++) {
          TMP += (A[I][J] * C[J]).abs();
        }
      } else if (CMODE == 0) {
        for (J = 1; J <= N; J++) {
          TMP += (A[I][J]).abs();
        }
      } else {
        for (J = 1; J <= N; J++) {
          TMP += (A[I][J] / C[J]).abs();
        }
      }
      WORK[2 * N + I] = TMP;
    }
  } else {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CMODE == 1) {
        for (J = 1; J <= N; J++) {
          TMP += (A[J][I] * C[J]).abs();
        }
      } else if (CMODE == 0) {
        for (J = 1; J <= N; J++) {
          TMP += (A[J][I]).abs();
        }
      } else {
        for (J = 1; J <= N; J++) {
          TMP += (A[J][I] / C[J]).abs();
        }
      }
      WORK[2 * N + I] = TMP;
    }
  }

  // Estimate the norm of inv(op(A)).

  AINVNM.value = 0.0;

  KASE.value = 0;
  while (true) {
    dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;

    if (KASE.value == 2) {
      // Multiply by R.

      for (I = 1; I <= N; I++) {
        WORK[I] *= WORK[2 * N + I];
      }

      if (NOTRANS) {
        dgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
      } else {
        dgetrs('Transpose', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
      }

      // Multiply by inv(C).

      if (CMODE == 1) {
        for (I = 1; I <= N; I++) {
          WORK[I] /= C[I];
        }
      } else if (CMODE == -1) {
        for (I = 1; I <= N; I++) {
          WORK[I] *= C[I];
        }
      }
    } else {
      // Multiply by inv(C**T).

      if (CMODE == 1) {
        for (I = 1; I <= N; I++) {
          WORK[I] /= C[I];
        }
      } else if (CMODE == -1) {
        for (I = 1; I <= N; I++) {
          WORK[I] *= C[I];
        }
      }

      if (NOTRANS) {
        dgetrs('Transpose', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
      } else {
        dgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
      }

      // Multiply by R.

      for (I = 1; I <= N; I++) {
        WORK[I] *= WORK[2 * N + I];
      }
    }
  }

  // Compute the estimate of the reciprocal condition number.

  return (AINVNM.value != 0.0) ? (1.0 / AINVNM.value) : 0;
}
