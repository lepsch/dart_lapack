import 'dart:math';

import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgtsv.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsytrs_aa(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final IPIV = IPIV_.dim();
  final B = B_.dim(LDB);
  final WORK = WORK_.dim();
  const ONE = 1.0;
  bool LQUERY, UPPER;
  int K, KP, LWKMIN;

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);
  if (min(N, NRHS) == 0) {
    LWKMIN = 1;
  } else {
    LWKMIN = 3 * N - 2;
  }

  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  } else if (LWORK < LWKMIN && !LQUERY) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('DSYTRS_AA', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = LWKMIN.toDouble();
    return;
  }

  // Quick return if possible

  if (min(N, NRHS) == 0) return;

  if (UPPER) {
    // Solve A*X = B, where A = U**T*T*U.

    // 1) Forward substitution with U**T

    if (N > 1) {
      // Pivot, P**T * B -> B

      for (K = 1; K <= N; K++) {
        KP = IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
      }

      // Compute U**T \ B -> B    [ (U**T \P**T * B) ]

      dtrsm('L', 'U', 'T', 'U', N - 1, NRHS, ONE, A(1, 2), LDA, B(2, 1), LDB);
    }

    // 2) Solve with triangular matrix T

    // Compute T \ B -> B   [ T \ (U**T \P**T * B) ]

    dlacpy('F', 1, N, A(1, 1), LDA + 1, WORK(N).asMatrix(1), 1);
    if (N > 1) {
      dlacpy('F', 1, N - 1, A(1, 2), LDA + 1, WORK(1).asMatrix(1), 1);
      dlacpy('F', 1, N - 1, A(1, 2), LDA + 1, WORK(2 * N).asMatrix(1), 1);
    }
    dgtsv(N, NRHS, WORK(1), WORK(N), WORK(2 * N), B, LDB, INFO);

    // 3) Backward substitution with U

    if (N > 1) {
      // Compute U \ B -> B   [ U \ (T \ (U**T \P**T * B) ) ]

      dtrsm('L', 'U', 'N', 'U', N - 1, NRHS, ONE, A(1, 2), LDA, B(2, 1), LDB);

      // Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]

      for (K = N; K >= 1; K--) {
        KP = IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
      }
    }
  } else {
    // Solve A*X = B, where A = L*T*L**T.

    // 1) Forward substitution with L

    if (N > 1) {
      // Pivot, P**T * B -> B

      for (K = 1; K <= N; K++) {
        KP = IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
      }

      // Compute L \ B -> B    [ (L \P**T * B) ]

      dtrsm('L', 'L', 'N', 'U', N - 1, NRHS, ONE, A(2, 1), LDA, B(2, 1), LDB);
    }

    // 2) Solve with triangular matrix T

    // Compute T \ B -> B   [ T \ (L \P**T * B) ]

    dlacpy('F', 1, N, A(1, 1), LDA + 1, WORK(N).asMatrix(1), 1);
    if (N > 1) {
      dlacpy('F', 1, N - 1, A(2, 1), LDA + 1, WORK(1).asMatrix(1), 1);
      dlacpy('F', 1, N - 1, A(2, 1), LDA + 1, WORK(2 * N).asMatrix(1), 1);
    }
    dgtsv(N, NRHS, WORK(1), WORK(N), WORK(2 * N), B, LDB, INFO);

    // 3) Backward substitution with L**T

    if (N > 1) {
      // Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]

      dtrsm('L', 'L', 'T', 'U', N - 1, NRHS, ONE, A(2, 1), LDA, B(2, 1), LDB);

      // Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]

      for (K = N; K >= 1; K--) {
        KP = IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
      }
    }
  }
}
