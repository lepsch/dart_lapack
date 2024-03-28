import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

double zla_herpvgrw(
  final String UPLO,
  final int N,
  final int INFO,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
  final Array<int> IPIV_,
  final Array<double> WORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  int NCOLS, I, J, K, KP;
  double AMAX, UMAX, RPVGRW, TMP;
  bool UPPER;

  UPPER = lsame('Upper', UPLO);
  if (INFO == 0) {
    if (UPPER) {
      NCOLS = 1;
    } else {
      NCOLS = N;
    }
  } else {
    NCOLS = INFO;
  }

  RPVGRW = 1.0;
  for (I = 1; I <= 2 * N; I++) {
    WORK[I] = 0.0;
  }

  // Find the max magnitude entry of each column of A.  Compute the max
  // for all N columns so we can apply the pivot permutation while
  // looping below.  Assume a full factorization is the common case.

  if (UPPER) {
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= J; I++) {
        WORK[N + I] = max(A[I][J].cabs1(), WORK[N + I]);
        WORK[N + J] = max(A[I][J].cabs1(), WORK[N + J]);
      }
    }
  } else {
    for (J = 1; J <= N; J++) {
      for (I = J; I <= N; I++) {
        WORK[N + I] = max(A[I][J].cabs1(), WORK[N + I]);
        WORK[N + J] = max(A[I][J].cabs1(), WORK[N + J]);
      }
    }
  }

  // Now find the max magnitude entry of each column of U or L.  Also
  // permute the magnitudes of A above so they're in the same order as
  // the factor.

  // The iteration orders and permutations were copied from zsytrs.
  // Calls to SSWAP would be severe overkill.

  if (UPPER) {
    K = N;
    while (K < NCOLS && K > 0) {
      if (IPIV[K] > 0) {
        // 1x1 pivot
        KP = IPIV[K];
        if (KP != K) {
          TMP = WORK[N + K];
          WORK[N + K] = WORK[N + KP];
          WORK[N + KP] = TMP;
        }
        for (I = 1; I <= K; I++) {
          WORK[K] = max(AF[I][K].cabs1(), WORK[K]);
        }
        K--;
      } else {
        // 2x2 pivot
        KP = -IPIV[K];
        TMP = WORK[N + K - 1];
        WORK[N + K - 1] = WORK[N + KP];
        WORK[N + KP] = TMP;
        for (I = 1; I <= K - 1; I++) {
          WORK[K] = max(AF[I][K].cabs1(), WORK[K]);
          WORK[K - 1] = max(AF[I][K - 1].cabs1(), WORK[K - 1]);
        }
        WORK[K] = max(AF[K][K].cabs1(), WORK[K]);
        K -= 2;
      }
    }
    K = NCOLS;
    while (K <= N) {
      if (IPIV[K] > 0) {
        KP = IPIV[K];
        if (KP != K) {
          TMP = WORK[N + K];
          WORK[N + K] = WORK[N + KP];
          WORK[N + KP] = TMP;
        }
        K++;
      } else {
        KP = -IPIV[K];
        TMP = WORK[N + K];
        WORK[N + K] = WORK[N + KP];
        WORK[N + KP] = TMP;
        K += 2;
      }
    }
  } else {
    K = 1;
    while (K <= NCOLS) {
      if (IPIV[K] > 0) {
        // 1x1 pivot
        KP = IPIV[K];
        if (KP != K) {
          TMP = WORK[N + K];
          WORK[N + K] = WORK[N + KP];
          WORK[N + KP] = TMP;
        }
        for (I = K; I <= N; I++) {
          WORK[K] = max(AF[I][K].cabs1(), WORK[K]);
        }
        K++;
      } else {
        // 2x2 pivot
        KP = -IPIV[K];
        TMP = WORK[N + K + 1];
        WORK[N + K + 1] = WORK[N + KP];
        WORK[N + KP] = TMP;
        for (I = K + 1; I <= N; I++) {
          WORK[K] = max(AF[I][K].cabs1(), WORK[K]);
          WORK[K + 1] = max(AF[I][K + 1].cabs1(), WORK[K + 1]);
        }
        WORK[K] = max(AF[K][K].cabs1(), WORK[K]);
        K += 2;
      }
    }
    K = NCOLS;
    while (K >= 1) {
      if (IPIV[K] > 0) {
        KP = IPIV[K];
        if (KP != K) {
          TMP = WORK[N + K];
          WORK[N + K] = WORK[N + KP];
          WORK[N + KP] = TMP;
        }
        K--;
      } else {
        KP = -IPIV[K];
        TMP = WORK[N + K];
        WORK[N + K] = WORK[N + KP];
        WORK[N + KP] = TMP;
        K -= 2;
      }
    }
  }

  // Compute the *inverse* of the max element growth factor.  Dividing
  // by zero would imply the largest entry of the factor's column is
  // zero.  Than can happen when either the column of A is zero or
  // massive pivots made the factor underflow to zero.  Neither counts
  // as growth in itself, so simply ignore terms with zero
  // denominators.

  if (UPPER) {
    for (I = NCOLS; I <= N; I++) {
      UMAX = WORK[I];
      AMAX = WORK[N + I];
      if (UMAX != 0.0) {
        RPVGRW = min(AMAX / UMAX, RPVGRW);
      }
    }
  } else {
    for (I = 1; I <= NCOLS; I++) {
      UMAX = WORK[I];
      AMAX = WORK[N + I];
      if (UMAX != 0.0) {
        RPVGRW = min(AMAX / UMAX, RPVGRW);
      }
    }
  }

  return RPVGRW;
}
