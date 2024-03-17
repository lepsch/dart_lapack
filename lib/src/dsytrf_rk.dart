import 'dart:math';

import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlasyf_rk.dart';
import 'package:lapack/src/dsytf2_rk.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsytrf_rk(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> E_,
  final Array<int> IPIV_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  bool LQUERY, UPPER;
  int I, IP, IWS, K, LDWORK, LWKOPT = 0, NB = 0, NBMIN;
  final IINFO = Box(0), KB = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -8;
  }

  if (INFO.value == 0) {
    // Determine the block size

    NB = ilaenv(1, 'DSYTRF_RK', UPLO, N, -1, -1, -1);
    LWKOPT = max(1, N * NB);
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DSYTRF_RK', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  NBMIN = 2;
  LDWORK = N;
  if (NB > 1 && NB < N) {
    IWS = LDWORK * NB;
    if (LWORK < IWS) {
      NB = max(LWORK ~/ LDWORK, 1);
      NBMIN = max(2, ilaenv(2, 'DSYTRF_RK', UPLO, N, -1, -1, -1));
    }
  } else {
    IWS = 1;
  }
  if (NB < NBMIN) NB = N;

  if (UPPER) {
    // Factorize A as U*D*U**T using the upper triangle of A

    // K is the main loop index, decreasing from N to 1 in steps of
    // KB.value, where KB.value is the number of columns factorized by DLASYF_RK;
    // KB.value is either NB or NB-1, or K for the last block

    K = N;
    while (K >= 1) {
      if (K > NB) {
        // Factorize columns k-kb+1:k of A and use blocked code to
        // update columns 1:k-kb

        dlasyf_rk(UPLO, K, NB, KB, A, LDA, E, IPIV, WORK.asMatrix(LDWORK),
            LDWORK, IINFO);
      } else {
        // Use unblocked code to factorize columns 1:k of A

        dsytf2_rk(UPLO, K, A, LDA, E, IPIV, IINFO);
        KB.value = K;
      }

      // Set INFO.value on the first occurrence of a zero pivot

      if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value;

      // No need to adjust IPIV

      // Apply permutations to the leading panel 1:k-1

      // Read IPIV from the last block factored, i.e.
      // indices  k-kb+1:k and apply row permutations to the
      // last k+1 colunms k+1:N after that block
      // (We can do the simple loop over IPIV with decrement -1,
      // since the ABS value of IPIV( I ) represents the row index
      // of the interchange with row i in both 1x1 and 2x2 pivot cases)

      if (K < N) {
        for (I = K; I >= (K - KB.value + 1); I--) {
          IP = IPIV[I].abs();
          if (IP != I) {
            dswap(
                N - K, A(I, K + 1).asArray(), LDA, A(IP, K + 1).asArray(), LDA);
          }
        }
      }

      // Decrease K and return to the start of the main loop

      K -= KB.value;
    }
  } else {
    // Factorize A as L*D*L**T using the lower triangle of A

    // K is the main loop index, increasing from 1 to N in steps of
    // KB.value, where KB.value is the number of columns factorized by DLASYF_RK;
    // KB.value is either NB or NB-1, or N-K+1 for the last block

    K = 1;
    while (K <= N) {
      if (K <= N - NB) {
        // Factorize columns k:k+kb-1 of A and use blocked code to
        // update columns k+kb:n

        dlasyf_rk(UPLO, N - K + 1, NB, KB, A(K, K), LDA, E(K), IPIV(K),
            WORK.asMatrix(LDWORK), LDWORK, IINFO);
      } else {
        // Use unblocked code to factorize columns k:n of A

        dsytf2_rk(UPLO, N - K + 1, A(K, K), LDA, E(K), IPIV(K), IINFO);
        KB.value = N - K + 1;
      }

      // Set INFO.value on the first occurrence of a zero pivot

      if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value + K - 1;

      // Adjust IPIV

      for (I = K; I <= K + KB.value - 1; I++) {
        if (IPIV[I] > 0) {
          IPIV[I] += K - 1;
        } else {
          IPIV[I] -= K + 1;
        }
      }

      // Apply permutations to the leading panel 1:k-1

      // Read IPIV from the last block factored, i.e.
      // indices  k:k+kb-1 and apply row permutations to the
      // first k-1 colunms 1:k-1 before that block
      // (We can do the simple loop over IPIV with increment 1,
      // since the ABS value of IPIV( I ) represents the row index
      // of the interchange with row i in both 1x1 and 2x2 pivot cases)

      if (K > 1) {
        for (I = K; I <= (K + KB.value - 1); I++) {
          IP = IPIV[I].abs();
          if (IP != I) {
            dswap(K - 1, A(I, 1).asArray(), LDA, A(IP, 1).asArray(), LDA);
          }
        }
      }

      // Increase K and return to the start of the main loop

      K += KB.value;
    }
    // End Lower
  }

  WORK[1] = LWKOPT.toDouble();
}
