import 'dart:math';

import 'package:lapack/src/blas/dtrmm.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dpotrf.dart';
import 'package:lapack/src/dsyev_2stage.dart';
import 'package:lapack/src/dsygst.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsygv_2stage(
  final int ITYPE,
  final String JOBZ,
  final String UPLO,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> B,
  final int LDB,
  final Array<double> W,
  final Array<double> WORK,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  bool LQUERY, UPPER, WANTZ;
  String TRANS;
  int NEIG, LWMIN = 0, LHTRD, LWTRD, KD, IB;

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);

  INFO.value = 0;
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!(lsame(JOBZ, 'N'))) {
    INFO.value = -2;
  } else if (!(UPPER || lsame(UPLO, 'L'))) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  }

  if (INFO.value == 0) {
    KD = ilaenv2stage(1, 'DSYTRD_2STAGE', JOBZ, N, -1, -1, -1);
    IB = ilaenv2stage(2, 'DSYTRD_2STAGE', JOBZ, N, KD, -1, -1);
    LHTRD = ilaenv2stage(3, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1);
    LWTRD = ilaenv2stage(4, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1);
    LWMIN = 2 * N + LHTRD + LWTRD;
    WORK[1] = LWMIN.toDouble();

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -11;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSYGV_2STAGE ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a Cholesky factorization of B.

  dpotrf(UPLO, N, B, LDB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem and solve.

  dsygst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO);
  dsyev_2stage(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO);

  if (WANTZ) {
    // Backtransform eigenvectors to the original problem.

    NEIG = N;
    if (INFO.value > 0) NEIG = INFO.value - 1;
    if (ITYPE == 1 || ITYPE == 2) {
      // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
      // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

      if (UPPER) {
        TRANS = 'N';
      } else {
        TRANS = 'T';
      }

      dtrsm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA);
    } else if (ITYPE == 3) {
      // For B*A*x=(lambda)*x;
      // backtransform eigenvectors: x = L*y or U**T*y

      if (UPPER) {
        TRANS = 'T';
      } else {
        TRANS = 'N';
      }

      dtrmm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA);
    }
  }

  WORK[1] = LWMIN.toDouble();
}
