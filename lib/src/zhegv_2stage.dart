import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/cholesky/top/zpotrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zheev_2stage.dart';
import 'package:lapack/src/zhegst.dart';

void zhegv_2stage(
  final int ITYPE,
  final String JOBZ,
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> W_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final W = W_.dim();
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
    KD = ilaenv2stage(1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1);
    IB = ilaenv2stage(2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1);
    LHTRD = ilaenv2stage(3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1);
    LWTRD = ilaenv2stage(4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1);
    LWMIN = N + LHTRD + LWTRD;
    WORK[1] = LWMIN.toComplex();

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -11;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHEGV_2STAGE ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a Cholesky factorization of B.

  zpotrf(UPLO, N, B, LDB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem and solve.

  zhegst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO);
  zheev_2stage(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO);

  if (WANTZ) {
    // Backtransform eigenvectors to the original problem.

    NEIG = N;
    if (INFO.value > 0) NEIG = INFO.value - 1;
    if (ITYPE == 1 || ITYPE == 2) {
      // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
      // backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y

      if (UPPER) {
        TRANS = 'N';
      } else {
        TRANS = 'C';
      }

      ztrsm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, Complex.one, B, LDB, A,
          LDA);
    } else if (ITYPE == 3) {
      // For B*A*x=(lambda)*x;
      // backtransform eigenvectors: x = L*y or U**H *y

      if (UPPER) {
        TRANS = 'C';
      } else {
        TRANS = 'N';
      }

      ztrmm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, Complex.one, B, LDB, A,
          LDA);
    }
  }

  WORK[1] = LWMIN.toComplex();
}
