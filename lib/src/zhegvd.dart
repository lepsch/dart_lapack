import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zpotrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zheevd.dart';
import 'package:lapack/src/zhegst.dart';

void zhegvd(
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
  final int LRWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final W = W_.having();
  bool LQUERY, UPPER, WANTZ;
  String TRANS;
  int LIOPT, LIWMIN, LOPT, LROPT, LRWMIN, LWMIN;

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1 || LRWORK == -1 || LIWORK == -1);

  INFO.value = 0;
  if (N <= 1) {
    LWMIN = 1;
    LRWMIN = 1;
    LIWMIN = 1;
  } else if (WANTZ) {
    LWMIN = 2 * N + N * N;
    LRWMIN = 1 + 5 * N + 2 * N * N;
    LIWMIN = 3 + 5 * N;
  } else {
    LWMIN = N + 1;
    LRWMIN = N;
    LIWMIN = 1;
  }
  LOPT = LWMIN;
  LROPT = LRWMIN;
  LIOPT = LIWMIN;
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!(WANTZ || lsame(JOBZ, 'N'))) {
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
    WORK[1] = LOPT.toComplex();
    RWORK[1] = LROPT.toDouble();
    IWORK[1] = LIOPT;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -11;
    } else if (LRWORK < LRWMIN && !LQUERY) {
      INFO.value = -13;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -15;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHEGVD', -INFO.value);
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
  zheevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK,
      INFO);
  LOPT = max(LOPT, WORK[1].toDouble()).toInt();
  LROPT = max(LROPT, RWORK[1]).toInt();
  LIOPT = max(LIOPT, IWORK[1]).toInt();

  if (WANTZ && INFO.value == 0) {
    // Backtransform eigenvectors to the original problem.

    if (ITYPE == 1 || ITYPE == 2) {
      // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
      // backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y

      if (UPPER) {
        TRANS = 'N';
      } else {
        TRANS = 'C';
      }

      ztrsm('Left', UPLO, TRANS, 'Non-unit', N, N, Complex.one, B, LDB, A, LDA);
    } else if (ITYPE == 3) {
      // For B*A*x=(lambda)*x;
      // backtransform eigenvectors: x = L*y or U**H *y

      if (UPPER) {
        TRANS = 'C';
      } else {
        TRANS = 'N';
      }

      ztrmm('Left', UPLO, TRANS, 'Non-unit', N, N, Complex.one, B, LDB, A, LDA);
    }
  }

  WORK[1] = LOPT.toComplex();
  RWORK[1] = LROPT.toDouble();
  IWORK[1] = LIOPT;
}
