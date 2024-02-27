import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/cholesky/top/zpotrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zheevx.dart';
import 'package:lapack/src/zhegst.dart';

void zhegvx(
  final int ITYPE,
  final String JOBZ,
  final String RANGE,
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final double ABSTOL,
  final Box<int> M,
  final Array<double> W_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Array<int> IFAIL_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final Z = Z_.dim(LDZ);
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();
  final IFAIL = IFAIL_.dim();
  final W = W_.dim();
  bool ALLEIG, INDEIG, LQUERY, UPPER, VALEIG, WANTZ;
  String TRANS;
  int LWKOPT = 0, NB;

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');
  LQUERY = (LWORK == -1);

  INFO.value = 0;
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -2;
  } else if (!(ALLEIG || VALEIG || INDEIG)) {
    INFO.value = -3;
  } else if (!(UPPER || lsame(UPLO, 'L'))) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else {
    if (VALEIG) {
      if (N > 0 && VU <= VL) INFO.value = -11;
    } else if (INDEIG) {
      if (IL < 1 || IL > max(1, N)) {
        INFO.value = -12;
      } else if (IU < min(N, IL) || IU > N) {
        INFO.value = -13;
      }
    }
  }
  if (INFO.value == 0) {
    if (LDZ < 1 || (WANTZ && LDZ < N)) {
      INFO.value = -18;
    }
  }

  if (INFO.value == 0) {
    NB = ilaenv(1, 'ZHETRD', UPLO, N, -1, -1, -1);
    LWKOPT = max(1, (NB + 1) * N);
    WORK[1] = LWKOPT.toComplex();

    if (LWORK < max(1, 2 * N) && !LQUERY) {
      INFO.value = -20;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHEGVX', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  M.value = 0;
  if (N == 0) {
    return;
  }

  // Form a Cholesky factorization of B.

  zpotrf(UPLO, N, B, LDB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem and solve.

  zhegst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO);
  zheevx(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ,
      WORK, LWORK, RWORK, IWORK, IFAIL, INFO);

  if (WANTZ) {
    // Backtransform eigenvectors to the original problem.

    if (INFO.value > 0) M.value = INFO.value - 1;
    if (ITYPE == 1 || ITYPE == 2) {
      // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
      // backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y

      if (UPPER) {
        TRANS = 'N';
      } else {
        TRANS = 'C';
      }

      ztrsm('Left', UPLO, TRANS, 'Non-unit', N, M.value, Complex.one, B, LDB, Z,
          LDZ);
    } else if (ITYPE == 3) {
      // For B*A*x=(lambda)*x;
      // backtransform eigenvectors: x = L*y or U**H *y

      if (UPPER) {
        TRANS = 'C';
      } else {
        TRANS = 'N';
      }

      ztrmm('Left', UPLO, TRANS, 'Non-unit', N, M.value, Complex.one, B, LDB, Z,
          LDZ);
    }
  }

  // Set WORK(1) to optimal complex workspace size.

  WORK[1] = LWKOPT.toComplex();
}
