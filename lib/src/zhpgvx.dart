import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/ztpmv.dart';
import 'package:lapack/src/blas/ztpsv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhpevx.dart';
import 'package:lapack/src/zhpgst.dart';
import 'package:lapack/src/zpptrf.dart';

void zhpgvx(
  final int ITYPE,
  final String JOBZ,
  final String RANGE,
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Array<Complex> BP_,
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
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Array<int> IFAIL_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Z = Z_.dim(LDZ);
  final W = W_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();
  final IFAIL = IFAIL_.dim();
  final AP = AP_.dim();
  final BP = BP_.dim();
  bool ALLEIG, INDEIG, UPPER, VALEIG, WANTZ;
  String TRANS;
  int J;

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');

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
  } else {
    if (VALEIG) {
      if (N > 0 && VU <= VL) {
        INFO.value = -9;
      }
    } else if (INDEIG) {
      if (IL < 1) {
        INFO.value = -10;
      } else if (IU < min(N, IL) || IU > N) {
        INFO.value = -11;
      }
    }
  }
  if (INFO.value == 0) {
    if (LDZ < 1 || (WANTZ && LDZ < N)) {
      INFO.value = -16;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHPGVX', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a Cholesky factorization of B.

  zpptrf(UPLO, N, BP, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem and solve.

  zhpgst(ITYPE, UPLO, N, AP, BP, INFO);
  zhpevx(JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
      RWORK, IWORK, IFAIL, INFO);

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

      for (J = 1; J <= M.value; J++) {
        // 10
        ztpsv(UPLO, TRANS, 'Non-unit', N, BP, Z(1, J).asArray(), 1);
      } // 10
    } else if (ITYPE == 3) {
      // For B*A*x=(lambda)*x;
      // backtransform eigenvectors: x = L*y or U**H *y

      if (UPPER) {
        TRANS = 'C';
      } else {
        TRANS = 'N';
      }

      for (J = 1; J <= M.value; J++) {
        // 20
        ztpmv(UPLO, TRANS, 'Non-unit', N, BP, Z(1, J).asArray(), 1);
      } // 20
    }
  }
}
