import 'package:lapack/src/blas/dtpmv.dart';
import 'package:lapack/src/blas/dtpsv.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dpptrf.dart';
import 'package:lapack/src/dspev.dart';
import 'package:lapack/src/dspgst.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dspgv(
  final int ITYPE,
  final String JOBZ,
  final String UPLO,
  final int N,
  final Array<double> AP,
  final Array<double> BP,
  final Array<double> W,
  final Matrix<double> Z,
  final int LDZ,
  final Array<double> WORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool UPPER, WANTZ;
  String TRANS;
  int J, NEIG;

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');

  INFO.value = 0;
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -2;
  } else if (!(UPPER || lsame(UPLO, 'L'))) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DSPGV ', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a Cholesky factorization of B.

  dpptrf(UPLO, N, BP, INFO.value);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem and solve.

  dspgst(ITYPE, UPLO, N, AP, BP, INFO.value);
  dspev(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO.value);

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

      for (J = 1; J <= NEIG; J++) {
        dtpsv(UPLO, TRANS, 'Non-unit', N, BP, Z(1, J).asArray(), 1);
      }
    } else if (ITYPE == 3) {
      // For B*A*x=(lambda)*x;
      // backtransform eigenvectors: x = L*y or U**T*y

      if (UPPER) {
        TRANS = 'T';
      } else {
        TRANS = 'N';
      }

      for (J = 1; J <= NEIG; J++) {
        dtpmv(UPLO, TRANS, 'Non-unit', N, BP, Z(1, J).asArray(), 1);
      }
    }
  }
}
