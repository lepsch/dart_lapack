import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgtcon.dart';
import 'package:lapack/src/zgtrfs.dart';
import 'package:lapack/src/zgttrf.dart';
import 'package:lapack/src/zgttrs.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlangt.dart';

void zgtsvx(
  final String FACT,
  final String TRANS,
  final int N,
  final int NRHS,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Array<Complex> DLF_,
  final Array<Complex> DF_,
  final Array<Complex> DUF_,
  final Array<Complex> DU2_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> X_,
  final int LDX,
  final Box<double> RCOND,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DLF = DLF_.having();
  final DF = DF_.having();
  final DUF = DUF_.having();
  final DU2 = DU2_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  const ZERO = 0.0;
  bool NOFACT, NOTRAN;
  String NORM;
  double ANORM;

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  NOTRAN = lsame(TRANS, 'N');
  if (!NOFACT && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDB < max(1, N)) {
    INFO.value = -14;
  } else if (LDX < max(1, N)) {
    INFO.value = -16;
  }
  if (INFO.value != 0) {
    xerbla('ZGTSVX', -INFO.value);
    return;
  }

  if (NOFACT) {
    // Compute the LU factorization of A.

    zcopy(N, D, 1, DF, 1);
    if (N > 1) {
      zcopy(N - 1, DL, 1, DLF, 1);
      zcopy(N - 1, DU, 1, DUF, 1);
    }
    zgttrf(N, DLF, DF, DUF, DU2, IPIV, INFO);

    // Return if INFO is non-zero.

    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.

  if (NOTRAN) {
    NORM = '1';
  } else {
    NORM = 'I';
  }
  ANORM = zlangt(NORM, N, DL, D, DU);

  // Compute the reciprocal of the condition number of A.

  zgtcon(NORM, N, DLF, DF, DUF, DU2, IPIV, ANORM, RCOND, WORK, INFO);

  // Compute the solution vectors X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zgttrs(TRANS, N, NRHS, DLF, DF, DUF, DU2, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solutions and
  // compute error bounds and backward error estimates for them.

  zgtrfs(TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX,
      FERR, BERR, WORK, RWORK, INFO);

  // Set INFO = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;
}
