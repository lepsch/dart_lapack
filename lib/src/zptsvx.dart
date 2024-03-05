import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlanht.dart';
import 'package:lapack/src/zptcon.dart';
import 'package:lapack/src/zptrfs.dart';
import 'package:lapack/src/zpttrf.dart';
import 'package:lapack/src/zpttrs.dart';

void zptsvx(
  final String FACT,
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<Complex> E_,
  final Array<double> DF_,
  final Array<Complex> EF_,
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
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final E = E_.having();
  final EF = EF_.having();
  final D = D_.having();
  final DF = DF_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  bool NOFACT;
  double ANORM;

  // Test the input parameters.

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  if (!NOFACT && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDX < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('ZPTSVX', -INFO.value);
    return;
  }

  if (NOFACT) {
    // Compute the L*D*L**H (or U**H*D*U) factorization of A.

    dcopy(N, D, 1, DF, 1);
    if (N > 1) zcopy(N - 1, E, 1, EF, 1);
    zpttrf(N, DF, EF, INFO);

    // Return if INFO.value is non-zero.

    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.

  ANORM = zlanht('1', N, D, E);

  // Compute the reciprocal of the condition number of A.

  zptcon(N, DF, EF, ANORM, RCOND, RWORK, INFO);

  // Compute the solution vectors X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zpttrs('Lower', N, NRHS, DF, EF, X, LDX, INFO);

  // Use iterative refinement to improve the computed solutions and
  // compute error bounds and backward error estimates for them.

  zptrfs('Lower', N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK,
      RWORK, INFO);

  // Set INFO.value = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;
}
