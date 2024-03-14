import 'dart:math';

import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/ztbmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlantb.dart';

void ztbt02(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int KD,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0 or NRHS = 0

  if (N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Compute the 1-norm of op(A).

  final double ANORM;
  if (lsame(TRANS, 'N')) {
    ANORM = zlantb('1', UPLO, DIAG, N, KD, AB, LDAB, RWORK);
  } else {
    ANORM = zlantb('I', UPLO, DIAG, N, KD, AB, LDAB, RWORK);
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute the maximum over the number of right hand sides of
  //    norm(op(A)*x - b) / ( norm(op(A)) * norm(x) * EPS ).

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    zcopy(N, X(1, J).asArray(), 1, WORK, 1);
    ztbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1);
    zaxpy(N, Complex(-ONE), B(1, J).asArray(), 1, WORK, 1);
    final BNORM = dzasum(N, WORK, 1);
    final XNORM = dzasum(N, X(1, J).asArray(), 1);
    if (XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value = max(RESID.value, ((BNORM / ANORM) / XNORM) / EPS);
    }
  }
}
