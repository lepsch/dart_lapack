import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/ztbmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void ztbt03(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int KD,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final double SCALE,
  final Array<double> CNORM_,
  final double TSCAL,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final CNORM = CNORM_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Quick exit if N = 0

  if (N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }
  final EPS = dlamch('Epsilon');
  final SMLNUM = dlamch('Safe minimum');

  // Compute the norm of the triangular matrix A using the column
  // norms already computed by ZLATBS.

  var TNORM = ZERO;
  if (lsame(DIAG, 'N')) {
    if (lsame(UPLO, 'U')) {
      for (var J = 1; J <= N; J++) {
        TNORM = max(TNORM, TSCAL * AB[KD + 1][J].abs() + CNORM[J]);
      }
    } else {
      for (var J = 1; J <= N; J++) {
        TNORM = max(TNORM, TSCAL * AB[1][J].abs() + CNORM[J]);
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      TNORM = max(TNORM, TSCAL + CNORM[J]);
    }
  }

  // Compute the maximum over the number of right hand sides of
  //    norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    zcopy(N, X(1, J).asArray(), 1, WORK, 1);
    var IX = izamax(N, WORK, 1);
    var XNORM = max(ONE, X[IX][J].abs());
    final XSCAL = (ONE / XNORM) / (KD + 1);
    zdscal(N, XSCAL, WORK, 1);
    ztbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1);
    zaxpy(N, Complex(-SCALE * XSCAL), B(1, J).asArray(), 1, WORK, 1);
    IX = izamax(N, WORK, 1);
    var ERR = TSCAL * WORK[IX].abs();
    IX = izamax(N, X(1, J).asArray(), 1);
    XNORM = X[IX][J].abs();
    if (ERR * SMLNUM <= XNORM) {
      if (XNORM > ZERO) ERR /= XNORM;
    } else {
      if (ERR > ZERO) ERR = ONE / EPS;
    }
    if (ERR * SMLNUM <= TNORM) {
      if (TNORM > ZERO) ERR /= TNORM;
    } else {
      if (ERR > ZERO) ERR = ONE / EPS;
    }
    RESID.value = max(RESID.value, ERR);
  }
}
