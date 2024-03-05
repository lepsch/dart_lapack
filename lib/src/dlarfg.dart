import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlarfg(
  final int N,
  final Box<double> ALPHA,
  final Array<double> X_,
  final int INCX,
  final Box<double> TAU,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  const ONE = 1.0, ZERO = 0.0;
  int J, KNT = 0;
  double BETA = 0, RSAFMN = 0, SAFMIN = 0, XNORM;

  if (N <= 1) {
    TAU.value = ZERO;
    return;
  }

  XNORM = dnrm2(N - 1, X, INCX);

  if (XNORM == ZERO) {
    // H  =  I
    TAU.value = ZERO;
    return;
  }

  // general case

  BETA = -sign(dlapy2(ALPHA.value, XNORM), ALPHA.value).toDouble();
  SAFMIN = dlamch('S') / dlamch('E');
  KNT = 0;
  if ((BETA).abs() < SAFMIN) {
    // XNORM, BETA may be inaccurate; scale X and recompute them

    RSAFMN = ONE / SAFMIN;
    do {
      KNT = KNT + 1;
      dscal(N - 1, RSAFMN, X, INCX);
      BETA = BETA * RSAFMN;
      ALPHA.value = ALPHA.value * RSAFMN;
    } while (((BETA).abs() < SAFMIN) && (KNT < 20));

    // New BETA is at most 1, at least SAFMIN

    XNORM = dnrm2(N - 1, X, INCX);
    BETA = -sign(dlapy2(ALPHA.value, XNORM), ALPHA.value).toDouble();
  }
  TAU.value = (BETA - ALPHA.value) / BETA;
  dscal(N - 1, ONE / (ALPHA.value - BETA), X, INCX);

  // If ALPHA.value is subnormal, it may lose relative accuracy

  for (J = 1; J <= KNT; J++) {
    BETA = BETA * SAFMIN;
  }
  ALPHA.value = BETA;
}
