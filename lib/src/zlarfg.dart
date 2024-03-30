import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlapy3.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zladiv.dart';

void zlarfg(
  final int N,
  final Box<Complex> ALPHA,
  final Array<Complex> X_,
  final int INCX,
  final Box<Complex> TAU,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  const ONE = 1.0, ZERO = 0.0;
  int J, KNT = 0;
  double ALPHI, ALPHR, BETA = 0, RSAFMN, SAFMIN = 0, XNORM;

  if (N <= 0) {
    TAU.value = ZERO.toComplex();
    return;
  }

  XNORM = dznrm2(N - 1, X, INCX);
  ALPHR = ALPHA.value.real;
  ALPHI = ALPHA.value.imaginary;

  if (XNORM == ZERO && ALPHI == ZERO) {
    // H  =  I

    TAU.value = ZERO.toComplex();
  } else {
    // general case

    BETA = -sign(dlapy3(ALPHR, ALPHI, XNORM), ALPHR);
    SAFMIN = dlamch('S') / dlamch('E');
    RSAFMN = ONE / SAFMIN;

    KNT = 0;
    if (BETA.abs() < SAFMIN) {
      // XNORM, BETA may be inaccurate; scale X and recompute them

      do {
        KNT++;
        zdscal(N - 1, RSAFMN, X, INCX);
        BETA *= RSAFMN;
        ALPHI *= RSAFMN;
        ALPHR *= RSAFMN;
      } while ((BETA.abs() < SAFMIN) && (KNT < 20));

      // New BETA is at most 1, at least SAFMIN

      XNORM = dznrm2(N - 1, X, INCX);
      ALPHA.value = Complex(ALPHR, ALPHI);
      BETA = -sign(dlapy3(ALPHR, ALPHI, XNORM), ALPHR);
    }
    TAU.value = Complex((BETA - ALPHR) / BETA, -ALPHI / BETA);
    ALPHA.value = zladiv(Complex(ONE), ALPHA.value - BETA.toComplex());
    zscal(N - 1, ALPHA.value, X, INCX);

    // If ALPHA is subnormal, it may lose relative accuracy

    for (J = 1; J <= KNT; J++) {
      BETA *= SAFMIN;
    }
    ALPHA.value = BETA.toComplex();
  }
}
