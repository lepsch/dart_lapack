import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/dlapy3.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zladiv.dart';

void zlarfgp(
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
  const TWO = 2.0, ONE = 1.0, ZERO = 0.0;
  int J, KNT;
  double ALPHI, ALPHR, BETA, BIGNUM, EPS, SMLNUM, XNORM;
  Complex SAVEALPHA;

  if (N <= 0) {
    TAU.value = Complex.zero;
    return;
  }

  EPS = dlamch('Precision');
  XNORM = dznrm2(N - 1, X, INCX);
  ALPHR = ALPHA.value.real;
  ALPHI = ALPHA.value.imaginary;

  if (XNORM <= EPS * ALPHA.value.abs() && ALPHI == ZERO) {
    // H  =  [1-alpha/abs(alpha) 0; 0 I], sign chosen so ALPHA >= 0.

    if (ALPHR >= ZERO) {
      // When TAU == ZERO, the vector is special-cased to be
      // all zeros in the application routines.  We do not need
      // to clear it.
      TAU.value = Complex.zero;
    } else {
      // However, the application routines rely on explicit
      // zero checks when TAU != ZERO, and we must clear X.
      TAU.value = TWO.toComplex();
      for (J = 1; J <= N - 1; J++) {
        X[1 + (J - 1) * INCX] = Complex.zero;
      }
      ALPHA.value = -ALPHA.value;
    }
  } else {
    // general case

    BETA = sign(dlapy3(ALPHR, ALPHI, XNORM), ALPHR);
    SMLNUM = dlamch('S') / dlamch('E');
    BIGNUM = ONE / SMLNUM;

    KNT = 0;
    if (BETA.abs() < SMLNUM) {
      // XNORM, BETA may be inaccurate; scale X and recompute them

      do {
        KNT++;
        zdscal(N - 1, BIGNUM, X, INCX);
        BETA *= BIGNUM;
        ALPHI *= BIGNUM;
        ALPHR *= BIGNUM;
      } while ((BETA.abs() < SMLNUM) && (KNT < 20));

      // New BETA is at most 1, at least SMLNUM

      XNORM = dznrm2(N - 1, X, INCX);
      ALPHA.value = Complex(ALPHR, ALPHI);
      BETA = sign(dlapy3(ALPHR, ALPHI, XNORM), ALPHR);
    }
    SAVEALPHA = ALPHA.value;
    ALPHA.value += BETA.toComplex();
    if (BETA < ZERO) {
      BETA = -BETA;
      TAU.value = -ALPHA.value / BETA.toComplex();
    } else {
      ALPHR = ALPHI * (ALPHI / ALPHA.value.real);
      ALPHR += XNORM * (XNORM / ALPHA.value.real);
      TAU.value = Complex(ALPHR / BETA, -ALPHI / BETA);
      ALPHA.value = Complex(-ALPHR, ALPHI);
    }
    ALPHA.value = zladiv(Complex.one, ALPHA.value);

    if (TAU.value.abs() <= SMLNUM) {
      // In the case where the computed TAU ends up being a denormalized number,
      // it loses relative accuracy. This is a BIG problem. Solution: flush TAU
      // to ZERO (or TWO or whatever makes a nonnegative real number for BETA).

      // (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
      // (Thanks Pat. Thanks MathWorks.)

      ALPHR = SAVEALPHA.real;
      ALPHI = SAVEALPHA.imaginary;
      if (ALPHI == ZERO) {
        if (ALPHR >= ZERO) {
          TAU.value = Complex.zero;
        } else {
          TAU.value = TWO.toComplex();
          for (J = 1; J <= N - 1; J++) {
            X[1 + (J - 1) * INCX] = Complex.zero;
          }
          BETA = (-SAVEALPHA).real;
        }
      } else {
        XNORM = dlapy2(ALPHR, ALPHI);
        TAU.value = Complex(ONE - ALPHR / XNORM, -ALPHI / XNORM);
        for (J = 1; J <= N - 1; J++) {
          X[1 + (J - 1) * INCX] = Complex.zero;
        }
        BETA = XNORM;
      }
    } else {
      // This is the general case.

      zscal(N - 1, ALPHA.value, X, INCX);
    }

    // If BETA is subnormal, it may lose relative accuracy

    for (J = 1; J <= KNT; J++) {
      BETA *= SMLNUM;
    }
    ALPHA.value = BETA.toComplex();
  }
}
