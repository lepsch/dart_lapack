import 'dart:math';

import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zlaic1(
  final int JOB,
  final int J,
  final Array<Complex> X_,
  final double SEST,
  final Array<Complex> W_,
  final Complex GAMMA,
  final Box<double> SESTPR,
  final Box<Complex> S,
  final Box<Complex> C,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final W = W_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const HALF = 0.5, FOUR = 4.0;
  double ABSALP,
      ABSEST,
      ABSGAM,
      B,
      EPS,
      NORMA,
      S1,
      S2,
      SCL,
      T,
      TEST,
      TMP,
      ZETA1,
      ZETA2;
  Complex ALPHA, COSINE, SINE;

  EPS = dlamch('Epsilon');
  ALPHA = zdotc(J, X, 1, W, 1);

  ABSALP = (ALPHA).abs();
  ABSGAM = (GAMMA).abs();
  ABSEST = (SEST).abs();

  if (JOB == 1) {
    // Estimating largest singular value

    // special cases

    if (SEST == ZERO) {
      S1 = max(ABSGAM, ABSALP);
      if (S1 == ZERO) {
        S.value = Complex.zero;
        C.value = Complex.one;
        SESTPR.value = ZERO;
      } else {
        S.value = ALPHA / S1.toComplex();
        C.value = GAMMA / S1.toComplex();
        TMP = (S.value * S.value.conjugate() + C.value * C.value.conjugate())
            .sqrt()
            .toDouble();
        S.value = S.value / TMP.toComplex();
        C.value = C.value / TMP.toComplex();
        SESTPR.value = S1 * TMP;
      }
      return;
    } else if (ABSGAM <= EPS * ABSEST) {
      S.value = Complex.one;
      C.value = Complex.zero;
      TMP = max(ABSEST, ABSALP);
      S1 = ABSEST / TMP;
      S2 = ABSALP / TMP;
      SESTPR.value = TMP * sqrt(S1 * S1 + S2 * S2);
      return;
    } else if (ABSALP <= EPS * ABSEST) {
      S1 = ABSGAM;
      S2 = ABSEST;
      if (S1 <= S2) {
        S.value = Complex.one;
        C.value = Complex.zero;
        SESTPR.value = S2;
      } else {
        S.value = Complex.zero;
        C.value = Complex.one;
        SESTPR.value = S1;
      }
      return;
    } else if (ABSEST <= EPS * ABSALP || ABSEST <= EPS * ABSGAM) {
      S1 = ABSGAM;
      S2 = ABSALP;
      if (S1 <= S2) {
        TMP = S1 / S2;
        SCL = sqrt(ONE + TMP * TMP);
        SESTPR.value = S2 * SCL;
        S.value = (ALPHA / S2.toComplex()) / SCL.toComplex();
        C.value = (GAMMA / S2.toComplex()) / SCL.toComplex();
      } else {
        TMP = S2 / S1;
        SCL = sqrt(ONE + TMP * TMP);
        SESTPR.value = S1 * SCL;
        S.value = (ALPHA / S1.toComplex()) / SCL.toComplex();
        C.value = (GAMMA / S1.toComplex()) / SCL.toComplex();
      }
      return;
    } else {
      // normal case

      ZETA1 = ABSALP / ABSEST;
      ZETA2 = ABSGAM / ABSEST;

      B = (ONE - ZETA1 * ZETA1 - ZETA2 * ZETA2) * HALF;
      C.value = (ZETA1 * ZETA1).toComplex();
      if (B > ZERO) {
        T = (C.value / (B.toComplex() + ((B * B).toComplex() + C.value).sqrt()))
            .toDouble();
      } else {
        T = (((B * B).toComplex() + C.value).sqrt() - B.toComplex()).toDouble();
      }

      SINE = -(ALPHA / ABSEST.toComplex()) / T.toComplex();
      COSINE = -(GAMMA / ABSEST.toComplex()) / (ONE + T).toComplex();
      TMP = (SINE * SINE.conjugate() + COSINE * COSINE.conjugate())
          .sqrt()
          .toDouble();

      S.value = SINE / TMP.toComplex();
      C.value = COSINE / TMP.toComplex();
      SESTPR.value = sqrt(T + ONE) * ABSEST;
      return;
    }
  } else if (JOB == 2) {
    // Estimating smallest singular value

    // special cases

    if (SEST == ZERO) {
      SESTPR.value = ZERO;
      if (max(ABSGAM, ABSALP) == ZERO) {
        SINE = Complex.one;
        COSINE = Complex.zero;
      } else {
        SINE = -GAMMA.conjugate();
        COSINE = ALPHA.conjugate();
      }
      S1 = max((SINE).abs(), (COSINE).abs());
      S.value = SINE / S1.toComplex();
      C.value = COSINE / S1.toComplex();
      TMP = (S.value * S.value.conjugate() + C.value * C.value.conjugate())
          .sqrt()
          .toDouble();
      S.value = S.value / TMP.toComplex();
      C.value = C.value / TMP.toComplex();
      return;
    } else if (ABSGAM <= EPS * ABSEST) {
      S.value = Complex.zero;
      C.value = Complex.one;
      SESTPR.value = ABSGAM;
      return;
    } else if (ABSALP <= EPS * ABSEST) {
      S1 = ABSGAM;
      S2 = ABSEST;
      if (S1 <= S2) {
        S.value = Complex.zero;
        C.value = Complex.one;
        SESTPR.value = S1;
      } else {
        S.value = Complex.one;
        C.value = Complex.zero;
        SESTPR.value = S2;
      }
      return;
    } else if (ABSEST <= EPS * ABSALP || ABSEST <= EPS * ABSGAM) {
      S1 = ABSGAM;
      S2 = ABSALP;
      if (S1 <= S2) {
        TMP = S1 / S2;
        SCL = sqrt(ONE + TMP * TMP);
        SESTPR.value = ABSEST * (TMP / SCL);
        S.value =
            -(GAMMA.conjugate().conjugate() / S2.toComplex()) / SCL.toComplex();
        C.value = (ALPHA.conjugate() / S2.toComplex()) / SCL.toComplex();
      } else {
        TMP = S2 / S1;
        SCL = sqrt(ONE + TMP * TMP);
        SESTPR.value = ABSEST / SCL;
        S.value = -(GAMMA.conjugate() / S1.toComplex()) / SCL.toComplex();
        C.value = (ALPHA.conjugate() / S1.toComplex()) / SCL.toComplex();
      }
      return;
    } else {
      // normal case

      ZETA1 = ABSALP / ABSEST;
      ZETA2 = ABSGAM / ABSEST;

      NORMA = max(
          ONE + ZETA1 * ZETA1 + ZETA1 * ZETA2, ZETA1 * ZETA2 + ZETA2 * ZETA2);

      // See if root is closer to zero or to ONE

      TEST = ONE + TWO * (ZETA1 - ZETA2) * (ZETA1 + ZETA2);
      if (TEST >= ZERO) {
        // root is close to zero, compute directly

        B = (ZETA1 * ZETA1 + ZETA2 * ZETA2 + ONE) * HALF;
        C.value = (ZETA2 * ZETA2).toComplex();
        T = (C.value /
                (B.toComplex() +
                    ((B * B).toComplex() - C.value).sqrt().abs().toComplex()))
            .toDouble();
        SINE = (ALPHA / ABSEST.toComplex()) / (ONE - T).toComplex();
        COSINE = -(GAMMA / ABSEST.toComplex()) / T.toComplex();
        SESTPR.value = sqrt(T + FOUR * EPS * EPS * NORMA) * ABSEST;
      } else {
        // root is closer to ONE, shift by that amount

        B = (ZETA2 * ZETA2 + ZETA1 * ZETA1 - ONE) * HALF;
        C.value = (ZETA1 * ZETA1).toComplex();
        if (B >= ZERO) {
          T = (-C.value /
                  (B.toComplex() + ((B * B).toComplex() + C.value).sqrt()))
              .toDouble();
        } else {
          T = (B.toComplex() - ((B * B).toComplex() + C.value).sqrt())
              .toDouble();
        }
        SINE = -(ALPHA / ABSEST.toComplex()) / T.toComplex();
        COSINE = -(GAMMA / ABSEST.toComplex()) / (ONE + T).toComplex();
        SESTPR.value = sqrt(ONE + T + FOUR * EPS * EPS * NORMA) * ABSEST;
      }
      TMP = (SINE * SINE.conjugate() + COSINE * COSINE.conjugate())
          .sqrt()
          .toDouble();
      S.value = SINE / TMP.toComplex();
      C.value = COSINE / TMP.toComplex();
      return;
    }
  }
}
