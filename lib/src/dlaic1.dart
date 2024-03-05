import 'dart:math';

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';

void dlaic1(
  final int JOB,
  final int J,
  final Array<double> X_,
  final double SEST,
  final Array<double> W_,
  final double GAMMA,
  final Box<double> SESTPR,
  final Box<double> S,
  final Box<double> C,
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
      ALPHA,
      B,
      COSINE,
      EPS,
      NORMA,
      S1,
      S2,
      SINE,
      T,
      TEST,
      TMP,
      ZETA1,
      ZETA2;

  EPS = dlamch('Epsilon');
  ALPHA = ddot(J, X, 1, W, 1);

  ABSALP = (ALPHA).abs();
  ABSGAM = (GAMMA).abs();
  ABSEST = (SEST).abs();

  if (JOB == 1) {
    // Estimating largest singular value

    // special cases

    if (SEST == ZERO) {
      S1 = max(ABSGAM, ABSALP);
      if (S1 == ZERO) {
        S.value = ZERO;
        C.value = ONE;
        SESTPR.value = ZERO;
      } else {
        S.value = ALPHA / S1;
        C.value = GAMMA / S1;
        TMP = sqrt(S.value * S.value + C.value * C.value);
        S.value = S.value / TMP;
        C.value = C.value / TMP;
        SESTPR.value = S1 * TMP;
      }
      return;
    } else if (ABSGAM <= EPS * ABSEST) {
      S.value = ONE;
      C.value = ZERO;
      TMP = max(ABSEST, ABSALP);
      S1 = ABSEST / TMP;
      S2 = ABSALP / TMP;
      SESTPR.value = TMP * sqrt(S1 * S1 + S2 * S2);
      return;
    } else if (ABSALP <= EPS * ABSEST) {
      S1 = ABSGAM;
      S2 = ABSEST;
      if (S1 <= S2) {
        S.value = ONE;
        C.value = ZERO;
        SESTPR.value = S2;
      } else {
        S.value = ZERO;
        C.value = ONE;
        SESTPR.value = S1;
      }
      return;
    } else if (ABSEST <= EPS * ABSALP || ABSEST <= EPS * ABSGAM) {
      S1 = ABSGAM;
      S2 = ABSALP;
      if (S1 <= S2) {
        TMP = S1 / S2;
        S.value = sqrt(ONE + TMP * TMP);
        SESTPR.value = S2 * S.value;
        C.value = (GAMMA / S2) / S.value;
        S.value = sign(ONE, ALPHA) / S.value;
      } else {
        TMP = S2 / S1;
        C.value = sqrt(ONE + TMP * TMP);
        SESTPR.value = S1 * C.value;
        S.value = (ALPHA / S1) / C.value;
        C.value = sign(ONE, GAMMA) / C.value;
      }
      return;
    } else {
      // normal case

      ZETA1 = ALPHA / ABSEST;
      ZETA2 = GAMMA / ABSEST;

      B = (ONE - ZETA1 * ZETA1 - ZETA2 * ZETA2) * HALF;
      C.value = ZETA1 * ZETA1;
      if (B > ZERO) {
        T = C.value / (B + sqrt(B * B + C.value));
      } else {
        T = sqrt(B * B + C.value) - B;
      }

      SINE = -ZETA1 / T;
      COSINE = -ZETA2 / (ONE + T);
      TMP = sqrt(SINE * SINE + COSINE * COSINE);
      S.value = SINE / TMP;
      C.value = COSINE / TMP;
      SESTPR.value = sqrt(T + ONE) * ABSEST;
      return;
    }
  } else if (JOB == 2) {
    // Estimating smallest singular value

    // special cases

    if (SEST == ZERO) {
      SESTPR.value = ZERO;
      if (max(ABSGAM, ABSALP) == ZERO) {
        SINE = ONE;
        COSINE = ZERO;
      } else {
        SINE = -GAMMA;
        COSINE = ALPHA;
      }
      S1 = max((SINE).abs(), (COSINE).abs());
      S.value = SINE / S1;
      C.value = COSINE / S1;
      TMP = sqrt(S.value * S.value + C.value * C.value);
      S.value = S.value / TMP;
      C.value = C.value / TMP;
      return;
    } else if (ABSGAM <= EPS * ABSEST) {
      S.value = ZERO;
      C.value = ONE;
      SESTPR.value = ABSGAM;
      return;
    } else if (ABSALP <= EPS * ABSEST) {
      S1 = ABSGAM;
      S2 = ABSEST;
      if (S1 <= S2) {
        S.value = ZERO;
        C.value = ONE;
        SESTPR.value = S1;
      } else {
        S.value = ONE;
        C.value = ZERO;
        SESTPR.value = S2;
      }
      return;
    } else if (ABSEST <= EPS * ABSALP || ABSEST <= EPS * ABSGAM) {
      S1 = ABSGAM;
      S2 = ABSALP;
      if (S1 <= S2) {
        TMP = S1 / S2;
        C.value = sqrt(ONE + TMP * TMP);
        SESTPR.value = ABSEST * (TMP / C.value);
        S.value = -(GAMMA / S2) / C.value;
        C.value = sign(ONE, ALPHA) / C.value;
      } else {
        TMP = S2 / S1;
        S.value = sqrt(ONE + TMP * TMP);
        SESTPR.value = ABSEST / S.value;
        C.value = (ALPHA / S1) / S.value;
        S.value = -sign(ONE, GAMMA) / S.value;
      }
      return;
    } else {
      // normal case

      ZETA1 = ALPHA / ABSEST;
      ZETA2 = GAMMA / ABSEST;

      NORMA = max(ONE + ZETA1 * ZETA1 + (ZETA1 * ZETA2).abs(),
          (ZETA1 * ZETA2).abs() + ZETA2 * ZETA2);

      // See if root is closer to zero or to ONE

      TEST = ONE + TWO * (ZETA1 - ZETA2) * (ZETA1 + ZETA2);
      if (TEST >= ZERO) {
        // root is close to zero, compute directly

        B = (ZETA1 * ZETA1 + ZETA2 * ZETA2 + ONE) * HALF;
        C.value = ZETA2 * ZETA2;
        T = C.value / (B + sqrt((B * B - C.value).abs()));
        SINE = ZETA1 / (ONE - T);
        COSINE = -ZETA2 / T;
        SESTPR.value = sqrt(T + FOUR * EPS * EPS * NORMA) * ABSEST;
      } else {
        // root is closer to ONE, shift by that amount

        B = (ZETA2 * ZETA2 + ZETA1 * ZETA1 - ONE) * HALF;
        C.value = ZETA1 * ZETA1;
        if (B >= ZERO) {
          T = -C.value / (B + sqrt(B * B + C.value));
        } else {
          T = B - sqrt(B * B + C.value);
        }
        SINE = -ZETA1 / T;
        COSINE = -ZETA2 / (ONE + T);
        SESTPR.value = sqrt(ONE + T + FOUR * EPS * EPS * NORMA) * ABSEST;
      }
      TMP = sqrt(SINE * SINE + COSINE * COSINE);
      S.value = SINE / TMP;
      C.value = COSINE / TMP;
      return;
    }
  }
}
