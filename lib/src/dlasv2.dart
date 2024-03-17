import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';

void dlasv2(
  final double F,
  final double G,
  final double H,
  final Box<double> SSMIN,
  final Box<double> SSMAX,
  final Box<double> SNR,
  final Box<double> CSR,
  final Box<double> SNL,
  final Box<double> CSL,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  const HALF = 0.5;
  const ONE = 1.0;
  const TWO = 2.0;
  const FOUR = 4.0;
  bool GASMAL, SWAP;
  int PMAX;
  double A,
      CLT = 0,
      CRT = 0,
      D,
      FA,
      FT,
      GA,
      GT,
      HA,
      HT,
      L,
      M,
      MM,
      R,
      S,
      SLT = 0,
      SRT = 0,
      T,
      TEMP,
      TSIGN = 0,
      TT;

  FT = F;
  FA = FT.abs();
  HT = H;
  HA = H.abs();

  // PMAX points to the maximum absolute element of matrix
  //   PMAX = 1 if F largest in absolute values
  //   PMAX = 2 if G largest in absolute values
  //   PMAX = 3 if H largest in absolute values

  PMAX = 1;
  SWAP = (HA > FA);
  if (SWAP) {
    PMAX = 3;
    TEMP = FT;
    FT = HT;
    HT = TEMP;
    TEMP = FA;
    FA = HA;
    HA = TEMP;

    // Now FA >= HA
  }
  GT = G;
  GA = GT.abs();
  if (GA == ZERO) {
    // Diagonal matrix

    SSMIN.value = HA;
    SSMAX.value = FA;
    CLT = ONE;
    CRT = ONE;
    SLT = ZERO;
    SRT = ZERO;
  } else {
    GASMAL = true;
    if (GA > FA) {
      PMAX = 2;
      if ((FA / GA) < dlamch('EPS')) {
        // Case of very large GA

        GASMAL = false;
        SSMAX.value = GA;
        if (HA > ONE) {
          SSMIN.value = FA / (GA / HA);
        } else {
          SSMIN.value = (FA / GA) * HA;
        }
        CLT = ONE;
        SLT = HT / GT;
        SRT = ONE;
        CRT = FT / GT;
      }
    }
    if (GASMAL) {
      // Normal case

      D = FA - HA;
      if (D == FA) {
        // Copes with infinite F or H

        L = ONE;
      } else {
        L = D / FA;
      }

      // Note that 0 <= L <= 1

      M = GT / FT;

      // Note that abs(M) <= 1/macheps

      T = TWO - L;

      // Note that T >= 1

      MM = M * M;
      TT = T * T;
      S = sqrt(TT + MM);

      // Note that 1 <= S <= 1 + 1/macheps

      if (L == ZERO) {
        R = M.abs();
      } else {
        R = sqrt(L * L + MM);
      }

      // Note that 0 <= R <= 1 + 1/macheps

      A = HALF * (S + R);

      // Note that 1 <= A <= 1 + abs(M)

      SSMIN.value = HA / A;
      SSMAX.value = FA * A;
      if (MM == ZERO) {
        // Note that M is very tiny

        if (L == ZERO) {
          T = sign(TWO, FT) * sign(ONE, GT).toDouble();
        } else {
          T = GT / sign(D, FT) + M / T;
        }
      } else {
        T = (M / (S + T) + M / (R + L)) * (ONE + A);
      }
      L = sqrt(T * T + FOUR);
      CRT = TWO / L;
      SRT = T / L;
      CLT = (CRT + SRT * M) / A;
      SLT = (HT / FT) * SRT / A;
    }
  }
  if (SWAP) {
    CSL.value = SRT;
    SNL.value = CRT;
    CSR.value = SLT;
    SNR.value = CLT;
  } else {
    CSL.value = CLT;
    SNL.value = SLT;
    CSR.value = CRT;
    SNR.value = SRT;
  }

  // Correct signs of SSMAX and SSMIN

  if (PMAX == 1) {
    TSIGN =
        sign(ONE, CSR.value).toDouble() * sign(ONE, CSL.value) * sign(ONE, F);
  }
  if (PMAX == 2) {
    TSIGN =
        sign(ONE, SNR.value).toDouble() * sign(ONE, CSL.value) * sign(ONE, G);
  }
  if (PMAX == 3) {
    TSIGN =
        sign(ONE, SNR.value).toDouble() * sign(ONE, SNL.value) * sign(ONE, H);
  }
  SSMAX.value = sign(SSMAX.value, TSIGN).toDouble();
  SSMIN.value =
      sign(SSMIN.value, TSIGN * sign(ONE, F) * sign(ONE, H)).toDouble();
}
