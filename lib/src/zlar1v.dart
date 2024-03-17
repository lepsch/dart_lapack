import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zlar1v(
  final int N,
  final int B1,
  final int BN,
  final double LAMBDA,
  final Array<double> D_,
  final Array<double> L_,
  final Array<double> LD_,
  final Array<double> LLD_,
  final double PIVMIN,
  final double GAPTOL,
  final Array<Complex> Z_,
  final bool WANTNC,
  final Box<int> NEGCNT,
  final Box<double> ZTZ,
  final Box<double> MINGMA,
  final Box<int> R,
  final Array<int> ISUPPZ_,
  final Box<double> NRMINV,
  final Box<double> RESID,
  final Box<double> RQCORR,
  final Array<double> WORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Z = Z_.having();
  final D = D_.having();
  final L = L_.having();
  final LD = LD_.having();
  final LLD = LLD_.having();
  final ISUPPZ = ISUPPZ_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool SAWNAN1, SAWNAN2;
  int I, INDLPL, INDP, INDS, INDUMN, NEG1, NEG2, R1, R2;
  double DMINUS, DPLUS, EPS, S, TMP;

  EPS = dlamch('Precision');

  if (R.value == 0) {
    R1 = B1;
    R2 = BN;
  } else {
    R1 = R.value;
    R2 = R.value;
  }

  // Storage for LPLUS
  INDLPL = 0;
  // Storage for UMINUS
  INDUMN = N;
  INDS = 2 * N + 1;
  INDP = 3 * N + 1;

  if (B1 == 1) {
    WORK[INDS] = ZERO;
  } else {
    WORK[INDS + B1 - 1] = LLD[B1 - 1];
  }

  // Compute the stationary transform (using the differential form)
  // until the index R2.

  SAWNAN1 = false;
  NEG1 = 0;
  S = WORK[INDS + B1 - 1] - LAMBDA;
  for (I = B1; I <= R1 - 1; I++) {
    DPLUS = D[I] + S;
    WORK[INDLPL + I] = LD[I] / DPLUS;
    if (DPLUS < ZERO) NEG1++;
    WORK[INDS + I] = S * WORK[INDLPL + I] * L[I];
    S = WORK[INDS + I] - LAMBDA;
  }
  SAWNAN1 = disnan(S);
  if (!SAWNAN1) {
    for (I = R1; I <= R2 - 1; I++) {
      DPLUS = D[I] + S;
      WORK[INDLPL + I] = LD[I] / DPLUS;
      WORK[INDS + I] = S * WORK[INDLPL + I] * L[I];
      S = WORK[INDS + I] - LAMBDA;
    }
    SAWNAN1 = disnan(S);
  }

  if (SAWNAN1) {
    // Runs a slower version of the above loop if a NaN is detected
    NEG1 = 0;
    S = WORK[INDS + B1 - 1] - LAMBDA;
    for (I = B1; I <= R1 - 1; I++) {
      DPLUS = D[I] + S;
      if (DPLUS.abs() < PIVMIN) DPLUS = -PIVMIN;
      WORK[INDLPL + I] = LD[I] / DPLUS;
      if (DPLUS < ZERO) NEG1++;
      WORK[INDS + I] = S * WORK[INDLPL + I] * L[I];
      if (WORK[INDLPL + I] == ZERO) WORK[INDS + I] = LLD[I];
      S = WORK[INDS + I] - LAMBDA;
    }
    for (I = R1; I <= R2 - 1; I++) {
      DPLUS = D[I] + S;
      if (DPLUS.abs() < PIVMIN) DPLUS = -PIVMIN;
      WORK[INDLPL + I] = LD[I] / DPLUS;
      WORK[INDS + I] = S * WORK[INDLPL + I] * L[I];
      if (WORK[INDLPL + I] == ZERO) WORK[INDS + I] = LLD[I];
      S = WORK[INDS + I] - LAMBDA;
    }
  }

  // Compute the progressive transform (using the differential form)
  // until the index R1

  SAWNAN2 = false;
  NEG2 = 0;
  WORK[INDP + BN - 1] = D[BN] - LAMBDA;
  for (I = BN - 1; I >= R1; I--) {
    DMINUS = LLD[I] + WORK[INDP + I];
    TMP = D[I] / DMINUS;
    if (DMINUS < ZERO) NEG2++;
    WORK[INDUMN + I] = L[I] * TMP;
    WORK[INDP + I - 1] = WORK[INDP + I] * TMP - LAMBDA;
  }
  TMP = WORK[INDP + R1 - 1];
  SAWNAN2 = disnan(TMP);

  if (SAWNAN2) {
    // Runs a slower version of the above loop if a NaN is detected
    NEG2 = 0;
    for (I = BN - 1; I >= R1; I--) {
      DMINUS = LLD[I] + WORK[INDP + I];
      if (DMINUS.abs() < PIVMIN) DMINUS = -PIVMIN;
      TMP = D[I] / DMINUS;
      if (DMINUS < ZERO) NEG2++;
      WORK[INDUMN + I] = L[I] * TMP;
      WORK[INDP + I - 1] = WORK[INDP + I] * TMP - LAMBDA;
      if (TMP == ZERO) WORK[INDP + I - 1] = D[I] - LAMBDA;
    }
  }

  // Find the index (from R1 to R2) of the largest (in magnitude)
  // diagonal element of the inverse

  MINGMA.value = WORK[INDS + R1 - 1] + WORK[INDP + R1 - 1];
  if (MINGMA.value < ZERO) NEG1++;
  if (WANTNC) {
    NEGCNT.value = NEG1 + NEG2;
  } else {
    NEGCNT.value = -1;
  }
  if ((MINGMA.value).abs() == ZERO) MINGMA.value = EPS * WORK[INDS + R1 - 1];
  R.value = R1;
  for (I = R1; I <= R2 - 1; I++) {
    TMP = WORK[INDS + I] + WORK[INDP + I];
    if (TMP == ZERO) TMP = EPS * WORK[INDS + I];
    if (TMP.abs() <= (MINGMA.value).abs()) {
      MINGMA.value = TMP;
      R.value = I + 1;
    }
  }

  // Compute the FP vector: solve N^T v = e_r

  ISUPPZ[1] = B1;
  ISUPPZ[2] = BN;
  Z[R.value] = Complex.one;
  ZTZ.value = ONE;

  // Compute the FP vector upwards from R.value

  if (!SAWNAN1 && !SAWNAN2) {
    for (I = R.value - 1; I >= B1; I--) {
      Z[I] = -(WORK[INDLPL + I].toComplex() * Z[I + 1]);
      if ((Z[I].abs() + Z[I + 1].abs()) * LD[I].abs() < GAPTOL) {
        Z[I] = Complex.zero;
        ISUPPZ[1] = I + 1;
        break;
      }
      ZTZ.value += (Z[I] * Z[I]).toDouble();
    }
  } else {
    // Run slower loop if NaN occurred.
    for (I = R.value - 1; I >= B1; I--) {
      if (Z[I + 1] == Complex.zero) {
        Z[I] = -(LD[I + 1] / LD[I]).toComplex() * Z[I + 2];
      } else {
        Z[I] = -(WORK[INDLPL + I].toComplex() * Z[I + 1]);
      }
      if ((Z[I].abs() + Z[I + 1].abs()) * LD[I].abs() < GAPTOL) {
        Z[I] = Complex.zero;
        ISUPPZ[1] = I + 1;
        break;
      }
      ZTZ.value += (Z[I] * Z[I]).toDouble();
    }
  }

  // Compute the FP vector downwards from R.value in blocks of size BLKSIZ
  if (!SAWNAN1 && !SAWNAN2) {
    for (I = R.value; I <= BN - 1; I++) {
      Z[I + 1] = -(WORK[INDUMN + I].toComplex() * Z[I]);
      if ((Z[I].abs() + Z[I + 1].abs()) * LD[I].abs() < GAPTOL) {
        Z[I + 1] = Complex.zero;
        ISUPPZ[2] = I;
        break;
      }
      ZTZ.value += (Z[I + 1] * Z[I + 1]).toDouble();
    }
  } else {
    // Run slower loop if NaN occurred.
    for (I = R.value; I <= BN - 1; I++) {
      if (Z[I] == Complex.zero) {
        Z[I + 1] = -(LD[I - 1] / LD[I]).toComplex() * Z[I - 1];
      } else {
        Z[I + 1] = -(WORK[INDUMN + I].toComplex() * Z[I]);
      }
      if ((Z[I].abs() + Z[I + 1].abs()) * LD[I].abs() < GAPTOL) {
        Z[I + 1] = Complex.zero;
        ISUPPZ[2] = I;
        break;
      }
      ZTZ.value += (Z[I + 1] * Z[I + 1]).toDouble();
    }
  }

  // Compute quantities for convergence test

  TMP = ONE / ZTZ.value;
  NRMINV.value = sqrt(TMP);
  RESID.value = (MINGMA.value).abs() * NRMINV.value;
  RQCORR.value = MINGMA.value * TMP;
}
