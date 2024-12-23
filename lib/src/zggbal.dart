// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/ddot.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/intrinsics/log10.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zggbal(
  final String JOB,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> ILO,
  final Box<int> IHI,
  final Array<double> LSCALE_,
  final Array<double> RSCALE_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final LSCALE = LSCALE_.having();
  final RSCALE = RSCALE_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  const THREE = 3.0, SCLFAC = 1.0e+1;
  int I = 0,
      ICAB,
      IFLOW = 0,
      IP1 = 0,
      IR,
      IRAB,
      IT,
      J = 0,
      JC,
      JP1 = 0,
      K,
      KOUNT,
      L,
      LCAB = 0,
      LM1 = 0,
      LRAB,
      LSFMAX = 0,
      LSFMIN,
      M = 0,
      NR,
      NRP2;
  double ALPHA,
      BASL,
      BETA,
      CAB,
      CMAX,
      COEF,
      COEF2,
      COEF5,
      COR,
      EW,
      EWC,
      GAMMA,
      PGAMMA = 0,
      RAB,
      SFMAX,
      SFMIN,
      SUM,
      T,
      TA,
      TB,
      TC;

  // Test the input parameters

  INFO.value = 0;
  if (!lsame(JOB, 'N') &&
      !lsame(JOB, 'P') &&
      !lsame(JOB, 'S') &&
      !lsame(JOB, 'B')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LDB < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZGGBAL', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) {
    ILO.value = 1;
    IHI.value = N;
    return;
  }

  if (N == 1) {
    ILO.value = 1;
    IHI.value = N;
    LSCALE[1] = ONE;
    RSCALE[1] = ONE;
    return;
  }

  if (lsame(JOB, 'N')) {
    ILO.value = 1;
    IHI.value = N;
    for (I = 1; I <= N; I++) {
      LSCALE[I] = ONE;
      RSCALE[I] = ONE;
    }
    return;
  }

  K = 1;
  L = N;
  if (!lsame(JOB, 'S')) {
    var firstRow = true;

    // Permute the matrices A and B to isolate the eigenvalues.

    // Find row with one nonzero in columns 1 through L
    permute:
    do {
      if (!firstRow) {
        L = LM1;
        if (L == 1) {
          RSCALE[1] = 1;
          LSCALE[1] = 1;
          break permute;
        }
        firstRow = false;
      }
      var nonZeroRowFound = true;
      LM1 = L - 1;
      findNonZero:
      for (I = L; I >= 1; I--) {
        var found = false;
        for (J = 1; J <= LM1; J++) {
          JP1 = J + 1;
          if (A[I][J] != Complex.zero || B[I][J] != Complex.zero) {
            found = true;
            break;
          }
        }
        if (!found) {
          J = L;
        } else {
          for (J = JP1; J <= L; J++) {
            if (A[I][J] != Complex.zero || B[I][J] != Complex.zero) {
              continue findNonZero;
            }
          }
          J = JP1 - 1;
        }
        M = L;
        IFLOW = 1;
        nonZeroRowFound = false;
        break;
      }
      var firstColumn = true;

      // Find column with one nonzero in rows K through N

      do {
        if (nonZeroRowFound) {
          if (!firstColumn) {
            K++;
          }
          firstColumn = false;

          var nonZeroColumnFound = true;
          findNonZero:
          for (J = K; J <= L; J++) {
            var found = false;
            for (I = K; I <= LM1; I++) {
              IP1 = I + 1;
              if (A[I][J] != Complex.zero || B[I][J] != Complex.zero) {
                found = true;
                break;
              }
            }
            if (!found) {
              I = L;
            } else {
              for (I = IP1; I <= L; I++) {
                if (A[I][J] != Complex.zero || B[I][J] != Complex.zero) {
                  continue findNonZero;
                }
              }
              I = IP1 - 1;
            }

            M = K;
            IFLOW = 2;
            nonZeroColumnFound = false;
            break;
          }
          if (nonZeroColumnFound) break permute;
        }
        nonZeroRowFound = true;

        // Permute rows M and I

        LSCALE[M] = I.toDouble();
        if (I != M) {
          zswap(N - K + 1, A(I, K).asArray(), LDA, A(M, K).asArray(), LDA);
          zswap(N - K + 1, B(I, K).asArray(), LDB, B(M, K).asArray(), LDB);
        }
        // Permute columns M and J

        RSCALE[M] = J.toDouble();
        if (J != M) {
          zswap(L, A(1, J).asArray(), 1, A(1, M).asArray(), 1);
          zswap(L, B(1, J).asArray(), 1, B(1, M).asArray(), 1);
        }
      } while (IFLOW == 2);
    } while (IFLOW == 1);
  }

  ILO.value = K;
  IHI.value = L;

  if (lsame(JOB, 'P')) {
    for (I = ILO.value; I <= IHI.value; I++) {
      LSCALE[I] = ONE;
      RSCALE[I] = ONE;
    }
    return;
  }

  if (ILO.value == IHI.value) return;

  // Balance the submatrix in rows ILO to IHI.

  NR = IHI.value - ILO.value + 1;
  for (I = ILO.value; I <= IHI.value; I++) {
    RSCALE[I] = ZERO;
    LSCALE[I] = ZERO;

    WORK[I] = ZERO;
    WORK[I + N] = ZERO;
    WORK[I + 2 * N] = ZERO;
    WORK[I + 3 * N] = ZERO;
    WORK[I + 4 * N] = ZERO;
    WORK[I + 5 * N] = ZERO;
  }

  // Compute right side vector in resulting linear equations

  BASL = log10(SCLFAC);
  for (I = ILO.value; I <= IHI.value; I++) {
    for (J = ILO.value; J <= IHI.value; J++) {
      if (A[I][J] == Complex.zero) {
        TA = ZERO;
      } else {
        TA = log10(A[I][J].cabs1()) / BASL;
      }

      if (B[I][J] == Complex.zero) {
        TB = ZERO;
      } else {
        TB = log10(B[I][J].cabs1()) / BASL;
      }

      WORK[I + 4 * N] -= TA + TB;
      WORK[J + 5 * N] -= TA + TB;
    }
  }

  COEF = ONE / (2 * NR);
  COEF2 = COEF * COEF;
  COEF5 = HALF * COEF2;
  NRP2 = NR + 2;
  BETA = ZERO;
  IT = 1;

  // Start generalized conjugate gradient iteration

  do {
    GAMMA = ddot(NR, WORK(ILO.value + 4 * N), 1, WORK(ILO.value + 4 * N), 1) +
        ddot(NR, WORK(ILO.value + 5 * N), 1, WORK(ILO.value + 5 * N), 1);

    EW = ZERO;
    EWC = ZERO;
    for (I = ILO.value; I <= IHI.value; I++) {
      EW += WORK[I + 4 * N];
      EWC += WORK[I + 5 * N];
    }

    GAMMA = COEF * GAMMA -
        COEF2 * (pow(EW, 2) + pow(EWC, 2)) -
        COEF5 * pow(EW - EWC, 2);
    if (GAMMA == ZERO) break;
    if (IT != 1) BETA = GAMMA / PGAMMA;
    T = COEF5 * (EWC - THREE * EW);
    TC = COEF5 * (EW - THREE * EWC);

    dscal(NR, BETA, WORK(ILO.value), 1);
    dscal(NR, BETA, WORK(ILO.value + N), 1);

    daxpy(NR, COEF, WORK(ILO.value + 4 * N), 1, WORK(ILO.value + N), 1);
    daxpy(NR, COEF, WORK(ILO.value + 5 * N), 1, WORK(ILO.value), 1);

    for (I = ILO.value; I <= IHI.value; I++) {
      WORK[I] += TC;
      WORK[I + N] += T;
    }

    // Apply matrix to vector

    for (I = ILO.value; I <= IHI.value; I++) {
      KOUNT = 0;
      SUM = ZERO;
      for (J = ILO.value; J <= IHI.value; J++) {
        if (A[I][J] != Complex.zero) {
          KOUNT++;
          SUM += WORK[J];
        }
        if (B[I][J] != Complex.zero) {
          KOUNT++;
          SUM += WORK[J];
        }
      }
      WORK[I + 2 * N] = KOUNT * WORK[I + N] + SUM;
    }

    for (J = ILO.value; J <= IHI.value; J++) {
      KOUNT = 0;
      SUM = ZERO;
      for (I = ILO.value; I <= IHI.value; I++) {
        if (A[I][J] != Complex.zero) {
          KOUNT++;
          SUM += WORK[I + N];
        }
        if (B[I][J] == Complex.zero) continue;
        KOUNT++;
        SUM += WORK[I + N];
      }
      WORK[J + 3 * N] = KOUNT * WORK[J] + SUM;
    }

    SUM = ddot(NR, WORK(ILO.value + N), 1, WORK(ILO.value + 2 * N), 1) +
        ddot(NR, WORK(ILO.value), 1, WORK(ILO.value + 3 * N), 1);
    ALPHA = GAMMA / SUM;

    // Determine correction to current iteration

    CMAX = ZERO;
    for (I = ILO.value; I <= IHI.value; I++) {
      COR = ALPHA * WORK[I + N];
      if (COR.abs() > CMAX) CMAX = COR.abs();
      LSCALE[I] += COR;
      COR = ALPHA * WORK[I];
      if (COR.abs() > CMAX) CMAX = COR.abs();
      RSCALE[I] += COR;
    }
    if (CMAX < HALF) break;

    daxpy(NR, -ALPHA, WORK(ILO.value + 2 * N), 1, WORK(ILO.value + 4 * N), 1);
    daxpy(NR, -ALPHA, WORK(ILO.value + 3 * N), 1, WORK(ILO.value + 5 * N), 1);

    PGAMMA = GAMMA;
    IT++;
  } while (IT <= NRP2);

  // End generalized conjugate gradient iteration

  SFMIN = dlamch('S');
  SFMAX = ONE / SFMIN;
  LSFMIN = (log10(SFMIN) / BASL + ONE).toInt();
  LSFMAX = log10(SFMAX) ~/ BASL;
  for (I = ILO.value; I <= IHI.value; I++) {
    IRAB = izamax(N - ILO.value + 1, A(I, ILO.value).asArray(), LDA);
    RAB = A[I][IRAB + ILO.value - 1].abs();
    IRAB = izamax(N - ILO.value + 1, B(I, ILO.value).asArray(), LDB);
    RAB = max(RAB, B[I][IRAB + ILO.value - 1].abs());
    LRAB = (log10(RAB + SFMIN) ~/ BASL + ONE).toInt();
    IR = (LSCALE[I] + sign(HALF, LSCALE[I])).toInt();
    IR = min(max(IR, LSFMIN), min(LSFMAX, LSFMAX - LRAB));
    LSCALE[I] = pow(SCLFAC, IR).toDouble();
    ICAB = izamax(IHI.value, A(1, I).asArray(), 1);
    CAB = A[ICAB][I].abs();
    ICAB = izamax(IHI.value, B(1, I).asArray(), 1);
    CAB = max(CAB, B[ICAB][I].abs());
    LCAB = (log10(CAB + SFMIN) / BASL + ONE).toInt();
    JC = (RSCALE[I] + sign(HALF, RSCALE[I])).toInt();
    JC = min(max(JC, LSFMIN), min(LSFMAX, LSFMAX - LCAB));
    RSCALE[I] = pow(SCLFAC, JC).toDouble();
  }

  // Row scaling of matrices A and B

  for (I = ILO.value; I <= IHI.value; I++) {
    zdscal(N - ILO.value + 1, LSCALE[I], A(I, ILO.value).asArray(), LDA);
    zdscal(N - ILO.value + 1, LSCALE[I], B(I, ILO.value).asArray(), LDB);
  }

  // Column scaling of matrices A and B

  for (J = ILO.value; J <= IHI.value; J++) {
    zdscal(IHI.value, RSCALE[J], A(1, J).asArray(), 1);
    zdscal(IHI.value, RSCALE[J], B(1, J).asArray(), 1);
  }
}
