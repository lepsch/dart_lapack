// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zladiv.dart';

void ztgevc(
  final String SIDE,
  final String HOWMNY,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<Complex> S_,
  final int LDS,
  final Matrix<Complex> P_,
  final int LDP,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final int MM,
  final Box<int> M,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SELECT = SELECT_.having();
  final S = S_.having(ld: LDS);
  final P = P_.having(ld: LDP);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

  const ZERO = 0.0, ONE = 1.0;
  bool COMPL = false,
      COMPR = false,
      ILALL = false,
      ILBACK = false,
      ILBBAD,
      ILCOMP,
      LSA,
      LSB;
  int I, IBEG, IEIG, IEND, IHWMNY, IM, ISIDE, ISRC, J, JE, JR;
  double ACOEFA,
      ACOEFF,
      ANORM,
      ASCALE,
      BCOEFA,
      BIG,
      BIGNUM,
      BNORM,
      BSCALE,
      DMIN,
      SAFMIN,
      SBETA,
      SCALE,
      SMALL,
      TEMP,
      ULP,
      XMAX;
  Complex BCOEFF, CA, CB, D, SALPHA, SUM, SUMA, SUMB;

  // Decode and Test the input parameters

  if (lsame(HOWMNY, 'A')) {
    IHWMNY = 1;
    ILALL = true;
    ILBACK = false;
  } else if (lsame(HOWMNY, 'S')) {
    IHWMNY = 2;
    ILALL = false;
    ILBACK = false;
  } else if (lsame(HOWMNY, 'B')) {
    IHWMNY = 3;
    ILALL = true;
    ILBACK = true;
  } else {
    IHWMNY = -1;
  }

  if (lsame(SIDE, 'R')) {
    ISIDE = 1;
    COMPL = false;
    COMPR = true;
  } else if (lsame(SIDE, 'L')) {
    ISIDE = 2;
    COMPL = true;
    COMPR = false;
  } else if (lsame(SIDE, 'B')) {
    ISIDE = 3;
    COMPL = true;
    COMPR = true;
  } else {
    ISIDE = -1;
  }

  INFO.value = 0;
  if (ISIDE < 0) {
    INFO.value = -1;
  } else if (IHWMNY < 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDS < max(1, N)) {
    INFO.value = -6;
  } else if (LDP < max(1, N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZTGEVC', -INFO.value);
    return;
  }

  // Count the number of eigenvectors

  if (!ILALL) {
    IM = 0;
    for (J = 1; J <= N; J++) {
      if (SELECT[J]) IM++;
    }
  } else {
    IM = N;
  }

  // Check diagonal of B

  ILBBAD = false;
  for (J = 1; J <= N; J++) {
    if (P[J][J].imaginary != ZERO) ILBBAD = true;
  }

  if (ILBBAD) {
    INFO.value = -7;
  } else if (COMPL && LDVL < N || LDVL < 1) {
    INFO.value = -10;
  } else if (COMPR && LDVR < N || LDVR < 1) {
    INFO.value = -12;
  } else if (MM < IM) {
    INFO.value = -13;
  }
  if (INFO.value != 0) {
    xerbla('ZTGEVC', -INFO.value);
    return;
  }

  // Quick return if possible

  M.value = IM;
  if (N == 0) return;

  // Machine Constants

  SAFMIN = dlamch('Safe minimum');
  BIG = ONE / SAFMIN;
  ULP = dlamch('Epsilon') * dlamch('Base');
  SMALL = SAFMIN * N / ULP;
  BIG = ONE / SMALL;
  BIGNUM = ONE / (SAFMIN * N);

  // Compute the 1-norm of each column of the strictly upper triangular
  // part of A and B to check for possible overflow in the triangular
  // solver.

  ANORM = S[1][1].cabs1();
  BNORM = P[1][1].cabs1();
  RWORK[1] = ZERO;
  RWORK[N + 1] = ZERO;
  for (J = 2; J <= N; J++) {
    RWORK[J] = ZERO;
    RWORK[N + J] = ZERO;
    for (I = 1; I <= J - 1; I++) {
      RWORK[J] += S[I][J].cabs1();
      RWORK[N + J] += P[I][J].cabs1();
    }
    ANORM = max(ANORM, RWORK[J] + S[J][J].cabs1());
    BNORM = max(BNORM, RWORK[N + J] + P[J][J].cabs1());
  }

  ASCALE = ONE / max(ANORM, SAFMIN);
  BSCALE = ONE / max(BNORM, SAFMIN);

  // Left eigenvectors

  if (COMPL) {
    IEIG = 0;

    // Main loop over eigenvalues

    for (JE = 1; JE <= N; JE++) {
      if (ILALL) {
        ILCOMP = true;
      } else {
        ILCOMP = SELECT[JE];
      }
      if (ILCOMP) {
        IEIG++;

        if (S[JE][JE].cabs1() <= SAFMIN && P[JE][JE].real.abs() <= SAFMIN) {
          // Singular matrix pencil -- return unit eigenvector

          for (JR = 1; JR <= N; JR++) {
            VL[JR][IEIG] = Complex.zero;
          }
          VL[IEIG][IEIG] = Complex.one;
          continue;
        }

        // Non-singular eigenvalue:
        // Compute coefficients  a  and  b  in
        //      H
        //    y  ( a A - b B ) = 0

        TEMP = ONE /
            max(S[JE][JE].cabs1() * ASCALE,
                max(P[JE][JE].real.abs() * BSCALE, SAFMIN));
        SALPHA = TEMP.toComplex() * S[JE][JE] * ASCALE.toComplex();
        SBETA = (TEMP * P[JE][JE].real) * BSCALE;
        ACOEFF = SBETA * ASCALE;
        BCOEFF = SALPHA * BSCALE.toComplex();

        // Scale to avoid underflow

        LSA = SBETA.abs() >= SAFMIN && ACOEFF.abs() < SMALL;
        LSB = SALPHA.cabs1() >= SAFMIN && BCOEFF.cabs1() < SMALL;

        SCALE = ONE;
        if (LSA) SCALE = (SMALL / SBETA.abs()) * min(ANORM, BIG);
        if (LSB) SCALE = max(SCALE, (SMALL / SALPHA.cabs1()) * min(BNORM, BIG));
        if (LSA || LSB) {
          SCALE = min(SCALE,
              ONE / (SAFMIN * max(ONE, max(ACOEFF.abs(), BCOEFF.cabs1()))));
          if (LSA) {
            ACOEFF = ASCALE * (SCALE * SBETA);
          } else {
            ACOEFF = SCALE * ACOEFF;
          }
          if (LSB) {
            BCOEFF = BSCALE.toComplex() * (SCALE.toComplex() * SALPHA);
          } else {
            BCOEFF = SCALE.toComplex() * BCOEFF;
          }
        }

        ACOEFA = ACOEFF.abs();
        BCOEFA = BCOEFF.cabs1();
        XMAX = ONE;
        for (JR = 1; JR <= N; JR++) {
          WORK[JR] = Complex.zero;
        }
        WORK[JE] = Complex.one;
        DMIN = max(ULP * ACOEFA * ANORM, max(ULP * BCOEFA * BNORM, SAFMIN));

        // H
        // Triangular solve of  (a A - b B)  y = 0

        //                         H
        // (rowwise in  (a A - b B) , or columnwise in a A - b B)

        for (J = JE + 1; J <= N; J++) {
          // Compute
          //       j-1
          // SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
          //       k=je
          // (Scale if necessary)

          TEMP = ONE / XMAX;
          if (ACOEFA * RWORK[J] + BCOEFA * RWORK[N + J] > BIGNUM * TEMP) {
            for (JR = JE; JR <= J - 1; JR++) {
              WORK[JR] = TEMP.toComplex() * WORK[JR];
            }
            XMAX = ONE;
          }
          SUMA = Complex.zero;
          SUMB = Complex.zero;

          for (JR = JE; JR <= J - 1; JR++) {
            SUMA += S[JR][J].conjugate() * WORK[JR];
            SUMB += P[JR][J].conjugate() * WORK[JR];
          }
          SUM = ACOEFF.toComplex() * SUMA - BCOEFF.conjugate() * SUMB;

          // Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) )

          // with scaling and perturbation of the denominator

          D = (ACOEFF.toComplex() * S[J][J] - BCOEFF * P[J][J]).conjugate();
          if (D.cabs1() <= DMIN) D = DMIN.toComplex();

          if (D.cabs1() < ONE) {
            if (SUM.cabs1() >= BIGNUM * D.cabs1()) {
              TEMP = ONE / SUM.cabs1();
              for (JR = JE; JR <= J - 1; JR++) {
                WORK[JR] = TEMP.toComplex() * WORK[JR];
              }
              XMAX = TEMP * XMAX;
              SUM = TEMP.toComplex() * SUM;
            }
          }
          WORK[J] = zladiv(-SUM, D);
          XMAX = max(XMAX, WORK[J].cabs1());
        }

        // Back transform eigenvector if HOWMNY='B'.

        if (ILBACK) {
          zgemv('N', N, N + 1 - JE, Complex.one, VL(1, JE), LDVL, WORK(JE), 1,
              Complex.zero, WORK(N + 1), 1);
          ISRC = 2;
          IBEG = 1;
        } else {
          ISRC = 1;
          IBEG = JE;
        }

        // Copy and scale eigenvector into column of VL

        XMAX = ZERO;
        for (JR = IBEG; JR <= N; JR++) {
          XMAX = max(XMAX, WORK[(ISRC - 1) * N + JR].cabs1());
        }

        if (XMAX > SAFMIN) {
          TEMP = ONE / XMAX;
          for (JR = IBEG; JR <= N; JR++) {
            VL[JR][IEIG] = TEMP.toComplex() * WORK[(ISRC - 1) * N + JR];
          }
        } else {
          IBEG = N + 1;
        }

        for (JR = 1; JR <= IBEG - 1; JR++) {
          VL[JR][IEIG] = Complex.zero;
        }
      }
    }
  }

  // Right eigenvectors

  if (COMPR) {
    IEIG = IM + 1;

    // Main loop over eigenvalues

    for (JE = N; JE >= 1; JE--) {
      if (ILALL) {
        ILCOMP = true;
      } else {
        ILCOMP = SELECT[JE];
      }
      if (ILCOMP) {
        IEIG--;

        if (S[JE][JE].cabs1() <= SAFMIN && P[JE][JE].real.abs() <= SAFMIN) {
          // Singular matrix pencil -- return unit eigenvector

          for (JR = 1; JR <= N; JR++) {
            VR[JR][IEIG] = Complex.zero;
          }
          VR[IEIG][IEIG] = Complex.one;
          continue;
        }

        // Non-singular eigenvalue:
        // Compute coefficients  a  and  b  in

        // ( a A - b B ) x  = 0

        TEMP = ONE /
            max(S[JE][JE].cabs1() * ASCALE,
                max(P[JE][JE].real.abs() * BSCALE, SAFMIN));
        SALPHA = (TEMP.toComplex() * S[JE][JE]) * ASCALE.toComplex();
        SBETA = (TEMP * P[JE][JE].real) * BSCALE;
        ACOEFF = SBETA * ASCALE;
        BCOEFF = SALPHA * BSCALE.toComplex();

        // Scale to avoid underflow

        LSA = SBETA.abs() >= SAFMIN && ACOEFF.abs() < SMALL;
        LSB = SALPHA.cabs1() >= SAFMIN && BCOEFF.cabs1() < SMALL;

        SCALE = ONE;
        if (LSA) SCALE = (SMALL / SBETA.abs()) * min(ANORM, BIG);
        if (LSB) SCALE = max(SCALE, (SMALL / SALPHA.cabs1()) * min(BNORM, BIG));
        if (LSA || LSB) {
          SCALE = min(SCALE,
              ONE / (SAFMIN * max(ONE, max(ACOEFF.abs(), BCOEFF.cabs1()))));
          if (LSA) {
            ACOEFF = ASCALE * (SCALE * SBETA);
          } else {
            ACOEFF = SCALE * ACOEFF;
          }
          if (LSB) {
            BCOEFF = BSCALE.toComplex() * (SCALE.toComplex() * SALPHA);
          } else {
            BCOEFF = SCALE.toComplex() * BCOEFF;
          }
        }

        ACOEFA = ACOEFF.abs();
        BCOEFA = BCOEFF.cabs1();
        XMAX = ONE;
        for (JR = 1; JR <= N; JR++) {
          WORK[JR] = Complex.zero;
        }
        WORK[JE] = Complex.one;
        DMIN = max(ULP * ACOEFA * ANORM, max(ULP * BCOEFA * BNORM, SAFMIN));

        // Triangular solve of  (a A - b B) x = 0  (columnwise)

        // WORK(1:j-1) contains sums w,
        // WORK(j+1:JE) contains x

        for (JR = 1; JR <= JE - 1; JR++) {
          WORK[JR] = ACOEFF.toComplex() * S[JR][JE] - BCOEFF * P[JR][JE];
        }
        WORK[JE] = Complex.one;

        for (J = JE - 1; J >= 1; J--) {
          // Form x(j) := - w(j) / d
          // with scaling and perturbation of the denominator

          D = ACOEFF.toComplex() * S[J][J] - BCOEFF * P[J][J];
          if (D.cabs1() <= DMIN) D = DMIN.toComplex();

          if (D.cabs1() < ONE) {
            if (WORK[J].cabs1() >= BIGNUM * D.cabs1()) {
              TEMP = ONE / WORK[J].cabs1();
              for (JR = 1; JR <= JE; JR++) {
                WORK[JR] = TEMP.toComplex() * WORK[JR];
              }
            }
          }

          WORK[J] = zladiv(-WORK[J], D);

          if (J > 1) {
            // w += x(j)*(a S(*,j) - b P(*,j) ) with scaling

            if (WORK[J].cabs1() > ONE) {
              TEMP = ONE / WORK[J].cabs1();
              if (ACOEFA * RWORK[J] + BCOEFA * RWORK[N + J] >= BIGNUM * TEMP) {
                for (JR = 1; JR <= JE; JR++) {
                  WORK[JR] = TEMP.toComplex() * WORK[JR];
                }
              }
            }

            CA = ACOEFF.toComplex() * WORK[J];
            CB = BCOEFF * WORK[J];
            for (JR = 1; JR <= J - 1; JR++) {
              WORK[JR] += CA * S[JR][J] - CB * P[JR][J];
            }
          }
        }

        // Back transform eigenvector if HOWMNY='B'.

        if (ILBACK) {
          zgemv('N', N, JE, Complex.one, VR, LDVR, WORK, 1, Complex.zero,
              WORK(N + 1), 1);
          ISRC = 2;
          IEND = N;
        } else {
          ISRC = 1;
          IEND = JE;
        }

        // Copy and scale eigenvector into column of VR

        XMAX = ZERO;
        for (JR = 1; JR <= IEND; JR++) {
          XMAX = max(XMAX, WORK[(ISRC - 1) * N + JR].cabs1());
        }

        if (XMAX > SAFMIN) {
          TEMP = ONE / XMAX;
          for (JR = 1; JR <= IEND; JR++) {
            VR[JR][IEIG] = TEMP.toComplex() * WORK[(ISRC - 1) * N + JR];
          }
        } else {
          IEND = 0;
        }

        for (JR = IEND + 1; JR <= N; JR++) {
          VR[JR][IEIG] = Complex.zero;
        }
      }
    }
  }
}
