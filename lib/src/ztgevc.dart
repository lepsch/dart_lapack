import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
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

  double ABS1(Complex X) => X.toDouble().abs() + X.imaginary.abs();

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
      // 10
      if (SELECT[J]) IM = IM + 1;
    } // 10
  } else {
    IM = N;
  }

  // Check diagonal of B

  ILBBAD = false;
  for (J = 1; J <= N; J++) {
    // 20
    if (P[J][J].imaginary != ZERO) ILBBAD = true;
  } // 20

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

  ANORM = ABS1(S[1][1]);
  BNORM = ABS1(P[1][1]);
  RWORK[1] = ZERO;
  RWORK[N + 1] = ZERO;
  for (J = 2; J <= N; J++) {
    // 40
    RWORK[J] = ZERO;
    RWORK[N + J] = ZERO;
    for (I = 1; I <= J - 1; I++) {
      // 30
      RWORK[J] = RWORK[J] + ABS1(S[I][J]);
      RWORK[N + J] = RWORK[N + J] + ABS1(P[I][J]);
    } // 30
    ANORM = max(ANORM, RWORK[J] + ABS1(S[J][J]));
    BNORM = max(BNORM, RWORK[N + J] + ABS1(P[J][J]));
  } // 40

  ASCALE = ONE / max(ANORM, SAFMIN);
  BSCALE = ONE / max(BNORM, SAFMIN);

  // Left eigenvectors

  if (COMPL) {
    IEIG = 0;

    // Main loop over eigenvalues

    for (JE = 1; JE <= N; JE++) {
      // 140
      if (ILALL) {
        ILCOMP = true;
      } else {
        ILCOMP = SELECT[JE];
      }
      if (ILCOMP) {
        IEIG++;

        if (ABS1(S[JE][JE]) <= SAFMIN && P[JE][JE].toDouble().abs() <= SAFMIN) {
          // Singular matrix pencil -- return unit eigenvector

          for (JR = 1; JR <= N; JR++) {
            // 50
            VL[JR][IEIG] = Complex.zero;
          } // 50
          VL[IEIG][IEIG] = Complex.one;
          continue;
        }

        // Non-singular eigenvalue:
        // Compute coefficients  a  and  b  in
        //      H
        //    y  ( a A - b B ) = 0

        TEMP = ONE /
            max(ABS1(S[JE][JE]) * ASCALE,
                max(P[JE][JE].toDouble().abs() * BSCALE, SAFMIN));
        SALPHA = TEMP.toComplex() * S[JE][JE] * ASCALE.toComplex();
        SBETA = (TEMP * (P[JE][JE]).toDouble()) * BSCALE;
        ACOEFF = SBETA * ASCALE;
        BCOEFF = SALPHA * BSCALE.toComplex();

        // Scale to avoid underflow

        LSA = (SBETA).abs() >= SAFMIN && (ACOEFF).abs() < SMALL;
        LSB = ABS1(SALPHA) >= SAFMIN && ABS1(BCOEFF) < SMALL;

        SCALE = ONE;
        if (LSA) SCALE = (SMALL / (SBETA).abs()) * min(ANORM, BIG);
        if (LSB) SCALE = max(SCALE, (SMALL / ABS1(SALPHA)) * min(BNORM, BIG));
        if (LSA || LSB) {
          SCALE = min(SCALE,
              ONE / (SAFMIN * max(ONE, max(ACOEFF.abs(), ABS1(BCOEFF)))));
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

        ACOEFA = (ACOEFF).abs();
        BCOEFA = ABS1(BCOEFF);
        XMAX = ONE;
        for (JR = 1; JR <= N; JR++) {
          // 60
          WORK[JR] = Complex.zero;
        } // 60
        WORK[JE] = Complex.one;
        DMIN = max(ULP * ACOEFA * ANORM, max(ULP * BCOEFA * BNORM, SAFMIN));

        // H
        // Triangular solve of  (a A - b B)  y = 0

        //                         H
        // (rowwise in  (a A - b B) , or columnwise in a A - b B)

        for (J = JE + 1; J <= N; J++) {
          // 100

          // Compute
          //       j-1
          // SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
          //       k=je
          // (Scale if necessary)

          TEMP = ONE / XMAX;
          if (ACOEFA * RWORK[J] + BCOEFA * RWORK[N + J] > BIGNUM * TEMP) {
            for (JR = JE; JR <= J - 1; JR++) {
              // 70
              WORK[JR] = TEMP.toComplex() * WORK[JR];
            } // 70
            XMAX = ONE;
          }
          SUMA = Complex.zero;
          SUMB = Complex.zero;

          for (JR = JE; JR <= J - 1; JR++) {
            // 80
            SUMA += S[JR][J].conjugate() * WORK[JR];
            SUMB += P[JR][J].conjugate() * WORK[JR];
          } // 80
          SUM = ACOEFF.toComplex() * SUMA - BCOEFF.conjugate() * SUMB;

          // Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) )

          // with scaling and perturbation of the denominator

          D = (ACOEFF.toComplex() * S[J][J] - BCOEFF * P[J][J]).conjugate();
          if (ABS1(D) <= DMIN) D = DMIN.toComplex();

          if (ABS1(D) < ONE) {
            if (ABS1(SUM) >= BIGNUM * ABS1(D)) {
              TEMP = ONE / ABS1(SUM);
              for (JR = JE; JR <= J - 1; JR++) {
                // 90
                WORK[JR] = TEMP.toComplex() * WORK[JR];
              } // 90
              XMAX = TEMP * XMAX;
              SUM = TEMP.toComplex() * SUM;
            }
          }
          WORK[J] = zladiv(-SUM, D);
          XMAX = max(XMAX, ABS1(WORK[J]));
        } // 100

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
          // 110
          XMAX = max(XMAX, ABS1(WORK[(ISRC - 1) * N + JR]));
        } // 110

        if (XMAX > SAFMIN) {
          TEMP = ONE / XMAX;
          for (JR = IBEG; JR <= N; JR++) {
            // 120
            VL[JR][IEIG] = TEMP.toComplex() * WORK[(ISRC - 1) * N + JR];
          } // 120
        } else {
          IBEG = N + 1;
        }

        for (JR = 1; JR <= IBEG - 1; JR++) {
          // 130
          VL[JR][IEIG] = Complex.zero;
        } // 130
      }
    } // 140
  }

  // Right eigenvectors

  if (COMPR) {
    IEIG = IM + 1;

    // Main loop over eigenvalues

    for (JE = N; JE >= 1; JE--) {
      // 250
      if (ILALL) {
        ILCOMP = true;
      } else {
        ILCOMP = SELECT[JE];
      }
      if (ILCOMP) {
        IEIG--;

        if (ABS1(S[JE][JE]) <= SAFMIN && P[JE][JE].toDouble().abs() <= SAFMIN) {
          // Singular matrix pencil -- return unit eigenvector

          for (JR = 1; JR <= N; JR++) {
            // 150
            VR[JR][IEIG] = Complex.zero;
          } // 150
          VR[IEIG][IEIG] = Complex.one;
          continue;
        }

        // Non-singular eigenvalue:
        // Compute coefficients  a  and  b  in

        // ( a A - b B ) x  = 0

        TEMP = ONE /
            max(ABS1(S[JE][JE]) * ASCALE,
                max(P[JE][JE].toDouble().abs() * BSCALE, SAFMIN));
        SALPHA = (TEMP.toComplex() * S[JE][JE]) * ASCALE.toComplex();
        SBETA = (TEMP * (P[JE][JE]).toDouble()) * BSCALE;
        ACOEFF = SBETA * ASCALE;
        BCOEFF = SALPHA * BSCALE.toComplex();

        // Scale to avoid underflow

        LSA = (SBETA).abs() >= SAFMIN && (ACOEFF).abs() < SMALL;
        LSB = ABS1(SALPHA) >= SAFMIN && ABS1(BCOEFF) < SMALL;

        SCALE = ONE;
        if (LSA) SCALE = (SMALL / (SBETA).abs()) * min(ANORM, BIG);
        if (LSB) SCALE = max(SCALE, (SMALL / ABS1(SALPHA)) * min(BNORM, BIG));
        if (LSA || LSB) {
          SCALE = min(SCALE,
              ONE / (SAFMIN * max(ONE, max(ACOEFF.abs(), ABS1(BCOEFF)))));
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

        ACOEFA = (ACOEFF).abs();
        BCOEFA = ABS1(BCOEFF);
        XMAX = ONE;
        for (JR = 1; JR <= N; JR++) {
          // 160
          WORK[JR] = Complex.zero;
        } // 160
        WORK[JE] = Complex.one;
        DMIN = max(ULP * ACOEFA * ANORM, max(ULP * BCOEFA * BNORM, SAFMIN));

        // Triangular solve of  (a A - b B) x = 0  (columnwise)

        // WORK(1:j-1) contains sums w,
        // WORK(j+1:JE) contains x

        for (JR = 1; JR <= JE - 1; JR++) {
          // 170
          WORK[JR] = ACOEFF.toComplex() * S[JR][JE] - BCOEFF * P[JR][JE];
        } // 170
        WORK[JE] = Complex.one;

        for (J = JE - 1; J >= 1; J--) {
          // 210

          // Form x(j) := - w(j) / d
          // with scaling and perturbation of the denominator

          D = ACOEFF.toComplex() * S[J][J] - BCOEFF * P[J][J];
          if (ABS1(D) <= DMIN) D = DMIN.toComplex();

          if (ABS1(D) < ONE) {
            if (ABS1(WORK[J]) >= BIGNUM * ABS1(D)) {
              TEMP = ONE / ABS1(WORK[J]);
              for (JR = 1; JR <= JE; JR++) {
                // 180
                WORK[JR] = TEMP.toComplex() * WORK[JR];
              } // 180
            }
          }

          WORK[J] = zladiv(-WORK[J], D);

          if (J > 1) {
            // w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling

            if (ABS1(WORK[J]) > ONE) {
              TEMP = ONE / ABS1(WORK[J]);
              if (ACOEFA * RWORK[J] + BCOEFA * RWORK[N + J] >= BIGNUM * TEMP) {
                for (JR = 1; JR <= JE; JR++) {
                  // 190
                  WORK[JR] = TEMP.toComplex() * WORK[JR];
                } // 190
              }
            }

            CA = ACOEFF.toComplex() * WORK[J];
            CB = BCOEFF * WORK[J];
            for (JR = 1; JR <= J - 1; JR++) {
              // 200
              WORK[JR] = WORK[JR] + CA * S[JR][J] - CB * P[JR][J];
            } // 200
          }
        } // 210

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
          // 220
          XMAX = max(XMAX, ABS1(WORK[(ISRC - 1) * N + JR]));
        } // 220

        if (XMAX > SAFMIN) {
          TEMP = ONE / XMAX;
          for (JR = 1; JR <= IEND; JR++) {
            // 230
            VR[JR][IEIG] = TEMP.toComplex() * WORK[(ISRC - 1) * N + JR];
          } // 230
        } else {
          IEND = 0;
        }

        for (JR = IEND + 1; JR <= N; JR++) {
          // 240
          VR[JR][IEIG] = Complex.zero;
        } // 240
      }
    } // 250
  }
}
