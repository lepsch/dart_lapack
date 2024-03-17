import 'dart:math';

import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlag2.dart';
import 'package:lapack/src/dlaln2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dtgevc(
  final String SIDE,
  final String HOWMNY,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<double> S_,
  final int LDS,
  final Matrix<double> P_,
  final int LDP,
  final Matrix<double> VL_,
  final int LDVL,
  final Matrix<double> VR_,
  final int LDVR,
  final int MM,
  final Box<int> M,
  final Array<double> WORK_,
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
  const ZERO = 0.0, ONE = 1.0, SAFETY = 1.0e+2;
  bool COMPL = false,
      COMPR = false,
      IL2BY2,
      ILABAD,
      ILALL,
      ILBACK = false,
      ILBBAD,
      ILCOMP,
      ILCPLX,
      LSA,
      LSB;
  int I, IBEG, IEIG, IEND, IHWMNY, IM, ISIDE, J, JA, JC, JE, JR, JW, NA, NW;
  double ACOEFA,
      ANORM,
      ASCALE,
      BCOEFA,
      BIG,
      BIGNUM,
      BNORM,
      BSCALE,
      CIM2A,
      CIM2B,
      CIMAGA,
      CIMAGB,
      CRE2A,
      CRE2B,
      CREALA,
      CREALB,
      DMIN,
      SAFMIN = 0,
      SALFAR,
      SBETA,
      SMALL,
      TEMP2I,
      TEMP2R,
      ULP,
      XMAX,
      XSCALE;
  final IINFO = Box(0);
  final ACOEF = Box(0.0),
      TEMP = Box(0.0),
      TEMP2 = Box(0.0),
      BCOEFR = Box(0.0),
      BCOEFI = Box(0.0),
      SCALE = Box(0.0);
  final BDIAG = Array<double>(2);
  final SUM = Matrix<double>(2, 2),
      SUMS = Matrix<double>(2, 2),
      SUMP = Matrix<double>(2, 2);

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
    ILALL = true;
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
    xerbla('DTGEVC', -INFO.value);
    return;
  }

  // Count the number of eigenvectors to be computed

  if (!ILALL) {
    IM = 0;
    ILCPLX = false;
    for (J = 1; J <= N; J++) {
      if (ILCPLX) {
        ILCPLX = false;
        continue;
      }
      if (J < N) {
        if (S[J + 1][J] != ZERO) ILCPLX = true;
      }
      if (ILCPLX) {
        if (SELECT[J] || SELECT[J + 1]) IM += 2;
      } else {
        if (SELECT[J]) IM++;
      }
    }
  } else {
    IM = N;
  }

  // Check 2-by-2 diagonal blocks of A, B

  ILABAD = false;
  ILBBAD = false;
  for (J = 1; J <= N - 1; J++) {
    if (S[J + 1][J] != ZERO) {
      if (P[J][J] == ZERO || P[J + 1][J + 1] == ZERO || P[J][J + 1] != ZERO) {
        ILBBAD = true;
      }
      if (J < N - 1) {
        if (S[J + 2][J + 1] != ZERO) ILABAD = true;
      }
    }
  }

  if (ILABAD) {
    INFO.value = -5;
  } else if (ILBBAD) {
    INFO.value = -7;
  } else if (COMPL && LDVL < N || LDVL < 1) {
    INFO.value = -10;
  } else if (COMPR && LDVR < N || LDVR < 1) {
    INFO.value = -12;
  } else if (MM < IM) {
    INFO.value = -13;
  }
  if (INFO.value != 0) {
    xerbla('DTGEVC', -INFO.value);
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
  // part (i.e., excluding all elements belonging to the diagonal
  // blocks) of A and B to check for possible overflow in the
  // triangular solver.

  ANORM = S[1][1].abs();
  if (N > 1) ANORM += S[2][1].abs();
  BNORM = P[1][1].abs();
  WORK[1] = ZERO;
  WORK[N + 1] = ZERO;

  for (J = 2; J <= N; J++) {
    TEMP.value = ZERO;
    TEMP2.value = ZERO;
    if (S[J][J - 1] == ZERO) {
      IEND = J - 1;
    } else {
      IEND = J - 2;
    }
    for (I = 1; I <= IEND; I++) {
      TEMP.value += S[I][J].abs();
      TEMP2.value += P[I][J].abs();
    }
    WORK[J] = TEMP.value;
    WORK[N + J] = TEMP2.value;
    for (I = IEND + 1; I <= min(J + 1, N); I++) {
      TEMP.value += S[I][J].abs();
      TEMP2.value += P[I][J].abs();
    }
    ANORM = max(ANORM, TEMP.value);
    BNORM = max(BNORM, TEMP2.value);
  }

  ASCALE = ONE / max(ANORM, SAFMIN);
  BSCALE = ONE / max(BNORM, SAFMIN);

  // Left eigenvectors

  if (COMPL) {
    IEIG = 0;

    // Main loop over eigenvalues

    ILCPLX = false;
    for (JE = 1; JE <= N; JE++) {
      // Skip this iteration if (a) HOWMNY='S' and SELECT= false , or
      // (b) this would be the second of a complex pair.
      // Check for complex eigenvalue, so as to be sure of which
      // entry(-ies) of SELECT to look at.

      if (ILCPLX) {
        ILCPLX = false;
        continue;
      }
      NW = 1;
      if (JE < N) {
        if (S[JE + 1][JE] != ZERO) {
          ILCPLX = true;
          NW = 2;
        }
      }
      if (ILALL) {
        ILCOMP = true;
      } else if (ILCPLX) {
        ILCOMP = SELECT[JE] || SELECT[JE + 1];
      } else {
        ILCOMP = SELECT[JE];
      }
      if (!ILCOMP) continue;

      // Decide if (a) singular pencil, (b) real eigenvalue, or
      // (c) complex eigenvalue.

      if (!ILCPLX) {
        if (S[JE][JE].abs() <= SAFMIN && P[JE][JE].abs() <= SAFMIN) {
          // Singular matrix pencil -- return unit eigenvector

          IEIG++;
          for (JR = 1; JR <= N; JR++) {
            VL[JR][IEIG] = ZERO;
          }
          VL[IEIG][IEIG] = ONE;
          continue;
        }
      }

      // Clear vector

      for (JR = 1; JR <= NW * N; JR++) {
        WORK[2 * N + JR] = ZERO;
      }
      // T
      // Compute coefficients in  ( a A - b B )  y = 0
      // a  is  ACOEF.value
      // b  is  BCOEFR.value + i*BCOEFI.value

      if (!ILCPLX) {
        // Real eigenvalue

        TEMP.value = ONE /
            max(S[JE][JE].abs() * ASCALE,
                max(P[JE][JE].abs() * BSCALE, SAFMIN));
        SALFAR = (TEMP.value * S[JE][JE]) * ASCALE;
        SBETA = (TEMP.value * P[JE][JE]) * BSCALE;
        ACOEF.value = SBETA * ASCALE;
        BCOEFR.value = SALFAR * BSCALE;
        BCOEFI.value = ZERO;

        // Scale to avoid underflow

        SCALE.value = ONE;
        LSA = SBETA.abs() >= SAFMIN && (ACOEF.value).abs() < SMALL;
        LSB = SALFAR.abs() >= SAFMIN && (BCOEFR.value).abs() < SMALL;
        if (LSA) SCALE.value = (SMALL / SBETA.abs()) * min(ANORM, BIG);
        if (LSB) {
          SCALE.value =
              max(SCALE.value, (SMALL / SALFAR.abs()) * min(BNORM, BIG));
        }
        if (LSA || LSB) {
          SCALE.value = min(
              SCALE.value,
              ONE /
                  (SAFMIN *
                      max(ONE, max(ACOEF.value.abs(), BCOEFR.value.abs()))));
          if (LSA) {
            ACOEF.value = ASCALE * (SCALE.value * SBETA);
          } else {
            ACOEF.value = SCALE.value * ACOEF.value;
          }
          if (LSB) {
            BCOEFR.value = BSCALE * (SCALE.value * SALFAR);
          } else {
            BCOEFR.value = SCALE.value * BCOEFR.value;
          }
        }
        ACOEFA = (ACOEF.value).abs();
        BCOEFA = (BCOEFR.value).abs();

        // First component is 1

        WORK[2 * N + JE] = ONE;
        XMAX = ONE;
      } else {
        // Complex eigenvalue

        dlag2(S(JE, JE), LDS, P(JE, JE), LDP, SAFMIN * SAFETY, ACOEF, TEMP,
            BCOEFR, TEMP2, BCOEFI);
        BCOEFI.value = -BCOEFI.value;
        if (BCOEFI.value == ZERO) {
          INFO.value = JE;
          return;
        }

        // Scale to avoid over/underflow

        ACOEFA = (ACOEF.value).abs();
        BCOEFA = (BCOEFR.value).abs() + (BCOEFI.value).abs();
        SCALE.value = ONE;
        if (ACOEFA * ULP < SAFMIN && ACOEFA >= SAFMIN) {
          SCALE.value = (SAFMIN / ULP) / ACOEFA;
        }
        if (BCOEFA * ULP < SAFMIN && BCOEFA >= SAFMIN) {
          SCALE.value = max(SCALE.value, (SAFMIN / ULP) / BCOEFA);
        }
        if (SAFMIN * ACOEFA > ASCALE) SCALE.value = ASCALE / (SAFMIN * ACOEFA);
        if (SAFMIN * BCOEFA > BSCALE) {
          SCALE.value = min(SCALE.value, BSCALE / (SAFMIN * BCOEFA));
        }
        if (SCALE.value != ONE) {
          ACOEF.value = SCALE.value * ACOEF.value;
          ACOEFA = (ACOEF.value).abs();
          BCOEFR.value = SCALE.value * BCOEFR.value;
          BCOEFI.value = SCALE.value * BCOEFI.value;
          BCOEFA = (BCOEFR.value).abs() + (BCOEFI.value).abs();
        }

        // Compute first two components of eigenvector

        TEMP.value = ACOEF.value * S[JE + 1][JE];
        TEMP2R = ACOEF.value * S[JE][JE] - BCOEFR.value * P[JE][JE];
        TEMP2I = -BCOEFI.value * P[JE][JE];
        if ((TEMP.value).abs() > TEMP2R.abs() + TEMP2I.abs()) {
          WORK[2 * N + JE] = ONE;
          WORK[3 * N + JE] = ZERO;
          WORK[2 * N + JE + 1] = -TEMP2R / TEMP.value;
          WORK[3 * N + JE + 1] = -TEMP2I / TEMP.value;
        } else {
          WORK[2 * N + JE + 1] = ONE;
          WORK[3 * N + JE + 1] = ZERO;
          TEMP.value = ACOEF.value * S[JE][JE + 1];
          WORK[2 * N + JE] = (BCOEFR.value * P[JE + 1][JE + 1] -
                  ACOEF.value * S[JE + 1][JE + 1]) /
              TEMP.value;
          WORK[3 * N + JE] = BCOEFI.value * P[JE + 1][JE + 1] / TEMP.value;
        }
        XMAX = max(
          WORK[2 * N + JE].abs() + WORK[3 * N + JE].abs(),
          WORK[2 * N + JE + 1].abs() + (WORK[3 * N + JE + 1]),
        ).abs();
      }

      DMIN = max(ULP * ACOEFA * ANORM, max(ULP * BCOEFA * BNORM, SAFMIN));

      // T
      // Triangular solve of  (a A - b B)  y = 0

      // T
      // (rowwise in  (a A - b B) , or columnwise in (a A - b B) )

      IL2BY2 = false;

      for (J = JE + NW; J <= N; J++) {
        if (IL2BY2) {
          IL2BY2 = false;
          continue;
        }

        NA = 1;
        BDIAG[1] = P[J][J];
        if (J < N) {
          if (S[J + 1][J] != ZERO) {
            IL2BY2 = true;
            BDIAG[2] = P[J + 1][J + 1];
            NA = 2;
          }
        }

        // Check whether scaling is necessary for dot products

        XSCALE = ONE / max(ONE, XMAX);
        TEMP.value = max(
            WORK[J], max(WORK[N + J], ACOEFA * WORK[J] + BCOEFA * WORK[N + J]));
        if (IL2BY2) {
          TEMP.value = max(
              TEMP.value,
              max(
                WORK[J + 1],
                max(
                  WORK[N + J + 1],
                  ACOEFA * WORK[J + 1] + BCOEFA * WORK[N + J + 1],
                ),
              ));
        }
        if (TEMP.value > BIGNUM * XSCALE) {
          for (JW = 0; JW <= NW - 1; JW++) {
            for (JR = JE; JR <= J - 1; JR++) {
              WORK[(JW + 2) * N + JR] = XSCALE * WORK[(JW + 2) * N + JR];
            }
          }
          XMAX *= XSCALE;
        }

        // Compute dot products
        //
        //       j-1
        // SUM = sum  conjg( a*S[k][j] - b*P[k][j] )*x(k)
        //       k=je
        //
        // To reduce the op count, this is done as
        //
        // _        j-1                  _        j-1
        // a*conjg( sum  S[k][j]*x(k) ) - b*conjg( sum  P[k][j]*x(k) )
        //          k=je                          k=je
        //
        // which may cause underflow problems if A or B are close
        // to underflow.  (E.g., less than SMALL.)

        for (JW = 1; JW <= NW; JW++) {
          for (JA = 1; JA <= NA; JA++) {
            SUMS[JA][JW] = ZERO;
            SUMP[JA][JW] = ZERO;

            for (JR = JE; JR <= J - 1; JR++) {
              SUMS[JA][JW] =
                  SUMS[JA][JW] + S[JR][J + JA - 1] * WORK[(JW + 1) * N + JR];
              SUMP[JA][JW] =
                  SUMP[JA][JW] + P[JR][J + JA - 1] * WORK[(JW + 1) * N + JR];
            }
          }
        }

        for (JA = 1; JA <= NA; JA++) {
          if (ILCPLX) {
            SUM[JA][1] = -ACOEF.value * SUMS[JA][1] +
                BCOEFR.value * SUMP[JA][1] -
                BCOEFI.value * SUMP[JA][2];
            SUM[JA][2] = -ACOEF.value * SUMS[JA][2] +
                BCOEFR.value * SUMP[JA][2] +
                BCOEFI.value * SUMP[JA][1];
          } else {
            SUM[JA][1] =
                -ACOEF.value * SUMS[JA][1] + BCOEFR.value * SUMP[JA][1];
          }
        }

        // T
        // Solve  ( a A - b B )  y = SUM(,)
        // with scaling and perturbation of the denominator

        dlaln2(
            true,
            NA,
            NW,
            DMIN,
            ACOEF.value,
            S(J, J),
            LDS,
            BDIAG[1],
            BDIAG[2],
            SUM,
            2,
            BCOEFR.value,
            BCOEFI.value,
            WORK(2 * N + J).asMatrix(N),
            N,
            SCALE,
            TEMP,
            IINFO);
        if (SCALE.value < ONE) {
          for (JW = 0; JW <= NW - 1; JW++) {
            for (JR = JE; JR <= J - 1; JR++) {
              WORK[(JW + 2) * N + JR] = SCALE.value * WORK[(JW + 2) * N + JR];
            }
          }
          XMAX = SCALE.value * XMAX;
        }
        XMAX = max(XMAX, TEMP.value);
      }

      // Copy eigenvector to VL, back transforming if
      // HOWMNY='B'.

      IEIG++;
      if (ILBACK) {
        for (JW = 0; JW <= NW - 1; JW++) {
          dgemv('N', N, N + 1 - JE, ONE, VL(1, JE), LDVL,
              WORK((JW + 2) * N + JE), 1, ZERO, WORK((JW + 4) * N + 1), 1);
        }
        dlacpy(' ', N, NW, WORK(4 * N + 1).asMatrix(N), N, VL(1, JE), LDVL);
        IBEG = 1;
      } else {
        dlacpy(' ', N, NW, WORK(2 * N + 1).asMatrix(N), N, VL(1, IEIG), LDVL);
        IBEG = JE;
      }

      // Scale eigenvector

      XMAX = ZERO;
      if (ILCPLX) {
        for (J = IBEG; J <= N; J++) {
          XMAX = max(XMAX, VL[J][IEIG].abs() + VL[J][IEIG + 1].abs());
        }
      } else {
        for (J = IBEG; J <= N; J++) {
          XMAX = max(XMAX, VL[J][IEIG].abs());
        }
      }

      if (XMAX > SAFMIN) {
        XSCALE = ONE / XMAX;

        for (JW = 0; JW <= NW - 1; JW++) {
          for (JR = IBEG; JR <= N; JR++) {
            VL[JR][IEIG + JW] = XSCALE * VL[JR][IEIG + JW];
          }
        }
      }
      IEIG += NW - 1;
    }
  }

  // Right eigenvectors

  if (COMPR) {
    IEIG = IM + 1;

    // Main loop over eigenvalues

    ILCPLX = false;
    for (JE = N; JE >= 1; JE--) {
      // Skip this iteration if (a) HOWMNY='S' and SELECT= false , or
      // (b) this would be the second of a complex pair.
      // Check for complex eigenvalue, so as to be sure of which
      // entry(-ies) of SELECT to look at -- if complex, SELECT[JE]
      // or SELECT[JE-1].
      // If this is a complex pair, the 2-by-2 diagonal block
      // corresponding to the eigenvalue is in rows/columns JE-1:JE

      if (ILCPLX) {
        ILCPLX = false;
        continue;
      }
      NW = 1;
      if (JE > 1) {
        if (S[JE][JE - 1] != ZERO) {
          ILCPLX = true;
          NW = 2;
        }
      }
      if (ILALL) {
        ILCOMP = true;
      } else if (ILCPLX) {
        ILCOMP = SELECT[JE] || SELECT[JE - 1];
      } else {
        ILCOMP = SELECT[JE];
      }
      if (!ILCOMP) continue;

      // Decide if (a) singular pencil, (b) real eigenvalue, or
      // (c) complex eigenvalue.

      if (!ILCPLX) {
        if (S[JE][JE].abs() <= SAFMIN && P[JE][JE].abs() <= SAFMIN) {
          // Singular matrix pencil -- unit eigenvector

          IEIG--;
          for (JR = 1; JR <= N; JR++) {
            VR[JR][IEIG] = ZERO;
          }
          VR[IEIG][IEIG] = ONE;
          continue;
        }
      }

      // Clear vector

      for (JW = 0; JW <= NW - 1; JW++) {
        for (JR = 1; JR <= N; JR++) {
          WORK[(JW + 2) * N + JR] = ZERO;
        }
      }

      // Compute coefficients in  ( a A - b B ) x = 0
      // a  is  ACOEF.value
      // b  is  BCOEFR.value + i*BCOEFI.value

      if (!ILCPLX) {
        // Real eigenvalue

        TEMP.value = ONE /
            max(S[JE][JE].abs() * ASCALE,
                max(P[JE][JE].abs() * BSCALE, SAFMIN));
        SALFAR = (TEMP.value * S[JE][JE]) * ASCALE;
        SBETA = (TEMP.value * P[JE][JE]) * BSCALE;
        ACOEF.value = SBETA * ASCALE;
        BCOEFR.value = SALFAR * BSCALE;
        BCOEFI.value = ZERO;

        // Scale to avoid underflow

        SCALE.value = ONE;
        LSA = SBETA.abs() >= SAFMIN && (ACOEF.value).abs() < SMALL;
        LSB = SALFAR.abs() >= SAFMIN && (BCOEFR.value).abs() < SMALL;
        if (LSA) SCALE.value = (SMALL / SBETA.abs()) * min(ANORM, BIG);
        if (LSB) {
          SCALE.value =
              max(SCALE.value, (SMALL / SALFAR.abs()) * min(BNORM, BIG));
        }
        if (LSA || LSB) {
          SCALE.value = min(
              SCALE.value,
              ONE /
                  SAFMIN *
                  max(ONE, max(ACOEF.value.abs(), BCOEFR.value.abs())));
          if (LSA) {
            ACOEF.value = ASCALE * (SCALE.value * SBETA);
          } else {
            ACOEF.value = SCALE.value * ACOEF.value;
          }
          if (LSB) {
            BCOEFR.value = BSCALE * (SCALE.value * SALFAR);
          } else {
            BCOEFR.value = SCALE.value * BCOEFR.value;
          }
        }
        ACOEFA = (ACOEF.value).abs();
        BCOEFA = (BCOEFR.value).abs();

        // First component is 1

        WORK[2 * N + JE] = ONE;
        XMAX = ONE;

        // Compute contribution from column JE of A and B to sum
        // (See "Further Details", above.)

        for (JR = 1; JR <= JE - 1; JR++) {
          WORK[2 * N + JR] = BCOEFR.value * P[JR][JE] - ACOEF.value * S[JR][JE];
        }
      } else {
        // Complex eigenvalue

        dlag2(S(JE - 1, JE - 1), LDS, P(JE - 1, JE - 1), LDP, SAFMIN * SAFETY,
            ACOEF, TEMP, BCOEFR, TEMP2, BCOEFI);
        if (BCOEFI.value == ZERO) {
          INFO.value = JE - 1;
          return;
        }

        // Scale to avoid over/underflow

        ACOEFA = (ACOEF.value).abs();
        BCOEFA = (BCOEFR.value).abs() + (BCOEFI.value).abs();
        SCALE.value = ONE;
        if (ACOEFA * ULP < SAFMIN && ACOEFA >= SAFMIN) {
          SCALE.value = (SAFMIN / ULP) / ACOEFA;
        }
        if (BCOEFA * ULP < SAFMIN && BCOEFA >= SAFMIN) {
          SCALE.value = max(SCALE.value, (SAFMIN / ULP) / BCOEFA);
        }
        if (SAFMIN * ACOEFA > ASCALE) SCALE.value = ASCALE / (SAFMIN * ACOEFA);
        if (SAFMIN * BCOEFA > BSCALE) {
          SCALE.value = min(SCALE.value, BSCALE / (SAFMIN * BCOEFA));
        }
        if (SCALE.value != ONE) {
          ACOEF.value = SCALE.value * ACOEF.value;
          ACOEFA = (ACOEF.value).abs();
          BCOEFR.value = SCALE.value * BCOEFR.value;
          BCOEFI.value = SCALE.value * BCOEFI.value;
          BCOEFA = (BCOEFR.value).abs() + (BCOEFI.value).abs();
        }

        // Compute first two components of eigenvector
        // and contribution to sums

        TEMP.value = ACOEF.value * S[JE][JE - 1];
        TEMP2R = ACOEF.value * S[JE][JE] - BCOEFR.value * P[JE][JE];
        TEMP2I = -BCOEFI.value * P[JE][JE];
        if ((TEMP.value).abs() >= TEMP2R.abs() + TEMP2I.abs()) {
          WORK[2 * N + JE] = ONE;
          WORK[3 * N + JE] = ZERO;
          WORK[2 * N + JE - 1] = -TEMP2R / TEMP.value;
          WORK[3 * N + JE - 1] = -TEMP2I / TEMP.value;
        } else {
          WORK[2 * N + JE - 1] = ONE;
          WORK[3 * N + JE - 1] = ZERO;
          TEMP.value = ACOEF.value * S[JE - 1][JE];
          WORK[2 * N + JE] = (BCOEFR.value * P[JE - 1][JE - 1] -
                  ACOEF.value * S[JE - 1][JE - 1]) /
              TEMP.value;
          WORK[3 * N + JE] = BCOEFI.value * P[JE - 1][JE - 1] / TEMP.value;
        }

        XMAX = max(
          WORK[2 * N + JE].abs() + WORK[3 * N + JE].abs(),
          WORK[2 * N + JE - 1].abs() + (WORK[3 * N + JE - 1]),
        ).abs();

        // Compute contribution from columns JE and JE-1
        // of A and B to the sums.

        CREALA = ACOEF.value * WORK[2 * N + JE - 1];
        CIMAGA = ACOEF.value * WORK[3 * N + JE - 1];
        CREALB = BCOEFR.value * WORK[2 * N + JE - 1] -
            BCOEFI.value * WORK[3 * N + JE - 1];
        CIMAGB = BCOEFI.value * WORK[2 * N + JE - 1] +
            BCOEFR.value * WORK[3 * N + JE - 1];
        CRE2A = ACOEF.value * WORK[2 * N + JE];
        CIM2A = ACOEF.value * WORK[3 * N + JE];
        CRE2B =
            BCOEFR.value * WORK[2 * N + JE] - BCOEFI.value * WORK[3 * N + JE];
        CIM2B =
            BCOEFI.value * WORK[2 * N + JE] + BCOEFR.value * WORK[3 * N + JE];
        for (JR = 1; JR <= JE - 2; JR++) {
          WORK[2 * N + JR] = -CREALA * S[JR][JE - 1] +
              CREALB * P[JR][JE - 1] -
              CRE2A * S[JR][JE] +
              CRE2B * P[JR][JE];
          WORK[3 * N + JR] = -CIMAGA * S[JR][JE - 1] +
              CIMAGB * P[JR][JE - 1] -
              CIM2A * S[JR][JE] +
              CIM2B * P[JR][JE];
        }
      }

      DMIN = max(ULP * ACOEFA * ANORM, max(ULP * BCOEFA * BNORM, SAFMIN));

      // Columnwise triangular solve of  (a A - b B)  x = 0

      IL2BY2 = false;
      for (J = JE - NW; J >= 1; J--) {
        // If a 2-by-2 block, is in position j-1:j, wait until
        // next iteration to process it (when it will be j:j+1)

        if (!IL2BY2 && J > 1) {
          if (S[J][J - 1] != ZERO) {
            IL2BY2 = true;
            continue;
          }
        }
        BDIAG[1] = P[J][J];
        if (IL2BY2) {
          NA = 2;
          BDIAG[2] = P[J + 1][J + 1];
        } else {
          NA = 1;
        }

        // Compute x(j) (and x(j+1), if 2-by-2 block)

        dlaln2(
            false,
            NA,
            NW,
            DMIN,
            ACOEF.value,
            S(J, J),
            LDS,
            BDIAG[1],
            BDIAG[2],
            WORK(2 * N + J).asMatrix(N),
            N,
            BCOEFR.value,
            BCOEFI.value,
            SUM,
            2,
            SCALE,
            TEMP,
            IINFO);
        if (SCALE.value < ONE) {
          for (JW = 0; JW <= NW - 1; JW++) {
            for (JR = 1; JR <= JE; JR++) {
              WORK[(JW + 2) * N + JR] = SCALE.value * WORK[(JW + 2) * N + JR];
            }
          }
        }
        XMAX = max(SCALE.value * XMAX, TEMP.value);

        for (JW = 1; JW <= NW; JW++) {
          for (JA = 1; JA <= NA; JA++) {
            WORK[(JW + 1) * N + J + JA - 1] = SUM[JA][JW];
          }
        }

        // w += x(j)*(a S[*][j] - b P[*][j] ) with scaling

        if (J > 1) {
          // Check whether scaling is necessary for sum.

          XSCALE = ONE / max(ONE, XMAX);
          TEMP.value = ACOEFA * WORK[J] + BCOEFA * WORK[N + J];
          if (IL2BY2) {
            TEMP.value = max(
                TEMP.value, ACOEFA * WORK[J + 1] + BCOEFA * WORK[N + J + 1]);
          }
          TEMP.value = max(TEMP.value, max(ACOEFA, BCOEFA));
          if (TEMP.value > BIGNUM * XSCALE) {
            for (JW = 0; JW <= NW - 1; JW++) {
              for (JR = 1; JR <= JE; JR++) {
                WORK[(JW + 2) * N + JR] = XSCALE * WORK[(JW + 2) * N + JR];
              }
            }
            XMAX *= XSCALE;
          }

          // Compute the contributions of the off-diagonals of
          // column j (and j+1, if 2-by-2 block) of A and B to the
          // sums.

          for (JA = 1; JA <= NA; JA++) {
            if (ILCPLX) {
              CREALA = ACOEF.value * WORK[2 * N + J + JA - 1];
              CIMAGA = ACOEF.value * WORK[3 * N + J + JA - 1];
              CREALB = BCOEFR.value * WORK[2 * N + J + JA - 1] -
                  BCOEFI.value * WORK[3 * N + J + JA - 1];
              CIMAGB = BCOEFI.value * WORK[2 * N + J + JA - 1] +
                  BCOEFR.value * WORK[3 * N + J + JA - 1];
              for (JR = 1; JR <= J - 1; JR++) {
                WORK[2 * N + JR] -=
                    CREALA * S[JR][J + JA - 1] + CREALB * P[JR][J + JA - 1];
                WORK[3 * N + JR] -=
                    CIMAGA * S[JR][J + JA - 1] + CIMAGB * P[JR][J + JA - 1];
              }
            } else {
              CREALA = ACOEF.value * WORK[2 * N + J + JA - 1];
              CREALB = BCOEFR.value * WORK[2 * N + J + JA - 1];
              for (JR = 1; JR <= J - 1; JR++) {
                WORK[2 * N + JR] -=
                    CREALA * S[JR][J + JA - 1] + CREALB * P[JR][J + JA - 1];
              }
            }
          }
        }

        IL2BY2 = false;
      }

      // Copy eigenvector to VR, back transforming if
      // HOWMNY='B'.

      IEIG -= NW;
      if (ILBACK) {
        for (JW = 0; JW <= NW - 1; JW++) {
          for (JR = 1; JR <= N; JR++) {
            WORK[(JW + 4) * N + JR] = WORK[(JW + 2) * N + 1] * VR[JR][1];
          }

          // A series of compiler directives to defeat
          // vectorization for the next loop

          for (JC = 2; JC <= JE; JC++) {
            for (JR = 1; JR <= N; JR++) {
              WORK[(JW + 4) * N + JR] += WORK[(JW + 2) * N + JC] * VR[JR][JC];
            }
          }
        }

        for (JW = 0; JW <= NW - 1; JW++) {
          for (JR = 1; JR <= N; JR++) {
            VR[JR][IEIG + JW] = WORK[(JW + 4) * N + JR];
          }
        }

        IEND = N;
      } else {
        for (JW = 0; JW <= NW - 1; JW++) {
          for (JR = 1; JR <= N; JR++) {
            VR[JR][IEIG + JW] = WORK[(JW + 2) * N + JR];
          }
        }

        IEND = JE;
      }

      // Scale eigenvector

      XMAX = ZERO;
      if (ILCPLX) {
        for (J = 1; J <= IEND; J++) {
          XMAX = max(XMAX, VR[J][IEIG].abs() + VR[J][IEIG + 1].abs());
        }
      } else {
        for (J = 1; J <= IEND; J++) {
          XMAX = max(XMAX, VR[J][IEIG].abs());
        }
      }

      if (XMAX > SAFMIN) {
        XSCALE = ONE / XMAX;
        for (JW = 0; JW <= NW - 1; JW++) {
          for (JR = 1; JR <= IEND; JR++) {
            VR[JR][IEIG + JW] = XSCALE * VR[JR][IEIG + JW];
          }
        }
      }
    }
  }
}
