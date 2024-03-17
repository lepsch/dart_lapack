import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zladiv.dart';
import 'package:lapack/src/zlanhs.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zrot.dart';

void zhgeqz(
  final String JOB,
  final String COMPQ,
  final String COMPZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<Complex> H_,
  final int LDH,
  final Matrix<Complex> T_,
  final int LDT,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final H = H_.having(ld: LDH);
  final T = T_.having(ld: LDT);
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const HALF = 0.5;
  bool ILAZR2, ILAZRO, ILQ, ILSCHR, ILZ, LQUERY;
  int ICOMPQ,
      ICOMPZ,
      IFIRST = 0,
      IFRSTM,
      IITER,
      ILAST,
      ILASTM,
      IN,
      ISCHUR,
      ISTART = 0,
      J,
      JC,
      JCH,
      JITER,
      JR,
      MAXIT;
  double ABSB,
      ANORM,
      ASCALE,
      ATOL,
      BNORM,
      BSCALE,
      BTOL,
      SAFMIN,
      TEMP,
      TEMP2,
      TEMPR,
      ULP;
  Complex ABI22,
      AD11,
      AD12,
      AD21,
      AD22,
      CTEMP = Complex.zero,
      CTEMP2,
      ESHIFT,
      SHIFT,
      SIGNBC,
      U12,
      X,
      ABI12,
      Y;
  final C = Box(0.0);
  final S = Box(Complex.zero), CTEMP3 = Box(Complex.zero);

  double ABS1(Complex X) => X.toDouble().abs() + X.imaginary.abs();

  // Decode JOB, COMPQ, COMPZ

  if (lsame(JOB, 'E')) {
    ILSCHR = false;
    ISCHUR = 1;
  } else if (lsame(JOB, 'S.value')) {
    ILSCHR = true;
    ISCHUR = 2;
  } else {
    ILSCHR = true;
    ISCHUR = 0;
  }

  if (lsame(COMPQ, 'N')) {
    ILQ = false;
    ICOMPQ = 1;
  } else if (lsame(COMPQ, 'V')) {
    ILQ = true;
    ICOMPQ = 2;
  } else if (lsame(COMPQ, 'I')) {
    ILQ = true;
    ICOMPQ = 3;
  } else {
    ILQ = true;
    ICOMPQ = 0;
  }

  if (lsame(COMPZ, 'N')) {
    ILZ = false;
    ICOMPZ = 1;
  } else if (lsame(COMPZ, 'V')) {
    ILZ = true;
    ICOMPZ = 2;
  } else if (lsame(COMPZ, 'I')) {
    ILZ = true;
    ICOMPZ = 3;
  } else {
    ILZ = true;
    ICOMPZ = 0;
  }

  // Check Argument Values

  INFO.value = 0;
  WORK[1] = max(1, N).toComplex();
  LQUERY = (LWORK == -1);
  if (ISCHUR == 0) {
    INFO.value = -1;
  } else if (ICOMPQ == 0) {
    INFO.value = -2;
  } else if (ICOMPZ == 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (ILO < 1) {
    INFO.value = -5;
  } else if (IHI > N || IHI < ILO - 1) {
    INFO.value = -6;
  } else if (LDH < N) {
    INFO.value = -8;
  } else if (LDT < N) {
    INFO.value = -10;
  } else if (LDQ < 1 || (ILQ && LDQ < N)) {
    INFO.value = -14;
  } else if (LDZ < 1 || (ILZ && LDZ < N)) {
    INFO.value = -16;
  } else if (LWORK < max(1, N) && !LQUERY) {
    INFO.value = -18;
  }
  if (INFO.value != 0) {
    xerbla('ZHGEQZ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  // WORK( 1 ) = CMPLX( 1 )
  if (N <= 0) {
    WORK[1] = Complex.one;
    return;
  }

  // Initialize Q and Z

  if (ICOMPQ == 3) zlaset('Full', N, N, Complex.zero, Complex.one, Q, LDQ);
  if (ICOMPZ == 3) zlaset('Full', N, N, Complex.zero, Complex.one, Z, LDZ);

  // Machine Constants

  IN = IHI + 1 - ILO;
  SAFMIN = dlamch('S.value');
  ULP = dlamch('E') * dlamch('B');
  ANORM = zlanhs('F', IN, H(ILO, ILO), LDH, RWORK);
  BNORM = zlanhs('F', IN, T(ILO, ILO), LDT, RWORK);
  ATOL = max(SAFMIN, ULP * ANORM);
  BTOL = max(SAFMIN, ULP * BNORM);
  ASCALE = ONE / max(SAFMIN, ANORM);
  BSCALE = ONE / max(SAFMIN, BNORM);

  // Set Eigenvalues IHI+1:N

  for (J = IHI + 1; J <= N; J++) {
    ABSB = T[J][J].abs();
    if (ABSB > SAFMIN) {
      SIGNBC = (T[J][J] / ABSB.toComplex()).conjugate();
      T[J][J] = ABSB.toComplex();
      if (ILSCHR) {
        zscal(J - 1, SIGNBC, T(1, J).asArray(), 1);
        zscal(J, SIGNBC, H(1, J).asArray(), 1);
      } else {
        zscal(1, SIGNBC, H(J, J).asArray(), 1);
      }
      if (ILZ) zscal(N, SIGNBC, Z(1, J).asArray(), 1);
    } else {
      T[J][J] = Complex.zero;
    }
    ALPHA[J] = H[J][J];
    BETA[J] = T[J][J];
  }

  // If IHI < ILO, skip QZ steps

  if (IHI >= ILO) {
    // MAIN QZ ITERATION LOOP

    // Initialize dynamic indices

    // Eigenvalues ILAST+1:N have been found.
    //    Column operations modify rows IFRSTM:whatever
    //    Row operations modify columns whatever:ILASTM

    // If only eigenvalues are being computed, then
    //    IFRSTM is the row of the last splitting row above row ILAST;
    //    this is always at least ILO.
    // IITER counts iterations since the last eigenvalue was found,
    //    to tell when to use an extraordinary shift.
    // MAXIT is the maximum number of QZ sweeps allowed.

    ILAST = IHI;
    if (ILSCHR) {
      IFRSTM = 1;
      ILASTM = N;
    } else {
      IFRSTM = ILO;
      ILASTM = IHI;
    }
    IITER = 0;
    ESHIFT = Complex.zero;
    MAXIT = 30 * (IHI - ILO + 1);
    var exhausted = true;
    jiterLoop:
    for (JITER = 1; JITER <= MAXIT; JITER++) {
      // Check for too many iterations.

      if (JITER > MAXIT) break;

      // Split the matrix if possible.

      // Two tests:
      //    1: H(j,j-1)=0  or  j=ILO
      //    2: T(j,j)=0

      // Special case: j=ILAST

      var standardizeB = false;

      if (ILAST == ILO) {
        standardizeB = true;
      } else {
        if (ABS1(H[ILAST][ILAST - 1]) <=
            max(
                SAFMIN,
                ULP *
                    (ABS1(H[ILAST][ILAST]) + ABS1(H[ILAST - 1][ILAST - 1])))) {
          H[ILAST][ILAST - 1] = Complex.zero;
          standardizeB = true;
        }
      }

      if (!standardizeB) {
        var splitOff = false;

        if ((T[ILAST][ILAST]).abs() <= BTOL) {
          T[ILAST][ILAST] = Complex.zero;
        } else {
          // General case: j<ILAST
          var exhausted = true;
          dropThrough:
          for (J = ILAST - 1; J >= ILO; J--) {
            // Test 1: for H(j,j-1)=0 or j=ILO

            if (J == ILO) {
              ILAZRO = true;
            } else {
              if (ABS1(H[J][J - 1]) <=
                  max(SAFMIN, ULP * (ABS1(H[J][J]) + ABS1(H[J - 1][J - 1])))) {
                H[J][J - 1] = Complex.zero;
                ILAZRO = true;
              } else {
                ILAZRO = false;
              }
            }

            // Test 2: for T(j,j)=0

            if ((T[J][J]).abs() < BTOL) {
              T[J][J] = Complex.zero;

              // Test 1a: Check for 2 consecutive small subdiagonals in A

              ILAZR2 = false;
              if (!ILAZRO) {
                if (ABS1(H[J][J - 1]) * (ASCALE * ABS1(H[J + 1][J])) <=
                    ABS1(H[J][J]) * (ASCALE * ATOL)) ILAZR2 = true;
              }

              // If both tests pass (1 & 2), i.e., the leading diagonal
              // element of B in the block is zero, split a 1x1 block off
              // at the top. (I.e., at the J-th row/column) The leading
              // diagonal element of the remainder can also be zero, so
              // this may have to be done repeatedly.

              if (ILAZRO || ILAZR2) {
                for (JCH = J; JCH <= ILAST - 1; JCH++) {
                  CTEMP = H[JCH][JCH];
                  zlartg(CTEMP, H[JCH + 1][JCH], C, S, H.box(JCH, JCH));
                  H[JCH + 1][JCH] = Complex.zero;
                  zrot(ILASTM - JCH, H(JCH, JCH + 1).asArray(), LDH,
                      H(JCH + 1, JCH + 1).asArray(), LDH, C.value, S.value);
                  zrot(ILASTM - JCH, T(JCH, JCH + 1).asArray(), LDT,
                      T(JCH + 1, JCH + 1).asArray(), LDT, C.value, S.value);
                  if (ILQ) {
                    zrot(N, Q(1, JCH).asArray(), 1, Q(1, JCH + 1).asArray(), 1,
                        C.value, S.value.conjugate());
                  }
                  if (ILAZR2) {
                    H[JCH][JCH - 1] *= C.value.toComplex();
                  }
                  ILAZR2 = false;
                  if (ABS1(T[JCH + 1][JCH + 1]) >= BTOL) {
                    if (JCH + 1 >= ILAST) {
                      standardizeB = true;
                      exhausted = false;
                      break dropThrough;
                    } else {
                      IFIRST = JCH + 1;
                      exhausted = false;
                      break dropThrough;
                    }
                  }
                  T[JCH + 1][JCH + 1] = Complex.zero;
                }
                splitOff = true;
                exhausted = false;
                break dropThrough;
              } else {
                // Only test 2 passed -- chase the zero to T(ILAST,ILAST)
                // Then process as in the case T(ILAST,ILAST)=0

                for (JCH = J; JCH <= ILAST - 1; JCH++) {
                  CTEMP = T[JCH][JCH + 1];
                  zlartg(CTEMP, T[JCH + 1][JCH + 1], C, S, T.box(JCH, JCH + 1));
                  T[JCH + 1][JCH + 1] = Complex.zero;
                  if (JCH < ILASTM - 1) {
                    zrot(ILASTM - JCH - 1, T(JCH, JCH + 2).asArray(), LDT,
                        T(JCH + 1, JCH + 2).asArray(), LDT, C.value, S.value);
                  }
                  zrot(ILASTM - JCH + 2, H(JCH, JCH - 1).asArray(), LDH,
                      H(JCH + 1, JCH - 1).asArray(), LDH, C.value, S.value);
                  if (ILQ) {
                    zrot(N, Q(1, JCH).asArray(), 1, Q(1, JCH + 1).asArray(), 1,
                        C.value, S.value.conjugate());
                  }
                  CTEMP = H[JCH + 1][JCH];
                  zlartg(CTEMP, H[JCH + 1][JCH - 1], C, S, H.box(JCH + 1, JCH));
                  H[JCH + 1][JCH - 1] = Complex.zero;
                  zrot(JCH + 1 - IFRSTM, H(IFRSTM, JCH).asArray(), 1,
                      H(IFRSTM, JCH - 1).asArray(), 1, C.value, S.value);
                  zrot(JCH - IFRSTM, T(IFRSTM, JCH).asArray(), 1,
                      T(IFRSTM, JCH - 1).asArray(), 1, C.value, S.value);
                  if (ILZ) {
                    zrot(N, Z(1, JCH).asArray(), 1, Z(1, JCH - 1).asArray(), 1,
                        C.value, S.value);
                  }
                }
                splitOff = true;
                exhausted = false;
                break dropThrough;
              }
            } else if (ILAZRO) {
              // Only test 1 passed -- work on J:ILAST

              IFIRST = J;
              exhausted = false;
              break dropThrough;
            }

            // Neither test passed -- try next J
          }

          if (exhausted) {
            // (Drop-through is "impossible")

            INFO.value = 2 * N + 1;
            WORK[1] = Complex(N.toDouble());
            return;
          }
        }

        // T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
        // 1x1 block.

        if (splitOff) {
          CTEMP = H[ILAST][ILAST];
          zlartg(CTEMP, H[ILAST][ILAST - 1], C, S, H.box(ILAST, ILAST));
          H[ILAST][ILAST - 1] = Complex.zero;
          zrot(ILAST - IFRSTM, H(IFRSTM, ILAST).asArray(), 1,
              H(IFRSTM, ILAST - 1).asArray(), 1, C.value, S.value);
          zrot(ILAST - IFRSTM, T(IFRSTM, ILAST).asArray(), 1,
              T(IFRSTM, ILAST - 1).asArray(), 1, C.value, S.value);
          if (ILZ) {
            zrot(N, Z(1, ILAST).asArray(), 1, Z(1, ILAST - 1).asArray(), 1,
                C.value, S.value);
          }
        }
      }

      // H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA

      if (standardizeB) {
        ABSB = T[ILAST][ILAST].abs();
        if (ABSB > SAFMIN) {
          SIGNBC = (T[ILAST][ILAST] / ABSB.toComplex()).conjugate();
          T[ILAST][ILAST] = ABSB.toComplex();
          if (ILSCHR) {
            zscal(ILAST - IFRSTM, SIGNBC, T(IFRSTM, ILAST).asArray(), 1);
            zscal(ILAST + 1 - IFRSTM, SIGNBC, H(IFRSTM, ILAST).asArray(), 1);
          } else {
            zscal(1, SIGNBC, H(ILAST, ILAST).asArray(), 1);
          }
          if (ILZ) zscal(N, SIGNBC, Z(1, ILAST).asArray(), 1);
        } else {
          T[ILAST][ILAST] = Complex.zero;
        }
        ALPHA[ILAST] = H[ILAST][ILAST];
        BETA[ILAST] = T[ILAST][ILAST];

        // Go to next block -- exit if finished.

        ILAST--;
        if (ILAST < ILO) {
          exhausted = false;
          break jiterLoop;
        }

        // Reset counters

        IITER = 0;
        ESHIFT = Complex.zero;
        if (!ILSCHR) {
          ILASTM = ILAST;
          if (IFRSTM > ILAST) IFRSTM = ILO;
        }
        break;
      }

      // QZ step

      // This iteration only involves rows/columns IFIRST:ILAST.  We
      // assume IFIRST < ILAST, and that the diagonal of B is non-zero.
      IITER++;
      if (!ILSCHR) {
        IFRSTM = IFIRST;
      }

      // Compute the Shift.

      // At this point, IFIRST < ILAST, and the diagonal elements of
      // T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
      // magnitude)

      if ((IITER ~/ 10) * 10 != IITER) {
        // The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
        // the bottom-right 2x2 block of A inv(B) which is nearest to
        // the bottom-right element.

        // We factor B as U*D, where U has unit diagonals, and
        // compute (A*inv(D))*inv(U).

        U12 = (BSCALE.toComplex() * T[ILAST - 1][ILAST]) /
            (BSCALE.toComplex() * T[ILAST][ILAST]);
        AD11 = (ASCALE.toComplex() * H[ILAST - 1][ILAST - 1]) /
            (BSCALE.toComplex() * T[ILAST - 1][ILAST - 1]);
        AD21 = (ASCALE.toComplex() * H[ILAST][ILAST - 1]) /
            (BSCALE.toComplex() * T[ILAST - 1][ILAST - 1]);
        AD12 = (ASCALE.toComplex() * H[ILAST - 1][ILAST]) /
            (BSCALE.toComplex() * T[ILAST][ILAST]);
        AD22 = (ASCALE.toComplex() * H[ILAST][ILAST]) /
            (BSCALE.toComplex() * T[ILAST][ILAST]);
        ABI22 = AD22 - U12 * AD21;
        ABI12 = AD12 - U12 * AD11;

        SHIFT = ABI22;
        CTEMP = ABI12.sqrt() * AD21.sqrt();
        TEMP = ABS1(CTEMP);
        if (CTEMP != Complex.zero) {
          X = HALF.toComplex() * (AD11 - SHIFT);
          TEMP2 = ABS1(X);
          TEMP = max(TEMP, ABS1(X));
          Y = TEMP.toComplex() *
              ((X / TEMP.toComplex()).pow(2) +
                      (CTEMP / TEMP.toComplex()).pow(2))
                  .sqrt();
          if (TEMP2 > ZERO) {
            if ((X / TEMP2.toComplex()).toDouble() * Y.toDouble() +
                    (X / TEMP2.toComplex()).imaginary * Y.imaginary <
                ZERO) Y = -Y;
          }
          SHIFT -= CTEMP * zladiv(CTEMP, (X + Y));
        }
      } else {
        // Exceptional shift.  Chosen for no particularly good reason.

        if ((IITER / 20) * 20 == IITER &&
            BSCALE * ABS1(T[ILAST][ILAST]) > SAFMIN) {
          ESHIFT += (ASCALE.toComplex() * H[ILAST][ILAST]) /
              (BSCALE.toComplex() * T[ILAST][ILAST]);
        } else {
          ESHIFT += (ASCALE.toComplex() * H[ILAST][ILAST - 1]) /
              (BSCALE.toComplex() * T[ILAST - 1][ILAST - 1]);
        }
        SHIFT = ESHIFT;
      }

      // Now check for two consecutive small subdiagonals.
      var consecutive = false;
      for (J = ILAST - 1; J >= IFIRST + 1; J--) {
        ISTART = J;
        CTEMP = ASCALE.toComplex() * H[J][J] -
            SHIFT * (BSCALE.toComplex() * T[J][J]);
        TEMP = ABS1(CTEMP);
        TEMP2 = ASCALE * ABS1(H[J + 1][J]);
        TEMPR = max(TEMP, TEMP2);
        if (TEMPR < ONE && TEMPR != ZERO) {
          TEMP /= TEMPR;
          TEMP2 /= TEMPR;
        }
        if (ABS1(H[J][J - 1]) * TEMP2 <= TEMP * ATOL) {
          consecutive = true;
          break;
        }
      }
      if (!consecutive) {
        ISTART = IFIRST;
        CTEMP = ASCALE.toComplex() * H[IFIRST][IFIRST] -
            SHIFT * (BSCALE.toComplex() * T[IFIRST][IFIRST]);
      }

      // Do an implicit-shift QZ sweep.

      // Initial Q

      CTEMP2 = ASCALE.toComplex() * H[ISTART + 1][ISTART];
      zlartg(CTEMP, CTEMP2, C, S, CTEMP3);

      // Sweep

      for (J = ISTART; J <= ILAST - 1; J++) {
        if (J > ISTART) {
          CTEMP = H[J][J - 1];
          zlartg(CTEMP, H[J + 1][J - 1], C, S, H.box(J, J - 1));
          H[J + 1][J - 1] = Complex.zero;
        }

        for (JC = J; JC <= ILASTM; JC++) {
          CTEMP = C.value.toComplex() * H[J][JC] + S.value * H[J + 1][JC];
          H[J + 1][JC] = -S.value.conjugate() * H[J][JC] +
              C.value.toComplex() * H[J + 1][JC];
          H[J][JC] = CTEMP;
          CTEMP2 = C.value.toComplex() * T[J][JC] + S.value * T[J + 1][JC];
          T[J + 1][JC] = -S.value.conjugate() * T[J][JC] +
              C.value.toComplex() * T[J + 1][JC];
          T[J][JC] = CTEMP2;
        }
        if (ILQ) {
          for (JR = 1; JR <= N; JR++) {
            CTEMP = C.value.toComplex() * Q[JR][J] +
                S.value.conjugate() * Q[JR][J + 1];
            Q[JR][J + 1] =
                -S.value * Q[JR][J] + C.value.toComplex() * Q[JR][J + 1];
            Q[JR][J] = CTEMP;
          }
        }

        CTEMP = T[J + 1][J + 1];
        zlartg(CTEMP, T[J + 1][J], C, S, T.box(J + 1, J + 1));
        T[J + 1][J] = Complex.zero;

        for (JR = IFRSTM; JR <= min(J + 2, ILAST); JR++) {
          CTEMP = C.value.toComplex() * H[JR][J + 1] + S.value * H[JR][J];
          H[JR][J] = -S.value.conjugate() * H[JR][J + 1] +
              C.value.toComplex() * H[JR][J];
          H[JR][J + 1] = CTEMP;
        }
        for (JR = IFRSTM; JR <= J; JR++) {
          CTEMP = C.value.toComplex() * T[JR][J + 1] + S.value * T[JR][J];
          T[JR][J] = -S.value.conjugate() * T[JR][J + 1] +
              C.value.toComplex() * T[JR][J];
          T[JR][J + 1] = CTEMP;
        }
        if (ILZ) {
          for (JR = 1; JR <= N; JR++) {
            CTEMP = C.value.toComplex() * Z[JR][J + 1] + S.value * Z[JR][J];
            Z[JR][J] = -S.value.conjugate() * Z[JR][J + 1] +
                C.value.toComplex() * Z[JR][J];
            Z[JR][J + 1] = CTEMP;
          }
        }
      }
    }

    // Drop-through = non-convergence

    if (exhausted) {
      INFO.value = ILAST;
      WORK[1] = Complex(N.toDouble());
      return;
    }
  }

  // Successful completion of all QZ steps

  // Set Eigenvalues 1:ILO-1

  for (J = 1; J <= ILO - 1; J++) {
    ABSB = (T[J][J]).abs();
    if (ABSB > SAFMIN) {
      SIGNBC = T[J][J].conjugate() / ABSB.toComplex();
      T[J][J] = ABSB.toComplex();
      if (ILSCHR) {
        zscal(J - 1, SIGNBC, T(1, J).asArray(), 1);
        zscal(J, SIGNBC, H(1, J).asArray(), 1);
      } else {
        zscal(1, SIGNBC, H(J, J).asArray(), 1);
      }
      if (ILZ) zscal(N, SIGNBC, Z(1, J).asArray(), 1);
    } else {
      T[J][J] = Complex.zero;
    }
    ALPHA[J] = H[J][J];
    BETA[J] = T[J][J];
  }

  // Normal Termination

  INFO.value = 0;

  // Exit (other than argument error) -- return optimal workspace size

  WORK[1] = Complex(N.toDouble());
}
