// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/ddot.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlaln2.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtrevc3(
  final String SIDE,
  final String HOWMNY,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> VL_,
  final int LDVL,
  final Matrix<double> VR_,
  final int LDVR,
  final int MM,
  final Box<int> M,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final SELECT = SELECT_.having();
  final T = T_.having(ld: LDT);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const NBMIN = 8, NBMAX = 128;
  bool ALLV, BOTHV, LEFTV, LQUERY, OVER, PAIR, RIGHTV, SOMEV;
  int I, II, IP, IS, J, J1, J2, JNXT, K, KI, IV, MAXWRK, NB, KI2;
  double BETA,
      BIGNUM,
      EMAX,
      REC,
      REMAX = 0,
      SMIN,
      SMLNUM,
      ULP,
      UNFL,
      VCRIT,
      VMAX,
      WI,
      WR;
  final X = Matrix<double>(2, 2);
  final ISCOMPLEX = Array<int>(NBMAX);
  final IERR = Box(0);
  final SCALE = Box(0.0), XNORM = Box(0.0);

  // Decode and test the input parameters
  BOTHV = lsame(SIDE, 'B');
  RIGHTV = lsame(SIDE, 'R') || BOTHV;
  LEFTV = lsame(SIDE, 'L') || BOTHV;

  ALLV = lsame(HOWMNY, 'A');
  OVER = lsame(HOWMNY, 'B');
  SOMEV = lsame(HOWMNY, 'S');

  INFO.value = 0;
  NB = ilaenv(1, 'DTREVC', SIDE + HOWMNY, N, -1, -1, -1);
  MAXWRK = max(1, N + 2 * N * NB);
  WORK[1] = MAXWRK.toDouble();
  LQUERY = LWORK == -1;
  if (!RIGHTV && !LEFTV) {
    INFO.value = -1;
  } else if (!ALLV && !OVER && !SOMEV) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDT < max(1, N)) {
    INFO.value = -6;
  } else if (LDVL < 1 || (LEFTV && LDVL < N)) {
    INFO.value = -8;
  } else if (LDVR < 1 || (RIGHTV && LDVR < N)) {
    INFO.value = -10;
  } else if (LWORK < max(1, 3 * N) && !LQUERY) {
    INFO.value = -14;
  } else {
    // Set M to the number of columns required to store the selected
    // eigenvectors, standardize the array SELECT if necessary, and
    // test MM.
    if (SOMEV) {
      M.value = 0;
      PAIR = false;
      for (J = 1; J <= N; J++) {
        if (PAIR) {
          PAIR = false;
          SELECT[J] = false;
        } else {
          if (J < N) {
            if (T[J + 1][J] == ZERO) {
              if (SELECT[J]) M.value++;
            } else {
              PAIR = true;
              if (SELECT[J] || SELECT[J + 1]) {
                SELECT[J] = true;
                M.value += 2;
              }
            }
          } else {
            if (SELECT[N]) M.value++;
          }
        }
      }
    } else {
      M.value = N;
    }

    if (MM < M.value) {
      INFO.value = -11;
    }
  }
  if (INFO.value != 0) {
    xerbla('DTREVC3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible.
  if (N == 0) return;

  // Use blocked version of back-transformation if sufficient workspace.
  // Zero-out the workspace to avoid potential NaN propagation.
  if (OVER && LWORK >= N + 2 * N * NBMIN) {
    NB = (LWORK - N) ~/ (2 * N);
    NB = min(NB, NBMAX);
    dlaset('F', N, 1 + 2 * NB, ZERO, ZERO, WORK.asMatrix(N), N);
  } else {
    NB = 1;
  }

  // Set the constants to control overflow.
  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');
  SMLNUM = UNFL * (N / ULP);
  BIGNUM = (ONE - ULP) / SMLNUM;

  // Compute 1-norm of each column of strictly upper triangular
  // part of T to control overflow in triangular solver.
  WORK[1] = ZERO;
  for (J = 2; J <= N; J++) {
    WORK[J] = ZERO;
    for (I = 1; I <= J - 1; I++) {
      WORK[J] += T[I][J].abs();
    }
  }

  // Index IP is used to specify the real or complex eigenvalue:
  // IP = 0, real eigenvalue,
  // 1, first  of conjugate complex pair: (wr,wi)
  // -1, second of conjugate complex pair: (wr,wi)
  // ISCOMPLEX array stores IP for each column in current block.
  if (RIGHTV) {
    // Compute right eigenvectors.

    // IV is index of column in current block.
    // For complex right vector, uses IV-1 for real part and IV for complex part.
    // Non-blocked version always uses IV=2;
    // blocked     version starts with IV=NB, goes down to 1 or 2.
    // (Note the "0-th" column is used for 1-norms computed above.)
    IV = 2;
    if (NB > 2) {
      IV = NB;
    }

    IP = 0;
    IS = M.value;
    for (KI = N; KI >= 1; KI--) {
      if (IP == -1) {
        // previous iteration (ki+1) was second of conjugate pair,
        // so this ki is first of conjugate pair; skip to end of loop
        IP = 1;
        continue;
      } else if (KI == 1) {
        // last column, so this ki must be real eigenvalue
        IP = 0;
      } else if (T[KI][KI - 1] == ZERO) {
        // zero on sub-diagonal, so this ki is real eigenvalue
        IP = 0;
      } else {
        // non-zero on sub-diagonal, so this ki is second of conjugate pair
        IP = -1;
      }

      if (SOMEV) {
        if (IP == 0) {
          if (!SELECT[KI]) continue;
        } else {
          if (!SELECT[KI - 1]) continue;
        }
      }

      // Compute the KI-th eigenvalue (WR,WI).
      WR = T[KI][KI];
      WI = ZERO;
      if (IP != 0) {
        WI = sqrt(T[KI][KI - 1].abs()) * sqrt(T[KI - 1][KI].abs());
      }
      SMIN = max(ULP * (WR.abs() + WI.abs()), SMLNUM);

      if (IP == 0) {
        // Real right eigenvector
        WORK[KI + IV * N] = ONE;

        // Form right-hand side.
        for (K = 1; K <= KI - 1; K++) {
          WORK[K + IV * N] = -T[K][KI];
        }

        // Solve upper quasi-triangular system:
        // [ T[1:KI-1][1:KI-1] - WR ]*X = SCALE*WORK.
        JNXT = KI - 1;
        for (J = KI - 1; J >= 1; J--) {
          if (J > JNXT) continue;
          J1 = J;
          J2 = J;
          JNXT = J - 1;
          if (J > 1) {
            if (T[J][J - 1] != ZERO) {
              J1 = J - 1;
              JNXT = J - 2;
            }
          }

          if (J1 == J2) {
            // 1-by-1 diagonal block
            dlaln2(
                false,
                1,
                1,
                SMIN,
                ONE,
                T(J, J),
                LDT,
                ONE,
                ONE,
                WORK(J + IV * N).asMatrix(N),
                N,
                WR,
                ZERO,
                X,
                2,
                SCALE,
                XNORM,
                IERR);

            // Scale X[1][1] to avoid overflow when updating
            // the right-hand side.
            if (XNORM.value > ONE) {
              if (WORK[J] > BIGNUM / XNORM.value) {
                X[1][1] /= XNORM.value;
                SCALE.value /= XNORM.value;
              }
            }

            // Scale if necessary
            if (SCALE.value != ONE) dscal(KI, SCALE.value, WORK(1 + IV * N), 1);
            WORK[J + IV * N] = X[1][1];

            // Update right-hand side
            daxpy(J - 1, -X[1][1], T(1, J).asArray(), 1, WORK(1 + IV * N), 1);
          } else {
            // 2-by-2 diagonal block
            dlaln2(
                false,
                2,
                1,
                SMIN,
                ONE,
                T(J - 1, J - 1),
                LDT,
                ONE,
                ONE,
                WORK(J - 1 + IV * N).asMatrix(N),
                N,
                WR,
                ZERO,
                X,
                2,
                SCALE,
                XNORM,
                IERR);

            // Scale X[1][1] and X[2][1] to avoid overflow when
            // updating the right-hand side.
            if (XNORM.value > ONE) {
              BETA = max(WORK[J - 1], WORK[J]);
              if (BETA > BIGNUM / XNORM.value) {
                X[1][1] /= XNORM.value;
                X[2][1] /= XNORM.value;
                SCALE.value /= XNORM.value;
              }
            }

            // Scale if necessary
            if (SCALE.value != ONE) dscal(KI, SCALE.value, WORK(1 + IV * N), 1);
            WORK[J - 1 + IV * N] = X[1][1];
            WORK[J + IV * N] = X[2][1];

            // Update right-hand side
            daxpy(
                J - 2, -X[1][1], T(1, J - 1).asArray(), 1, WORK(1 + IV * N), 1);
            daxpy(J - 2, -X[2][1], T(1, J).asArray(), 1, WORK(1 + IV * N), 1);
          }
        }

        // Copy the vector x or Q*x to VR and normalize.
        if (!OVER) {
          // no back-transform: copy x to VR and normalize.
          dcopy(KI, WORK(1 + IV * N), 1, VR(1, IS).asArray(), 1);

          II = idamax(KI, VR(1, IS).asArray(), 1);
          REMAX = ONE / VR[II][IS].abs();
          dscal(KI, REMAX, VR(1, IS).asArray(), 1);

          for (K = KI + 1; K <= N; K++) {
            VR[K][IS] = ZERO;
          }
        } else if (NB == 1) {
          // version 1: back-transform each vector with GEMV, Q*x.
          if (KI > 1) {
            dgemv('N', N, KI - 1, ONE, VR, LDVR, WORK(1 + IV * N), 1,
                WORK[KI + IV * N], VR(1, KI).asArray(), 1);
          }

          II = idamax(N, VR(1, KI).asArray(), 1);
          REMAX = ONE / VR[II][KI].abs();
          dscal(N, REMAX, VR(1, KI).asArray(), 1);
        } else {
          // version 2: back-transform block of vectors with GEMM
          // zero out below vector
          for (K = KI + 1; K <= N; K++) {
            WORK[K + IV * N] = ZERO;
          }
          ISCOMPLEX[IV] = IP;
          // back-transform and normalization is done below
        }
      } else {
        // Complex right eigenvector.

        // Initial solve
        // [ ( T[KI-1][KI-1] T[KI-1][KI] ) - (WR + I*WI) ]*X = 0.
        // [ ( T[KI][  KI-1] T[KI][  KI] )               ]
        if (T[KI - 1][KI].abs() >= T[KI][KI - 1].abs()) {
          WORK[KI - 1 + (IV - 1) * N] = ONE;
          WORK[KI + IV * N] = WI / T[KI - 1][KI];
        } else {
          WORK[KI - 1 + (IV - 1) * N] = -WI / T[KI][KI - 1];
          WORK[KI + IV * N] = ONE;
        }
        WORK[KI + (IV - 1) * N] = ZERO;
        WORK[KI - 1 + IV * N] = ZERO;

        // Form right-hand side.
        for (K = 1; K <= KI - 2; K++) {
          WORK[K + (IV - 1) * N] = -WORK[KI - 1 + (IV - 1) * N] * T[K][KI - 1];
          WORK[K + IV * N] = -WORK[KI + IV * N] * T[K][KI];
        }

        // Solve upper quasi-triangular system:
        // [ T[1:KI-2][1:KI-2] - (WR+i*WI) ]*X = SCALE*(WORK+i*WORK2)
        JNXT = KI - 2;
        for (J = KI - 2; J >= 1; J--) {
          if (J > JNXT) continue;
          J1 = J;
          J2 = J;
          JNXT = J - 1;
          if (J > 1) {
            if (T[J][J - 1] != ZERO) {
              J1 = J - 1;
              JNXT = J - 2;
            }
          }

          if (J1 == J2) {
            // 1-by-1 diagonal block
            dlaln2(
                false,
                1,
                2,
                SMIN,
                ONE,
                T(J, J),
                LDT,
                ONE,
                ONE,
                WORK(J + (IV - 1) * N).asMatrix(N),
                N,
                WR,
                WI,
                X,
                2,
                SCALE,
                XNORM,
                IERR);

            // Scale X[1][1] and X[1][2] to avoid overflow when
            // updating the right-hand side.
            if (XNORM.value > ONE) {
              if (WORK[J] > BIGNUM / XNORM.value) {
                X[1][1] /= XNORM.value;
                X[1][2] /= XNORM.value;
                SCALE.value /= XNORM.value;
              }
            }

            // Scale if necessary
            if (SCALE.value != ONE) {
              dscal(KI, SCALE.value, WORK(1 + (IV - 1) * N), 1);
              dscal(KI, SCALE.value, WORK(1 + IV * N), 1);
            }
            WORK[J + (IV - 1) * N] = X[1][1];
            WORK[J + IV * N] = X[1][2];

            // Update the right-hand side
            daxpy(J - 1, -X[1][1], T(1, J).asArray(), 1, WORK(1 + (IV - 1) * N),
                1);
            daxpy(J - 1, -X[1][2], T(1, J).asArray(), 1, WORK(1 + IV * N), 1);
          } else {
            // 2-by-2 diagonal block
            dlaln2(
                false,
                2,
                2,
                SMIN,
                ONE,
                T(J - 1, J - 1),
                LDT,
                ONE,
                ONE,
                WORK(J - 1 + (IV - 1) * N).asMatrix(N),
                N,
                WR,
                WI,
                X,
                2,
                SCALE,
                XNORM,
                IERR);

            // Scale X to avoid overflow when updating
            // the right-hand side.
            if (XNORM.value > ONE) {
              BETA = max(WORK[J - 1], WORK[J]);
              if (BETA > BIGNUM / XNORM.value) {
                REC = ONE / XNORM.value;
                X[1][1] *= REC;
                X[1][2] *= REC;
                X[2][1] *= REC;
                X[2][2] *= REC;
                SCALE.value *= REC;
              }
            }

            // Scale if necessary
            if (SCALE.value != ONE) {
              dscal(KI, SCALE.value, WORK(1 + (IV - 1) * N), 1);
              dscal(KI, SCALE.value, WORK(1 + IV * N), 1);
            }
            WORK[J - 1 + (IV - 1) * N] = X[1][1];
            WORK[J + (IV - 1) * N] = X[2][1];
            WORK[J - 1 + IV * N] = X[1][2];
            WORK[J + IV * N] = X[2][2];

            // Update the right-hand side
            daxpy(J - 2, -X[1][1], T(1, J - 1).asArray(), 1,
                WORK(1 + (IV - 1) * N), 1);
            daxpy(J - 2, -X[2][1], T(1, J).asArray(), 1, WORK(1 + (IV - 1) * N),
                1);
            daxpy(
                J - 2, -X[1][2], T(1, J - 1).asArray(), 1, WORK(1 + IV * N), 1);
            daxpy(J - 2, -X[2][2], T(1, J).asArray(), 1, WORK(1 + IV * N), 1);
          }
        }

        // Copy the vector x or Q*x to VR and normalize.
        if (!OVER) {
          // no back-transform: copy x to VR and normalize.
          dcopy(KI, WORK(1 + (IV - 1) * N), 1, VR(1, IS - 1).asArray(), 1);
          dcopy(KI, WORK(1 + IV * N), 1, VR(1, IS).asArray(), 1);

          EMAX = ZERO;
          for (K = 1; K <= KI; K++) {
            EMAX = max(EMAX, VR[K][IS - 1].abs() + VR[K][IS].abs());
          }
          REMAX = ONE / EMAX;
          dscal(KI, REMAX, VR(1, IS - 1).asArray(), 1);
          dscal(KI, REMAX, VR(1, IS).asArray(), 1);

          for (K = KI + 1; K <= N; K++) {
            VR[K][IS - 1] = ZERO;
            VR[K][IS] = ZERO;
          }
        } else if (NB == 1) {
          // version 1: back-transform each vector with GEMV, Q*x.
          if (KI > 2) {
            dgemv('N', N, KI - 2, ONE, VR, LDVR, WORK(1 + (IV - 1) * N), 1,
                WORK[KI - 1 + (IV - 1) * N], VR(1, KI - 1).asArray(), 1);
            dgemv('N', N, KI - 2, ONE, VR, LDVR, WORK(1 + IV * N), 1,
                WORK[KI + IV * N], VR(1, KI).asArray(), 1);
          } else {
            dscal(N, WORK[KI - 1 + (IV - 1) * N], VR(1, KI - 1).asArray(), 1);
            dscal(N, WORK[KI + IV * N], VR(1, KI).asArray(), 1);
          }

          EMAX = ZERO;
          for (K = 1; K <= N; K++) {
            EMAX = max(EMAX, VR[K][KI - 1].abs() + VR[K][KI].abs());
          }
          REMAX = ONE / EMAX;
          dscal(N, REMAX, VR(1, KI - 1).asArray(), 1);
          dscal(N, REMAX, VR(1, KI).asArray(), 1);
        } else {
          // version 2: back-transform block of vectors with GEMM
          // zero out below vector
          for (K = KI + 1; K <= N; K++) {
            WORK[K + (IV - 1) * N] = ZERO;
            WORK[K + IV * N] = ZERO;
          }
          ISCOMPLEX[IV - 1] = -IP;
          ISCOMPLEX[IV] = IP;
          IV--;
          // back-transform and normalization is done below
        }
      }

      if (NB > 1) {
        // Blocked version of back-transform
        // For complex case, KI2 includes both vectors (KI-1 and KI)
        if (IP == 0) {
          KI2 = KI;
        } else {
          KI2 = KI - 1;
        }

        // Columns IV:NB of work are valid vectors.
        // When the number of vectors stored reaches NB-1 or NB,
        // or if this was last vector, do the GEMM
        if ((IV <= 2) || (KI2 == 1)) {
          dgemm(
              'N',
              'N',
              N,
              NB - IV + 1,
              KI2 + NB - IV,
              ONE,
              VR,
              LDVR,
              WORK(1 + IV * N).asMatrix(N),
              N,
              ZERO,
              WORK(1 + (NB + IV) * N).asMatrix(N),
              N);
          // normalize vectors
          for (K = IV; K <= NB; K++) {
            if (ISCOMPLEX[K] == 0) {
              // real eigenvector
              II = idamax(N, WORK(1 + (NB + K) * N), 1);
              REMAX = ONE / WORK[II + (NB + K) * N].abs();
            } else if (ISCOMPLEX[K] == 1) {
              // first eigenvector of conjugate pair
              EMAX = ZERO;
              for (II = 1; II <= N; II++) {
                EMAX = max(
                    EMAX,
                    WORK[II + (NB + K) * N].abs() +
                        WORK[II + (NB + K + 1) * N].abs());
              }
              REMAX = ONE / EMAX;
              // else if ISCOMPLEX[K] == -1
              // second eigenvector of conjugate pair
              // reuse same REMAX as previous K
            }
            dscal(N, REMAX, WORK(1 + (NB + K) * N), 1);
          }
          dlacpy('F', N, NB - IV + 1, WORK(1 + (NB + IV) * N).asMatrix(N), N,
              VR(1, KI2), LDVR);
          IV = NB;
        } else {
          IV--;
        }
      } // ! blocked back-transform;

      IS--;
      if (IP != 0) IS--;
    }
  }

  if (LEFTV) {
    // Compute left eigenvectors.

    // IV is index of column in current block.
    // For complex left vector, uses IV for real part and IV+1 for complex part.
    // Non-blocked version always uses IV=1;
    // blocked     version starts with IV=1, goes up to NB-1 or NB.
    // (Note the "0-th" column is used for 1-norms computed above.)
    IV = 1;
    IP = 0;
    IS = 1;
    for (KI = 1; KI <= N; KI++) {
      if (IP == 1) {
        // previous iteration (ki-1) was first of conjugate pair,
        // so this ki is second of conjugate pair; skip to end of loop
        IP = -1;
        continue;
      } else if (KI == N) {
        // last column, so this ki must be real eigenvalue
        IP = 0;
      } else if (T[KI + 1][KI] == ZERO) {
        // zero on sub-diagonal, so this ki is real eigenvalue
        IP = 0;
      } else {
        // non-zero on sub-diagonal, so this ki is first of conjugate pair
        IP = 1;
      }

      if (SOMEV) {
        if (!SELECT[KI]) continue;
      }

      // Compute the KI-th eigenvalue (WR,WI).
      WR = T[KI][KI];
      WI = ZERO;
      if (IP != 0) {
        WI = sqrt(T[KI][KI + 1].abs()) * sqrt(T[KI + 1][KI].abs());
      }
      SMIN = max(ULP * (WR.abs() + WI.abs()), SMLNUM);

      if (IP == 0) {
        // Real left eigenvector
        WORK[KI + IV * N] = ONE;

        // Form right-hand side.
        for (K = KI + 1; K <= N; K++) {
          WORK[K + IV * N] = -T[KI][K];
        }

        // Solve transposed quasi-triangular system:
        // [ T[KI+1:N][KI+1:N] - WR ]**T * X = SCALE*WORK
        VMAX = ONE;
        VCRIT = BIGNUM;

        JNXT = KI + 1;
        for (J = KI + 1; J <= N; J++) {
          if (J < JNXT) continue;
          J1 = J;
          J2 = J;
          JNXT = J + 1;
          if (J < N) {
            if (T[J + 1][J] != ZERO) {
              J2 = J + 1;
              JNXT = J + 2;
            }
          }

          if (J1 == J2) {
            // 1-by-1 diagonal block

            // Scale if necessary to avoid overflow when forming
            // the right-hand side.
            if (WORK[J] > VCRIT) {
              REC = ONE / VMAX;
              dscal(N - KI + 1, REC, WORK(KI + IV * N), 1);
              VMAX = ONE;
              VCRIT = BIGNUM;
            }

            WORK[J + IV * N] -= ddot(J - KI - 1, T(KI + 1, J).asArray(), 1,
                WORK(KI + 1 + IV * N), 1);

            // Solve [ T[J][J] - WR ]**T * X = WORK
            dlaln2(
                false,
                1,
                1,
                SMIN,
                ONE,
                T(J, J),
                LDT,
                ONE,
                ONE,
                WORK(J + IV * N).asMatrix(N),
                N,
                WR,
                ZERO,
                X,
                2,
                SCALE,
                XNORM,
                IERR);

            // Scale if necessary
            if (SCALE.value != ONE) {
              dscal(N - KI + 1, SCALE.value, WORK(KI + IV * N), 1);
            }
            WORK[J + IV * N] = X[1][1];
            VMAX = max(WORK[J + IV * N].abs(), VMAX);
            VCRIT = BIGNUM / VMAX;
          } else {
            // 2-by-2 diagonal block

            // Scale if necessary to avoid overflow when forming
            // the right-hand side.
            BETA = max(WORK[J], WORK[J + 1]);
            if (BETA > VCRIT) {
              REC = ONE / VMAX;
              dscal(N - KI + 1, REC, WORK(KI + IV * N), 1);
              VMAX = ONE;
              VCRIT = BIGNUM;
            }

            WORK[J + IV * N] -= ddot(J - KI - 1, T(KI + 1, J).asArray(), 1,
                WORK(KI + 1 + IV * N), 1);

            WORK[J + 1 + IV * N] -= ddot(J - KI - 1, T(KI + 1, J + 1).asArray(),
                1, WORK(KI + 1 + IV * N), 1);

            // Solve
            // [ T[J][J]-WR   T[J][J+1]      ]**T * X = SCALE*( WORK1 )
            // [ T[J+1][J]    T[J+1][J+1]-WR ]                ( WORK2 )
            dlaln2(
                true,
                2,
                1,
                SMIN,
                ONE,
                T(J, J),
                LDT,
                ONE,
                ONE,
                WORK(J + IV * N).asMatrix(N),
                N,
                WR,
                ZERO,
                X,
                2,
                SCALE,
                XNORM,
                IERR);

            // Scale if necessary
            if (SCALE.value != ONE) {
              dscal(N - KI + 1, SCALE.value, WORK(KI + IV * N), 1);
            }
            WORK[J + IV * N] = X[1][1];
            WORK[J + 1 + IV * N] = X[2][1];

            VMAX = max(
                WORK[J + IV * N].abs(), max(WORK[J + 1 + IV * N].abs(), VMAX));
            VCRIT = BIGNUM / VMAX;
          }
        }

        // Copy the vector x or Q*x to VL and normalize.
        if (!OVER) {
          // no back-transform: copy x to VL and normalize.
          dcopy(N - KI + 1, WORK(KI + IV * N), 1, VL(KI, IS).asArray(), 1);

          II = idamax(N - KI + 1, VL(KI, IS).asArray(), 1) + KI - 1;
          REMAX = ONE / VL[II][IS].abs();
          dscal(N - KI + 1, REMAX, VL(KI, IS).asArray(), 1);

          for (K = 1; K <= KI - 1; K++) {
            VL[K][IS] = ZERO;
          }
        } else if (NB == 1) {
          // version 1: back-transform each vector with GEMV, Q*x.
          if (KI < N) {
            dgemv(
                'N',
                N,
                N - KI,
                ONE,
                VL(1, KI + 1),
                LDVL,
                WORK(KI + 1 + IV * N),
                1,
                WORK[KI + IV * N],
                VL(1, KI).asArray(),
                1);
          }

          II = idamax(N, VL(1, KI).asArray(), 1);
          REMAX = ONE / VL[II][KI].abs();
          dscal(N, REMAX, VL(1, KI).asArray(), 1);
        } else {
          // version 2: back-transform block of vectors with GEMM
          // zero out above vector
          // could go from KI-NV+1 to KI-1
          for (K = 1; K <= KI - 1; K++) {
            WORK[K + IV * N] = ZERO;
          }
          ISCOMPLEX[IV] = IP;
          // back-transform and normalization is done below
        }
      } else {
        // Complex left eigenvector.

        // Initial solve:
        // [ ( T[KI][KI]    T[KI][KI+1]  )**T - (WR - I* WI) ]*X = 0.
        // [ ( T[KI+1][KI] T[KI+1][KI+1] )                   ]
        if (T[KI][KI + 1].abs() >= T[KI + 1][KI].abs()) {
          WORK[KI + IV * N] = WI / T[KI][KI + 1];
          WORK[KI + 1 + (IV + 1) * N] = ONE;
        } else {
          WORK[KI + IV * N] = ONE;
          WORK[KI + 1 + (IV + 1) * N] = -WI / T[KI + 1][KI];
        }
        WORK[KI + 1 + IV * N] = ZERO;
        WORK[KI + (IV + 1) * N] = ZERO;

        // Form right-hand side.
        for (K = KI + 2; K <= N; K++) {
          WORK[K + IV * N] = -WORK[KI + IV * N] * T[KI][K];
          WORK[K + (IV + 1) * N] = -WORK[KI + 1 + (IV + 1) * N] * T[KI + 1][K];
        }

        // Solve transposed quasi-triangular system:
        // [ T[KI+2:N][KI+2:N]**T - (WR-i*WI) ]*X = WORK1+i*WORK2
        VMAX = ONE;
        VCRIT = BIGNUM;

        JNXT = KI + 2;
        for (J = KI + 2; J <= N; J++) {
          if (J < JNXT) continue;
          J1 = J;
          J2 = J;
          JNXT = J + 1;
          if (J < N) {
            if (T[J + 1][J] != ZERO) {
              J2 = J + 1;
              JNXT = J + 2;
            }
          }

          if (J1 == J2) {
            // 1-by-1 diagonal block

            // Scale if necessary to avoid overflow when
            // forming the right-hand side elements.
            if (WORK[J] > VCRIT) {
              REC = ONE / VMAX;
              dscal(N - KI + 1, REC, WORK(KI + IV * N), 1);
              dscal(N - KI + 1, REC, WORK(KI + (IV + 1) * N), 1);
              VMAX = ONE;
              VCRIT = BIGNUM;
            }

            WORK[J + IV * N] -= ddot(J - KI - 2, T(KI + 2, J).asArray(), 1,
                WORK(KI + 2 + IV * N), 1);
            WORK[J + (IV + 1) * N] -= ddot(J - KI - 2, T(KI + 2, J).asArray(),
                1, WORK(KI + 2 + (IV + 1) * N), 1);

            // Solve [ T[J][J]-(WR-i*WI) ]*(X11+i*X12)= WK+I*WK2
            dlaln2(
                false,
                1,
                2,
                SMIN,
                ONE,
                T(J, J),
                LDT,
                ONE,
                ONE,
                WORK(J + IV * N).asMatrix(N),
                N,
                WR,
                -WI,
                X,
                2,
                SCALE,
                XNORM,
                IERR);

            // Scale if necessary
            if (SCALE.value != ONE) {
              dscal(N - KI + 1, SCALE.value, WORK(KI + IV * N), 1);
              dscal(N - KI + 1, SCALE.value, WORK(KI + (IV + 1) * N), 1);
            }
            WORK[J + IV * N] = X[1][1];
            WORK[J + (IV + 1) * N] = X[1][2];
            VMAX = max(WORK[J + IV * N].abs(),
                max(WORK[J + (IV + 1) * N].abs(), VMAX));
            VCRIT = BIGNUM / VMAX;
          } else {
            // 2-by-2 diagonal block

            // Scale if necessary to avoid overflow when forming
            // the right-hand side elements.
            BETA = max(WORK[J], WORK[J + 1]);
            if (BETA > VCRIT) {
              REC = ONE / VMAX;
              dscal(N - KI + 1, REC, WORK(KI + IV * N), 1);
              dscal(N - KI + 1, REC, WORK(KI + (IV + 1) * N), 1);
              VMAX = ONE;
              VCRIT = BIGNUM;
            }

            WORK[J + IV * N] -= ddot(J - KI - 2, T(KI + 2, J).asArray(), 1,
                WORK(KI + 2 + IV * N), 1);

            WORK[J + (IV + 1) * N] -= ddot(J - KI - 2, T(KI + 2, J).asArray(),
                1, WORK(KI + 2 + (IV + 1) * N), 1);

            WORK[J + 1 + IV * N] -= ddot(J - KI - 2, T(KI + 2, J + 1).asArray(),
                1, WORK(KI + 2 + IV * N), 1);

            WORK[J + 1 + (IV + 1) * N] -= ddot(J - KI - 2,
                T(KI + 2, J + 1).asArray(), 1, WORK(KI + 2 + (IV + 1) * N), 1);

            // Solve 2-by-2 complex linear equation
            // [ (T[j][j]   T[j][j+1]  )**T - (wr-i*wi)*I ]*X = SCALE*B
            // [ (T[j+1][j] T[j+1][j+1])                  ]
            dlaln2(
                true,
                2,
                2,
                SMIN,
                ONE,
                T(J, J),
                LDT,
                ONE,
                ONE,
                WORK(J + IV * N).asMatrix(N),
                N,
                WR,
                -WI,
                X,
                2,
                SCALE,
                XNORM,
                IERR);

            // Scale if necessary
            if (SCALE.value != ONE) {
              dscal(N - KI + 1, SCALE.value, WORK(KI + IV * N), 1);
              dscal(N - KI + 1, SCALE.value, WORK(KI + (IV + 1) * N), 1);
            }
            WORK[J + IV * N] = X[1][1];
            WORK[J + (IV + 1) * N] = X[1][2];
            WORK[J + 1 + IV * N] = X[2][1];
            WORK[J + 1 + (IV + 1) * N] = X[2][2];
            VMAX = max(
                X[1][1].abs(),
                max(
                  max(X[1][2].abs(), X[2][1].abs()),
                  max(X[2][2].abs(), VMAX),
                ));
            VCRIT = BIGNUM / VMAX;
          }
        }

        // Copy the vector x or Q*x to VL and normalize.
        if (!OVER) {
          // no back-transform: copy x to VL and normalize.
          dcopy(N - KI + 1, WORK(KI + IV * N), 1, VL(KI, IS).asArray(), 1);
          dcopy(N - KI + 1, WORK(KI + (IV + 1) * N), 1,
              VL(KI, IS + 1).asArray(), 1);

          EMAX = ZERO;
          for (K = KI; K <= N; K++) {
            EMAX = max(EMAX, VL[K][IS].abs() + VL[K][IS + 1].abs());
          }
          REMAX = ONE / EMAX;
          dscal(N - KI + 1, REMAX, VL(KI, IS).asArray(), 1);
          dscal(N - KI + 1, REMAX, VL(KI, IS + 1).asArray(), 1);

          for (K = 1; K <= KI - 1; K++) {
            VL[K][IS] = ZERO;
            VL[K][IS + 1] = ZERO;
          }
        } else if (NB == 1) {
          // version 1: back-transform each vector with GEMV, Q*x.
          if (KI < N - 1) {
            dgemv(
                'N',
                N,
                N - KI - 1,
                ONE,
                VL(1, KI + 2),
                LDVL,
                WORK(KI + 2 + IV * N),
                1,
                WORK[KI + IV * N],
                VL(1, KI).asArray(),
                1);
            dgemv(
                'N',
                N,
                N - KI - 1,
                ONE,
                VL(1, KI + 2),
                LDVL,
                WORK(KI + 2 + (IV + 1) * N),
                1,
                WORK[KI + 1 + (IV + 1) * N],
                VL(1, KI + 1).asArray(),
                1);
          } else {
            dscal(N, WORK[KI + IV * N], VL(1, KI).asArray(), 1);
            dscal(N, WORK[KI + 1 + (IV + 1) * N], VL(1, KI + 1).asArray(), 1);
          }

          EMAX = ZERO;
          for (K = 1; K <= N; K++) {
            EMAX = max(EMAX, VL[K][KI].abs() + VL[K][KI + 1].abs());
          }
          REMAX = ONE / EMAX;
          dscal(N, REMAX, VL(1, KI).asArray(), 1);
          dscal(N, REMAX, VL(1, KI + 1).asArray(), 1);
        } else {
          // version 2: back-transform block of vectors with GEMM
          // zero out above vector
          // could go from KI-NV+1 to KI-1
          for (K = 1; K <= KI - 1; K++) {
            WORK[K + IV * N] = ZERO;
            WORK[K + (IV + 1) * N] = ZERO;
          }
          ISCOMPLEX[IV] = IP;
          ISCOMPLEX[IV + 1] = -IP;
          IV++;
          // back-transform and normalization is done below
        }
      }

      if (NB > 1) {
        // Blocked version of back-transform
        // For complex case, KI2 includes both vectors (KI and KI+1)
        if (IP == 0) {
          KI2 = KI;
        } else {
          KI2 = KI + 1;
        }

        // Columns 1:IV of work are valid vectors.
        // When the number of vectors stored reaches NB-1 or NB,
        // or if this was last vector, do the GEMM
        if ((IV >= NB - 1) || (KI2 == N)) {
          dgemm(
              'N',
              'N',
              N,
              IV,
              N - KI2 + IV,
              ONE,
              VL(1, KI2 - IV + 1),
              LDVL,
              WORK(KI2 - IV + 1 + 1 * N).asMatrix(N),
              N,
              ZERO,
              WORK(1 + (NB + 1) * N).asMatrix(N),
              N);
          // normalize vectors
          for (K = 1; K <= IV; K++) {
            if (ISCOMPLEX[K] == 0) {
              // real eigenvector
              II = idamax(N, WORK(1 + (NB + K) * N), 1);
              REMAX = ONE / WORK[II + (NB + K) * N].abs();
            } else if (ISCOMPLEX[K] == 1) {
              // first eigenvector of conjugate pair
              EMAX = ZERO;
              for (II = 1; II <= N; II++) {
                EMAX = max(
                    EMAX,
                    WORK[II + (NB + K) * N].abs() +
                        WORK[II + (NB + K + 1) * N].abs());
              }
              REMAX = ONE / EMAX;
              // else if ISCOMPLEX[K] == -1
              // second eigenvector of conjugate pair
              // reuse same REMAX as previous K
            }
            dscal(N, REMAX, WORK(1 + (NB + K) * N), 1);
          }
          dlacpy('F', N, IV, WORK(1 + (NB + 1) * N).asMatrix(N), N,
              VL(1, KI2 - IV + 1), LDVL);
          IV = 1;
        } else {
          IV++;
        }
      } // ! blocked back-transform;

      IS++;
      if (IP != 0) IS++;
    }
  }
}
