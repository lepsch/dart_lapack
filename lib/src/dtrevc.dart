import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaln2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtrevc(
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
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SELECT = SELECT_.having();
  final T = T_.having(ld: LDT);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool ALLV, BOTHV, LEFTV, OVER, PAIR, RIGHTV, SOMEV;
  int I, II, IP, IS, J, J1, J2, JNXT, K, KI, N2;
  double BETA,
      BIGNUM,
      EMAX,
      // OVFL,
      REC,
      REMAX,
      SMIN,
      SMLNUM,
      ULP,
      UNFL,
      VCRIT,
      VMAX,
      WI,
      WR;
  final X = Matrix<double>(2, 2);
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
  } else {
    // Set M.value to the number of columns required to store the selected
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
              if (SELECT[J]) M.value = M.value + 1;
            } else {
              PAIR = true;
              if (SELECT[J] || SELECT[J + 1]) {
                SELECT[J] = true;
                M.value = M.value + 2;
              }
            }
          } else {
            if (SELECT[N]) M.value = M.value + 1;
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
    xerbla('DTREVC', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  // Set the constants to control overflow.

  UNFL = dlamch('Safe minimum');
  // OVFL = ONE / UNFL;
  ULP = dlamch('Precision');
  SMLNUM = UNFL * (N / ULP);
  BIGNUM = (ONE - ULP) / SMLNUM;

  // Compute 1-norm of each column of strictly upper triangular
  // part of T to control overflow in triangular solver.

  WORK[1] = ZERO;
  for (J = 2; J <= N; J++) {
    WORK[J] = ZERO;
    for (I = 1; I <= J - 1; I++) {
      WORK[J] = WORK[J] + (T[I][J]).abs();
    }
  }

  // Index IP is used to specify the real or complex eigenvalue:
  // IP = 0, real eigenvalue,
  // 1, first of conjugate complex pair: (wr,wi)
  // -1, second of conjugate complex pair: (wr,wi)

  N2 = 2 * N;

  if (RIGHTV) {
    // Compute right eigenvectors.

    IP = 0;
    IS = M.value;
    for (KI = N; KI >= 1; KI--) {
      while (true) {
        if (IP == 1) break;
        if (KI != 1 && T[KI][KI - 1] != ZERO) {
          IP = -1;
        }

        if (SOMEV) {
          if (IP == 0) {
            if (!SELECT[KI]) break;
          } else {
            if (!SELECT[KI - 1]) break;
          }
        }

        // Compute the KI-th eigenvalue (WR,WI).

        WR = T[KI][KI];
        WI = ZERO;
        if (IP != 0) {
          WI = sqrt((T[KI][KI - 1]).abs()) * sqrt((T[KI - 1][KI]).abs());
        }
        SMIN = max(ULP * ((WR).abs() + (WI).abs()), SMLNUM);

        if (IP == 0) {
          // Real right eigenvector

          WORK[KI + N] = ONE;

          // Form right-hand side

          for (K = 1; K <= KI - 1; K++) {
            WORK[K + N] = -T[K][KI];
          }

          // Solve the upper quasi-triangular system:
          // (T[1:KI-1][1:KI-1] - WR)*X = SCALE.value*WORK.

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
                  WORK(J + N).asMatrix(N),
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
                  X[1][1] = X[1][1] / XNORM.value;
                  SCALE.value = SCALE.value / XNORM.value;
                }
              }

              // Scale if necessary

              if (SCALE.value != ONE) dscal(KI, SCALE.value, WORK(1 + N), 1);
              WORK[J + N] = X[1][1];

              // Update right-hand side

              daxpy(J - 1, -X[1][1], T(1, J).asArray(), 1, WORK(1 + N), 1);
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
                  WORK(J - 1 + N).asMatrix(N),
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
                  X[1][1] = X[1][1] / XNORM.value;
                  X[2][1] = X[2][1] / XNORM.value;
                  SCALE.value = SCALE.value / XNORM.value;
                }
              }

              // Scale if necessary

              if (SCALE.value != ONE) dscal(KI, SCALE.value, WORK(1 + N), 1);
              WORK[J - 1 + N] = X[1][1];
              WORK[J + N] = X[2][1];

              // Update right-hand side

              daxpy(J - 2, -X[1][1], T(1, J - 1).asArray(), 1, WORK(1 + N), 1);
              daxpy(J - 2, -X[2][1], T(1, J).asArray(), 1, WORK(1 + N), 1);
            }
          }

          // Copy the vector x or Q*x to VR and normalize.

          if (!OVER) {
            dcopy(KI, WORK(1 + N), 1, VR(1, IS).asArray(), 1);

            II = idamax(KI, VR(1, IS).asArray(), 1);
            REMAX = ONE / (VR[II][IS]).abs();
            dscal(KI, REMAX, VR(1, IS).asArray(), 1);

            for (K = KI + 1; K <= N; K++) {
              VR[K][IS] = ZERO;
            }
          } else {
            if (KI > 1) {
              dgemv('N', N, KI - 1, ONE, VR, LDVR, WORK(1 + N), 1, WORK[KI + N],
                  VR(1, KI).asArray(), 1);
            }

            II = idamax(N, VR(1, KI).asArray(), 1);
            REMAX = ONE / (VR[II][KI]).abs();
            dscal(N, REMAX, VR(1, KI).asArray(), 1);
          }
        } else {
          // Complex right eigenvector.

          // Initial solve
          // [ (T[KI-1][KI-1] T[KI-1][KI] ) - (WR + I* WI)]*X = 0.
          // [ (T[KI][KI-1]   T[KI][KI]   )               ]

          if ((T[KI - 1][KI]).abs() >= (T[KI][KI - 1]).abs()) {
            WORK[KI - 1 + N] = ONE;
            WORK[KI + N2] = WI / T[KI - 1][KI];
          } else {
            WORK[KI - 1 + N] = -WI / T[KI][KI - 1];
            WORK[KI + N2] = ONE;
          }
          WORK[KI + N] = ZERO;
          WORK[KI - 1 + N2] = ZERO;

          // Form right-hand side

          for (K = 1; K <= KI - 2; K++) {
            WORK[K + N] = -WORK[KI - 1 + N] * T[K][KI - 1];
            WORK[K + N2] = -WORK[KI + N2] * T[K][KI];
          }

          // Solve upper quasi-triangular system:
          // (T[1:KI-2][1:KI-2] - (WR+i*WI))*X = SCALE.value*(WORK+i*WORK2)

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

              dlaln2(false, 1, 2, SMIN, ONE, T(J, J), LDT, ONE, ONE,
                  WORK(J + N).asMatrix(N), N, WR, WI, X, 2, SCALE, XNORM, IERR);

              // Scale X[1][1] and X[1][2] to avoid overflow when
              // updating the right-hand side.

              if (XNORM.value > ONE) {
                if (WORK[J] > BIGNUM / XNORM.value) {
                  X[1][1] = X[1][1] / XNORM.value;
                  X[1][2] = X[1][2] / XNORM.value;
                  SCALE.value = SCALE.value / XNORM.value;
                }
              }

              // Scale if necessary

              if (SCALE.value != ONE) {
                dscal(KI, SCALE.value, WORK(1 + N), 1);
                dscal(KI, SCALE.value, WORK(1 + N2), 1);
              }
              WORK[J + N] = X[1][1];
              WORK[J + N2] = X[1][2];

              // Update the right-hand side

              daxpy(J - 1, -X[1][1], T(1, J).asArray(), 1, WORK(1 + N), 1);
              daxpy(J - 1, -X[1][2], T(1, J).asArray(), 1, WORK(1 + N2), 1);
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
                  WORK(J - 1 + N).asMatrix(N),
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
                  X[1][1] = X[1][1] * REC;
                  X[1][2] = X[1][2] * REC;
                  X[2][1] = X[2][1] * REC;
                  X[2][2] = X[2][2] * REC;
                  SCALE.value = SCALE.value * REC;
                }
              }

              // Scale if necessary

              if (SCALE.value != ONE) {
                dscal(KI, SCALE.value, WORK(1 + N), 1);
                dscal(KI, SCALE.value, WORK(1 + N2), 1);
              }
              WORK[J - 1 + N] = X[1][1];
              WORK[J + N] = X[2][1];
              WORK[J - 1 + N2] = X[1][2];
              WORK[J + N2] = X[2][2];

              // Update the right-hand side

              daxpy(J - 2, -X[1][1], T(1, J - 1).asArray(), 1, WORK(1 + N), 1);
              daxpy(J - 2, -X[2][1], T(1, J).asArray(), 1, WORK(1 + N), 1);
              daxpy(J - 2, -X[1][2], T(1, J - 1).asArray(), 1, WORK(1 + N2), 1);
              daxpy(J - 2, -X[2][2], T(1, J).asArray(), 1, WORK(1 + N2), 1);
            }
          }

          // Copy the vector x or Q*x to VR and normalize.

          if (!OVER) {
            dcopy(KI, WORK(1 + N), 1, VR(1, IS - 1).asArray(), 1);
            dcopy(KI, WORK(1 + N2), 1, VR(1, IS).asArray(), 1);

            EMAX = ZERO;
            for (K = 1; K <= KI; K++) {
              EMAX = max(EMAX, (VR[K][IS - 1]).abs() + (VR[K][IS]).abs());
            }

            REMAX = ONE / EMAX;
            dscal(KI, REMAX, VR(1, IS - 1).asArray(), 1);
            dscal(KI, REMAX, VR(1, IS).asArray(), 1);

            for (K = KI + 1; K <= N; K++) {
              VR[K][IS - 1] = ZERO;
              VR[K][IS] = ZERO;
            }
          } else {
            if (KI > 2) {
              dgemv('N', N, KI - 2, ONE, VR, LDVR, WORK(1 + N), 1,
                  WORK[KI - 1 + N], VR(1, KI - 1).asArray(), 1);
              dgemv('N', N, KI - 2, ONE, VR, LDVR, WORK(1 + N2), 1,
                  WORK[KI + N2], VR(1, KI).asArray(), 1);
            } else {
              dscal(N, WORK[KI - 1 + N], VR(1, KI - 1).asArray(), 1);
              dscal(N, WORK[KI + N2], VR(1, KI).asArray(), 1);
            }

            EMAX = ZERO;
            for (K = 1; K <= N; K++) {
              EMAX = max(EMAX, (VR[K][KI - 1]).abs() + (VR[K][KI]).abs());
            }
            REMAX = ONE / EMAX;
            dscal(N, REMAX, VR(1, KI - 1).asArray(), 1);
            dscal(N, REMAX, VR(1, KI).asArray(), 1);
          }
        }

        IS = IS - 1;
        if (IP != 0) IS = IS - 1;
        break;
      }
      if (IP == 1) IP = 0;
      if (IP == -1) IP = 1;
    }
  }

  if (!LEFTV) return;

  // Compute left eigenvectors.

  IP = 0;
  IS = 1;
  for (KI = 1; KI <= N; KI++) {
    while (true) {
      if (IP == -1) break;
      if (KI != N && T[KI + 1][KI] != ZERO) {
        IP = 1;
      }
      if (SOMEV) {
        if (!SELECT[KI]) break;
      }

      // Compute the KI-th eigenvalue (WR,WI).

      WR = T[KI][KI];
      WI = ZERO;
      if (IP != 0) {
        WI = sqrt((T[KI][KI + 1]).abs()) * sqrt((T[KI + 1][KI]).abs());
      }
      SMIN = max(ULP * ((WR).abs() + (WI).abs()), SMLNUM);

      if (IP == 0) {
        // Real left eigenvector.

        WORK[KI + N] = ONE;

        // Form right-hand side

        for (K = KI + 1; K <= N; K++) {
          WORK[K + N] = -T[KI][K];
        }

        // Solve the quasi-triangular system:
        // (T[KI+1:N][KI+1:N] - WR)**T*X = SCALE.value*WORK

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
              dscal(N - KI + 1, REC, WORK(KI + N), 1);
              VMAX = ONE;
              VCRIT = BIGNUM;
            }

            WORK[J + N] = WORK[J + N] -
                ddot(
                    J - KI - 1, T(KI + 1, J).asArray(), 1, WORK(KI + 1 + N), 1);

            // Solve (T[J][J]-WR)**T*X = WORK

            dlaln2(false, 1, 1, SMIN, ONE, T(J, J), LDT, ONE, ONE,
                WORK(J + N).asMatrix(N), N, WR, ZERO, X, 2, SCALE, XNORM, IERR);

            // Scale if necessary

            if (SCALE.value != ONE) {
              dscal(N - KI + 1, SCALE.value, WORK(KI + N), 1);
            }
            WORK[J + N] = X[1][1];
            VMAX = max((WORK[J + N]).abs(), VMAX);
            VCRIT = BIGNUM / VMAX;
          } else {
            // 2-by-2 diagonal block

            // Scale if necessary to avoid overflow when forming
            // the right-hand side.

            BETA = max(WORK[J], WORK[J + 1]);
            if (BETA > VCRIT) {
              REC = ONE / VMAX;
              dscal(N - KI + 1, REC, WORK(KI + N), 1);
              VMAX = ONE;
              VCRIT = BIGNUM;
            }

            WORK[J + N] = WORK[J + N] -
                ddot(
                    J - KI - 1, T(KI + 1, J).asArray(), 1, WORK(KI + 1 + N), 1);

            WORK[J + 1 + N] = WORK[J + 1 + N] -
                ddot(J - KI - 1, T(KI + 1, J + 1).asArray(), 1,
                    WORK(KI + 1 + N), 1);

            // Solve
            // [T[J][J]-WR   T[J][J+1]     ]**T * X = SCALE.value*( WORK1 )
            // [T[J+1][J]    T[J+1][J+1]-WR]                ( WORK2 )

            dlaln2(true, 2, 1, SMIN, ONE, T(J, J), LDT, ONE, ONE,
                WORK(J + N).asMatrix(N), N, WR, ZERO, X, 2, SCALE, XNORM, IERR);

            // Scale if necessary

            if (SCALE.value != ONE) {
              dscal(N - KI + 1, SCALE.value, WORK(KI + N), 1);
            }
            WORK[J + N] = X[1][1];
            WORK[J + 1 + N] = X[2][1];

            VMAX = max((WORK[J + N]).abs(), max((WORK[J + 1 + N]).abs(), VMAX));
            VCRIT = BIGNUM / VMAX;
          }
        }

        // Copy the vector x or Q*x to VL and normalize.

        if (!OVER) {
          dcopy(N - KI + 1, WORK(KI + N), 1, VL(KI, IS).asArray(), 1);

          II = idamax(N - KI + 1, VL(KI, IS).asArray(), 1) + KI - 1;
          REMAX = ONE / (VL[II][IS]).abs();
          dscal(N - KI + 1, REMAX, VL(KI, IS).asArray(), 1);

          for (K = 1; K <= KI - 1; K++) {
            VL[K][IS] = ZERO;
          }
        } else {
          if (KI < N) {
            dgemv('N', N, N - KI, ONE, VL(1, KI + 1), LDVL, WORK(KI + 1 + N), 1,
                WORK[KI + N], VL(1, KI).asArray(), 1);
          }

          II = idamax(N, VL(1, KI).asArray(), 1);
          REMAX = ONE / (VL[II][KI]).abs();
          dscal(N, REMAX, VL(1, KI).asArray(), 1);
        }
      } else {
        // Complex left eigenvector.

        // Initial solve:
        // ((T[KI][KI]    T[KI][KI+1] )**T - (WR - I* WI))*X = 0.
        // ((T[KI+1][KI] T[KI+1][KI+1])                )

        if ((T[KI][KI + 1]).abs() >= (T[KI + 1][KI]).abs()) {
          WORK[KI + N] = WI / T[KI][KI + 1];
          WORK[KI + 1 + N2] = ONE;
        } else {
          WORK[KI + N] = ONE;
          WORK[KI + 1 + N2] = -WI / T[KI + 1][KI];
        }
        WORK[KI + 1 + N] = ZERO;
        WORK[KI + N2] = ZERO;

        // Form right-hand side

        for (K = KI + 2; K <= N; K++) {
          WORK[K + N] = -WORK[KI + N] * T[KI][K];
          WORK[K + N2] = -WORK[KI + 1 + N2] * T[KI + 1][K];
        }

        // Solve complex quasi-triangular system:
        // ( T[KI+2,N:KI+2][N] - (WR-i*WI) )*X = WORK1+i*WORK2

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
              dscal(N - KI + 1, REC, WORK(KI + N), 1);
              dscal(N - KI + 1, REC, WORK(KI + N2), 1);
              VMAX = ONE;
              VCRIT = BIGNUM;
            }

            WORK[J + N] = WORK[J + N] -
                ddot(
                    J - KI - 2, T(KI + 2, J).asArray(), 1, WORK(KI + 2 + N), 1);
            WORK[J + N2] = WORK[J + N2] -
                ddot(J - KI - 2, T(KI + 2, J).asArray(), 1, WORK(KI + 2 + N2),
                    1);

            // Solve (T[J][J]-(WR-i*WI))*(X11+i*X12)= WK+I*WK2

            dlaln2(false, 1, 2, SMIN, ONE, T(J, J), LDT, ONE, ONE,
                WORK(J + N).asMatrix(N), N, WR, -WI, X, 2, SCALE, XNORM, IERR);

            // Scale if necessary

            if (SCALE.value != ONE) {
              dscal(N - KI + 1, SCALE.value, WORK(KI + N), 1);
              dscal(N - KI + 1, SCALE.value, WORK(KI + N2), 1);
            }
            WORK[J + N] = X[1][1];
            WORK[J + N2] = X[1][2];
            VMAX = max((WORK[J + N]).abs(), max((WORK[J + N2]).abs(), VMAX));
            VCRIT = BIGNUM / VMAX;
          } else {
            // 2-by-2 diagonal block

            // Scale if necessary to avoid overflow when forming
            // the right-hand side elements.

            BETA = max(WORK[J], WORK[J + 1]);
            if (BETA > VCRIT) {
              REC = ONE / VMAX;
              dscal(N - KI + 1, REC, WORK(KI + N), 1);
              dscal(N - KI + 1, REC, WORK(KI + N2), 1);
              VMAX = ONE;
              VCRIT = BIGNUM;
            }

            WORK[J + N] = WORK[J + N] -
                ddot(
                    J - KI - 2, T(KI + 2, J).asArray(), 1, WORK(KI + 2 + N), 1);

            WORK[J + N2] = WORK[J + N2] -
                ddot(J - KI - 2, T(KI + 2, J).asArray(), 1, WORK(KI + 2 + N2),
                    1);

            WORK[J + 1 + N] = WORK[J + 1 + N] -
                ddot(J - KI - 2, T(KI + 2, J + 1).asArray(), 1,
                    WORK(KI + 2 + N), 1);

            WORK[J + 1 + N2] = WORK[J + 1 + N2] -
                ddot(J - KI - 2, T(KI + 2, J + 1).asArray(), 1,
                    WORK(KI + 2 + N2), 1);

            // Solve 2-by-2 complex linear equation
            // ([T[j][j]   T[j][j+1]  ]**T-(wr-i*wi)*I)*X = SCALE.value*B
            // ([T[j+1][j] T[j+1][j+1]]               )

            dlaln2(true, 2, 2, SMIN, ONE, T(J, J), LDT, ONE, ONE,
                WORK(J + N).asMatrix(N), N, WR, -WI, X, 2, SCALE, XNORM, IERR);

            // Scale if necessary

            if (SCALE.value != ONE) {
              dscal(N - KI + 1, SCALE.value, WORK(KI + N), 1);
              dscal(N - KI + 1, SCALE.value, WORK(KI + N2), 1);
            }
            WORK[J + N] = X[1][1];
            WORK[J + N2] = X[1][2];
            WORK[J + 1 + N] = X[2][1];
            WORK[J + 1 + N2] = X[2][2];
            VMAX = max(
                (X[1][1]).abs(),
                max(
                  max((X[1][2]).abs(), (X[2][1]).abs()),
                  max((X[2][2]).abs(), VMAX),
                ));
            VCRIT = BIGNUM / VMAX;
          }
        }

        // Copy the vector x or Q*x to VL and normalize.

        if (!OVER) {
          dcopy(N - KI + 1, WORK(KI + N), 1, VL(KI, IS).asArray(), 1);
          dcopy(N - KI + 1, WORK(KI + N2), 1, VL(KI, IS + 1).asArray(), 1);

          EMAX = ZERO;
          for (K = KI; K <= N; K++) {
            EMAX = max(EMAX, (VL[K][IS]).abs() + (VL[K][IS + 1]).abs());
          }
          REMAX = ONE / EMAX;
          dscal(N - KI + 1, REMAX, VL(KI, IS).asArray(), 1);
          dscal(N - KI + 1, REMAX, VL(KI, IS + 1).asArray(), 1);

          for (K = 1; K <= KI - 1; K++) {
            VL[K][IS] = ZERO;
            VL[K][IS + 1] = ZERO;
          }
        } else {
          if (KI < N - 1) {
            dgemv('N', N, N - KI - 1, ONE, VL(1, KI + 2), LDVL,
                WORK(KI + 2 + N), 1, WORK[KI + N], VL(1, KI).asArray(), 1);
            dgemv(
                'N',
                N,
                N - KI - 1,
                ONE,
                VL(1, KI + 2),
                LDVL,
                WORK(KI + 2 + N2),
                1,
                WORK[KI + 1 + N2],
                VL(1, KI + 1).asArray(),
                1);
          } else {
            dscal(N, WORK[KI + N], VL(1, KI).asArray(), 1);
            dscal(N, WORK[KI + 1 + N2], VL(1, KI + 1).asArray(), 1);
          }

          EMAX = ZERO;
          for (K = 1; K <= N; K++) {
            EMAX = max(EMAX, (VL[K][KI]).abs() + (VL[K][KI + 1]).abs());
          }
          REMAX = ONE / EMAX;
          dscal(N, REMAX, VL(1, KI).asArray(), 1);
          dscal(N, REMAX, VL(1, KI + 1).asArray(), 1);
        }
      }

      IS = IS + 1;
      if (IP != 0) IS = IS + 1;
      break;
    }
    if (IP == -1) IP = 0;
    if (IP == 1) IP = -1;
  }
}
