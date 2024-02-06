import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dladiv.dart';
import 'package:lapack/src/dlaln2.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlaqtr(
  final bool LTRAN,
  final bool LREAL,
  final int N,
  final Matrix<double> T,
  final int LDT,
  final Array<double> B,
  final double W,
  final Box<double> SCALE,
  final Array<double> X,
  final Array<double> WORK,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool NOTRAN;
  int I, J, J1, J2, JNEXT, K, N1, N2;
  double BIGNUM, EPS, REC, SMIN, SMINW, SMLNUM, TJJ, TMP, XJ, XMAX, Z;
  final D = Matrix<double>(2, 2), V = Matrix<double>(2, 2);
  final IERR = Box(0);
  final SCALOC = Box(0.0), XNORM = Box(0.0), SR = Box(0.0), SI = Box(0.0);

  // Do not test the input parameters for errors

  NOTRAN = !LTRAN;
  INFO.value = 0;

  // Quick return if possible

  if (N == 0) return;

  // Set constants to control overflow

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  XNORM.value = dlange('M', N, N, T, LDT, D.asArray());
  if (!LREAL) {
    XNORM.value = max(
      XNORM.value,
      max(W.abs(), dlange('M', N, 1, B.asMatrix(N), N, D.asArray())),
    );
  }
  SMIN = max(SMLNUM, EPS * XNORM.value);

  // Compute 1-norm of each column of strictly upper triangular
  // part of T to control overflow in triangular solver.

  WORK[1] = ZERO;
  for (J = 2; J <= N; J++) {
    WORK[J] = dasum(J - 1, T(1, J).asArray(), 1);
  }

  if (!LREAL) {
    for (I = 2; I <= N; I++) {
      WORK[I] = WORK[I] + B[I].abs();
    }
  }

  N2 = 2 * N;
  N1 = N;
  if (!LREAL) N1 = N2;
  K = idamax(N1, X, 1);
  XMAX = X[K].abs();
  SCALE.value = ONE;

  if (XMAX > BIGNUM) {
    SCALE.value = BIGNUM / XMAX;
    dscal(N1, SCALE.value, X, 1);
    XMAX = BIGNUM;
  }

  if (LREAL) {
    if (NOTRAN) {
      // Solve T*p = scale*c

      JNEXT = N;
      for (J = N; J >= 1; J--) {
        if (J > JNEXT) continue;
        J1 = J;
        J2 = J;
        JNEXT = J - 1;
        if (J > 1) {
          if (T[J][J - 1] != ZERO) {
            J1 = J - 1;
            JNEXT = J - 2;
          }
        }

        if (J1 == J2) {
          // Meet 1 by 1 diagonal block

          // Scale to avoid overflow when computing
          // x[j] = b[j]/T[j][j]

          XJ = X[J1].abs();
          TJJ = T[J1][J1].abs();
          TMP = T[J1][J1];
          if (TJJ < SMIN) {
            TMP = SMIN;
            TJJ = SMIN;
            INFO.value = 1;
          }

          if (XJ == ZERO) continue;

          if (TJJ < ONE) {
            if (XJ > BIGNUM * TJJ) {
              REC = ONE / XJ;
              dscal(N, REC, X, 1);
              SCALE.value = SCALE.value * REC;
              XMAX = XMAX * REC;
            }
          }
          X[J1] = X[J1] / TMP;
          XJ = X[J1].abs();

          // Scale x if necessary to avoid overflow when adding a
          // multiple of column j1 of T.

          if (XJ > ONE) {
            REC = ONE / XJ;
            if (WORK[J1] > (BIGNUM - XMAX) * REC) {
              dscal(N, REC, X, 1);
              SCALE.value = SCALE.value * REC;
            }
          }
          if (J1 > 1) {
            daxpy(J1 - 1, -X[J1], T(1, J1).asArray(), 1, X, 1);
            K = idamax(J1 - 1, X, 1);
            XMAX = X[K].abs();
          }
        } else {
          // Meet 2 by 2 diagonal block

          // Call 2 by 2 linear system solve, to take
          // care of possible overflow by scaling factor.

          D[1][1] = X[J1];
          D[2][1] = X[J2];
          dlaln2(
            false,
            2,
            1,
            SMIN,
            ONE,
            T(J1, J1),
            LDT,
            ONE,
            ONE,
            D,
            2,
            ZERO,
            ZERO,
            V,
            2,
            SCALOC,
            XNORM,
            IERR,
          );
          if (IERR.value != 0) INFO.value = 2;

          if (SCALOC.value != ONE) {
            dscal(N, SCALOC.value, X, 1);
            SCALE.value = SCALE.value * SCALOC.value;
          }
          X[J1] = V[1][1];
          X[J2] = V[2][1];

          // Scale V[1][1] (= X[J1]) and/or V[2][1] (=X[J2])
          // to avoid overflow in updating right-hand side.

          XJ = max((V[1][1]).abs(), (V[2][1]).abs());
          if (XJ > ONE) {
            REC = ONE / XJ;
            if (max(WORK[J1], WORK[J2]) > (BIGNUM - XMAX) * REC) {
              dscal(N, REC, X, 1);
              SCALE.value = SCALE.value * REC;
            }
          }

          // Update right-hand side

          if (J1 > 1) {
            daxpy(J1 - 1, -X[J1], T(1, J1).asArray(), 1, X, 1);
            daxpy(J1 - 1, -X[J2], T(1, J2).asArray(), 1, X, 1);
            K = idamax(J1 - 1, X, 1);
            XMAX = (X[K]).abs();
          }
        }
      }
    } else {
      // Solve T**T*p = scale*c

      JNEXT = 1;
      for (J = 1; J <= N; J++) {
        if (J < JNEXT) continue;
        J1 = J;
        J2 = J;
        JNEXT = J + 1;
        if (J < N) {
          if (T[J + 1][J] != ZERO) {
            J2 = J + 1;
            JNEXT = J + 2;
          }
        }

        if (J1 == J2) {
          // 1 by 1 diagonal block

          // Scale if necessary to avoid overflow in forming the
          // right-hand side element by inner product.

          XJ = (X[J1]).abs();
          if (XMAX > ONE) {
            REC = ONE / XMAX;
            if (WORK[J1] > (BIGNUM - XJ) * REC) {
              dscal(N, REC, X, 1);
              SCALE.value = SCALE.value * REC;
              XMAX = XMAX * REC;
            }
          }

          X[J1] = X[J1] - ddot(J1 - 1, T(1, J1).asArray(), 1, X, 1);

          XJ = X[J1].abs();
          TJJ = T[J1][J1].abs();
          TMP = T[J1][J1];
          if (TJJ < SMIN) {
            TMP = SMIN;
            TJJ = SMIN;
            INFO.value = 1;
          }

          if (TJJ < ONE) {
            if (XJ > BIGNUM * TJJ) {
              REC = ONE / XJ;
              dscal(N, REC, X, 1);
              SCALE.value = SCALE.value * REC;
              XMAX = XMAX * REC;
            }
          }
          X[J1] = X[J1] / TMP;
          XMAX = max(XMAX, (X[J1]).abs());
        } else {
          // 2 by 2 diagonal block

          // Scale if necessary to avoid overflow in forming the
          // right-hand side elements by inner product.

          XJ = max((X[J1]).abs(), (X[J2]).abs());
          if (XMAX > ONE) {
            REC = ONE / XMAX;
            if (max(WORK[J2], WORK[J1]) > (BIGNUM - XJ) * REC) {
              dscal(N, REC, X, 1);
              SCALE.value = SCALE.value * REC;
              XMAX = XMAX * REC;
            }
          }

          D[1][1] = X[J1] - ddot(J1 - 1, T(1, J1).asArray(), 1, X, 1);
          D[2][1] = X[J2] - ddot(J1 - 1, T(1, J2).asArray(), 1, X, 1);

          dlaln2(
            true,
            2,
            1,
            SMIN,
            ONE,
            T(J1, J1),
            LDT,
            ONE,
            ONE,
            D,
            2,
            ZERO,
            ZERO,
            V,
            2,
            SCALOC,
            XNORM,
            IERR,
          );
          if (IERR.value != 0) INFO.value = 2;

          if (SCALOC.value != ONE) {
            dscal(N, SCALOC.value, X, 1);
            SCALE.value = SCALE.value * SCALOC.value;
          }
          X[J1] = V[1][1];
          X[J2] = V[2][1];
          XMAX = max(X[J1].abs(), max(X[J2].abs(), XMAX));
        }
      }
    }
  } else {
    SMINW = max(EPS * (W).abs(), SMIN);
    if (NOTRAN) {
      // Solve (T + iB)*(p+iq) = c+id

      JNEXT = N;
      for (J = N; J >= 1; J--) {
        if (J > JNEXT) continue;
        J1 = J;
        J2 = J;
        JNEXT = J - 1;
        if (J > 1) {
          if (T[J][J - 1] != ZERO) {
            J1 = J - 1;
            JNEXT = J - 2;
          }
        }

        if (J1 == J2) {
          // 1 by 1 diagonal block

          // Scale if necessary to avoid overflow in division

          Z = W;
          if (J1 == 1) Z = B[1];
          XJ = X[J1].abs() + X[N + J1].abs();
          TJJ = T[J1][J1].abs() + Z.abs();
          TMP = T[J1][J1];
          if (TJJ < SMINW) {
            TMP = SMINW;
            TJJ = SMINW;
            INFO.value = 1;
          }

          if (XJ == ZERO) continue;

          if (TJJ < ONE) {
            if (XJ > BIGNUM * TJJ) {
              REC = ONE / XJ;
              dscal(N2, REC, X, 1);
              SCALE.value = SCALE.value * REC;
              XMAX = XMAX * REC;
            }
          }
          dladiv(X[J1], X[N + J1], TMP, Z, SR.value, SI.value);
          X[J1] = SR.value;
          X[N + J1] = SI.value;
          XJ = X[J1].abs() + X[N + J1].abs();

          // Scale x if necessary to avoid overflow when adding a
          // multiple of column j1 of T.

          if (XJ > ONE) {
            REC = ONE / XJ;
            if (WORK[J1] > (BIGNUM - XMAX) * REC) {
              dscal(N2, REC, X, 1);
              SCALE.value = SCALE.value * REC;
            }
          }

          if (J1 > 1) {
            daxpy(J1 - 1, -X[J1], T(1, J1).asArray(), 1, X, 1);
            daxpy(J1 - 1, -X[N + J1], T(1, J1).asArray(), 1, X(N + 1), 1);

            X[1] = X[1] + B[J1] * X[N + J1];
            X[N + 1] = X[N + 1] - B[J1] * X[J1];

            XMAX = ZERO;
            for (K = 1; K <= J1 - 1; K++) {
              XMAX = max(XMAX, (X[K]).abs() + (X[K + N]).abs());
            }
          }
        } else {
          // Meet 2 by 2 diagonal block

          D[1][1] = X[J1];
          D[2][1] = X[J2];
          D[1][2] = X[N + J1];
          D[2][2] = X[N + J2];
          dlaln2(
            false,
            2,
            2,
            SMINW,
            ONE,
            T(J1, J1),
            LDT,
            ONE,
            ONE,
            D,
            2,
            ZERO,
            -W,
            V,
            2,
            SCALOC,
            XNORM,
            IERR,
          );
          if (IERR.value != 0) INFO.value = 2;

          if (SCALOC.value != ONE) {
            dscal(2 * N, SCALOC.value, X, 1);
            SCALE.value = SCALOC.value * SCALE.value;
          }
          X[J1] = V[1][1];
          X[J2] = V[2][1];
          X[N + J1] = V[1][2];
          X[N + J2] = V[2][2];

          // Scale X[J1], .... to avoid overflow in
          // updating right hand side.

          XJ =
              max(V[1][1].abs() + V[1][2].abs(), V[2][1].abs() + V[2][2].abs());
          if (XJ > ONE) {
            REC = ONE / XJ;
            if (max(WORK[J1], WORK[J2]) > (BIGNUM - XMAX) * REC) {
              dscal(N2, REC, X, 1);
              SCALE.value = SCALE.value * REC;
            }
          }

          // Update the right-hand side.

          if (J1 > 1) {
            daxpy(J1 - 1, -X[J1], T(1, J1).asArray(), 1, X, 1);
            daxpy(J1 - 1, -X[J2], T(1, J2).asArray(), 1, X, 1);

            daxpy(J1 - 1, -X[N + J1], T(1, J1).asArray(), 1, X(N + 1), 1);
            daxpy(J1 - 1, -X[N + J2], T(1, J2).asArray(), 1, X(N + 1), 1);

            X[1] = X[1] + B[J1] * X[N + J1] + B[J2] * X[N + J2];
            X[N + 1] = X[N + 1] - B[J1] * X[J1] - B[J2] * X[J2];

            XMAX = ZERO;
            for (K = 1; K <= J1 - 1; K++) {
              XMAX = max((X[K]).abs() + (X[K + N]).abs(), XMAX);
            }
          }
        }
      }
    } else {
      // Solve (T + iB)**T*(p+iq) = c+id

      JNEXT = 1;
      for (J = 1; J <= N; J++) {
        if (J < JNEXT) continue;
        J1 = J;
        J2 = J;
        JNEXT = J + 1;
        if (J < N) {
          if (T[J + 1][J] != ZERO) {
            J2 = J + 1;
            JNEXT = J + 2;
          }
        }

        if (J1 == J2) {
          // 1 by 1 diagonal block

          // Scale if necessary to avoid overflow in forming the
          // right-hand side element by inner product.

          XJ = (X[J1]).abs() + (X[J1 + N]).abs();
          if (XMAX > ONE) {
            REC = ONE / XMAX;
            if (WORK[J1] > (BIGNUM - XJ) * REC) {
              dscal(N2, REC, X, 1);
              SCALE.value = SCALE.value * REC;
              XMAX = XMAX * REC;
            }
          }

          X[J1] = X[J1] - ddot(J1 - 1, T(1, J1).asArray(), 1, X, 1);
          X[N + J1] =
              X[N + J1] - ddot(J1 - 1, T(1, J1).asArray(), 1, X(N + 1), 1);
          if (J1 > 1) {
            X[J1] = X[J1] - B[J1] * X[N + 1];
            X[N + J1] = X[N + J1] + B[J1] * X[1];
          }
          XJ = (X[J1]).abs() + (X[J1 + N]).abs();

          Z = W;
          if (J1 == 1) Z = B[1];

          // Scale if necessary to avoid overflow in
          // complex division

          TJJ = (T[J1][J1]).abs() + (Z).abs();
          TMP = T[J1][J1];
          if (TJJ < SMINW) {
            TMP = SMINW;
            TJJ = SMINW;
            INFO.value = 1;
          }

          if (TJJ < ONE) {
            if (XJ > BIGNUM * TJJ) {
              REC = ONE / XJ;
              dscal(N2, REC, X, 1);
              SCALE.value = SCALE.value * REC;
              XMAX = XMAX * REC;
            }
          }
          dladiv(X[J1], X[N + J1], TMP, -Z, SR.value, SI.value);
          X[J1] = SR.value;
          X[J1 + N] = SI.value;
          XMAX = max((X[J1]).abs() + (X[J1 + N]).abs(), XMAX);
        } else {
          // 2 by 2 diagonal block

          // Scale if necessary to avoid overflow in forming the
          // right-hand side element by inner product.

          XJ = max(
            (X[J1]).abs() + (X[N + J1]).abs(),
            (X[J2]).abs() + (X[N + J2]),
          ).abs();
          if (XMAX > ONE) {
            REC = ONE / XMAX;
            if (max(WORK[J1], WORK[J2]) > (BIGNUM - XJ) / XMAX) {
              dscal(N2, REC, X, 1);
              SCALE.value = SCALE.value * REC;
              XMAX = XMAX * REC;
            }
          }

          D[1][1] = X[J1] - ddot(J1 - 1, T(1, J1).asArray(), 1, X, 1);
          D[2][1] = X[J2] - ddot(J1 - 1, T(1, J2).asArray(), 1, X, 1);
          D[1][2] =
              X[N + J1] - ddot(J1 - 1, T(1, J1).asArray(), 1, X(N + 1), 1);
          D[2][2] =
              X[N + J2] - ddot(J1 - 1, T(1, J2).asArray(), 1, X(N + 1), 1);
          D[1][1] = D[1][1] - B[J1] * X[N + 1];
          D[2][1] = D[2][1] - B[J2] * X[N + 1];
          D[1][2] = D[1][2] + B[J1] * X[1];
          D[2][2] = D[2][2] + B[J2] * X[1];

          dlaln2(
            true,
            2,
            2,
            SMINW,
            ONE,
            T(J1, J1),
            LDT,
            ONE,
            ONE,
            D,
            2,
            ZERO,
            W,
            V,
            2,
            SCALOC,
            XNORM,
            IERR,
          );
          if (IERR.value != 0) INFO.value = 2;

          if (SCALOC.value != ONE) {
            dscal(N2, SCALOC.value, X, 1);
            SCALE.value = SCALOC.value * SCALE.value;
          }
          X[J1] = V[1][1];
          X[J2] = V[2][1];
          X[N + J1] = V[1][2];
          X[N + J2] = V[2][2];
          XMAX = max(
            (X[J1]).abs() + (X[N + J1]).abs(),
            max((X[J2]).abs() + (X[N + J2]).abs(), XMAX),
          );
        }
      }
    }
  }
}
