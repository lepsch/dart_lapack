// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/drot.dart';
import 'package:dart_lapack/src/blas/drotg.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarnv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/matrix.dart';

import '../matgen/dlarnd.dart';
import '../matgen/dlatms.dart';
import 'dlatb4.dart';

void dlattp(
  final int IMAT,
  final String UPLO,
  final String TRANS,
  final Box<String> DIAG,
  final Array<int> ISEED_,
  final int N,
  final Array<double> A_,
  final Array<double> B_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having();
  final B = B_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, TWO = 2.0, ZERO = 0.0;
  final PATH = '${'Double precision'[0]}TP';
  final UNFL = dlamch('Safe minimum');
  final ULP = dlamch('Epsilon') * dlamch('Base');
  final SMLNUM = UNFL;
  final BIGNUM = (ONE - ULP) / SMLNUM;
  if ((IMAT >= 7 && IMAT <= 10) || IMAT == 18) {
    DIAG.value = 'U';
  } else {
    DIAG.value = 'N';
  }
  INFO.value = 0;

  // Quick return if N <= 0.

  if (N <= 0) return;

  // Call DLATB4 to set parameters for DLATMS.

  final UPPER = lsame(UPLO, 'U');
  final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
      dlatb4(PATH, UPPER ? IMAT : -IMAT, N, N);
  final PACKIT = UPPER ? 'C' : 'R';

  // IMAT <= 6:  Non-unit triangular matrix

  if (IMAT <= 6) {
    dlatms(N, N, DIST, ISEED, TYPE, B, MODE, CNDNUM, ANORM, KL, KU, PACKIT,
        A.asMatrix(), N, WORK, INFO);

    // IMAT > 6:  Unit triangular matrix
    // The diagonal is deliberately set to something other than 1.

    // IMAT = 7:  Matrix is the identity
  } else if (IMAT == 7) {
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J - 1; I++) {
          A[JC + I - 1] = ZERO;
        }
        A[JC + J - 1] = J.toDouble();
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        A[JC] = J.toDouble();
        for (var I = J + 1; I <= N; I++) {
          A[JC + I - J] = ZERO;
        }
        JC += N - J + 1;
      }
    }

    // IMAT > 7:  Non-trivial unit triangular matrix

    // Generate a unit triangular matrix T with condition CNDNUM by
    // forming a triangular matrix with known singular values and
    // filling in the zero entries with Givens rotations.
  } else if (IMAT <= 10) {
    if (UPPER) {
      var JC = 0;
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J - 1; I++) {
          A[JC + I] = ZERO;
        }
        A[JC + J] = J.toDouble();
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        A[JC] = J.toDouble();
        for (var I = J + 1; I <= N; I++) {
          A[JC + I - J] = ZERO;
        }
        JC += N - J + 1;
      }
    }

    // Since the trace of a unit triangular matrix is 1, the product
    // of its singular values must be 1.  Let s = sqrt(CNDNUM),
    // x = sqrt(s) - 1/sqrt(s), y = sqrt(2/(n-2))*x, and z = x**2.
    // The following triangular matrix has singular values s, 1, 1,
    // ..., 1, 1/s:

    // 1  y  y  y  ...  y  y  z
    //    1  0  0  ...  0  0  y
    //       1  0  ...  0  0  y
    //          .  ...  .  .  .
    //              .   .  .  .
    //                  1  0  y
    //                     1  y
    //                        1

    // To fill in the zeros, we first multiply by a matrix with small
    // condition number of the form

    // 1  0  0  0  0  ...
    //    1  +  *  0  0  ...
    //       1  +  0  0  0
    //          1  +  *  0  0
    //             1  +  0  0
    //                ...
    //                   1  +  0
    //                      1  0
    //                         1

    // Each element marked with a '*' is formed by taking the product
    // of the adjacent elements marked with '+'.  The '*'s can be
    // chosen freely, and the '+'s are chosen so that the inverse of
    // T will have elements of the same magnitude as T.  If the *'s in
    // both T and inv(T) have small magnitude, T is well conditioned.
    // The two offdiagonals of T are stored in WORK.

    // The product of these two matrices has the form

    // 1  y  y  y  y  y  .  y  y  z
    //    1  +  *  0  0  .  0  0  y
    //       1  +  0  0  .  0  0  y
    //          1  +  *  .  .  .  .
    //             1  +  .  .  .  .
    //                .  .  .  .  .
    //                   .  .  .  .
    //                      1  +  y
    //                         1  y
    //                            1

    // Now we multiply by Givens rotations, using the fact that

    // [  c   s ] [  1   w ] [ -c  -s ] =  [  1  -w ]
    // [ -s   c ] [  0   1 ] [  s  -c ]    [  0   1 ]
    // and
    //       [ -c  -s ] [  1   0 ] [  c   s ] =  [  1   0 ]
    //       [  s  -c ] [  w   1 ] [ -s   c ]    [ -w   1 ]

    // where c = w / sqrt(w**2+4) and s = 2 / sqrt(w**2+4).

    var STAR1 = 0.25;
    var SFAC = 0.5;
    var PLUS1 = SFAC;
    for (var J = 1; J <= N; J += 2) {
      final PLUS2 = STAR1 / PLUS1;
      WORK[J] = PLUS1;
      WORK[N + J] = STAR1;
      if (J + 1 <= N) {
        WORK[J + 1] = PLUS2;
        WORK[N + J + 1] = ZERO;
        PLUS1 = STAR1 / PLUS2;
        final REXP = dlarnd(2, ISEED);
        STAR1 *= pow(SFAC, REXP);
        if (REXP < ZERO) {
          STAR1 = -pow(SFAC, ONE - REXP).toDouble();
        } else {
          STAR1 = pow(SFAC, ONE + REXP).toDouble();
        }
      }
    }

    final X = sqrt(CNDNUM) - ONE / sqrt(CNDNUM);
    final Y = N > 2 ? sqrt(TWO / (N - 2)) * X : ZERO;
    final Z = X * X;

    if (UPPER) {
      // Set the upper triangle of A with a unit triangular matrix
      // of known condition number.

      var JC = 1;
      for (var J = 2; J <= N; J++) {
        A[JC + 1] = Y;
        if (J > 2) A[JC + J - 1] = WORK[J - 2];
        if (J > 3) A[JC + J - 2] = WORK[N + J - 3];
        JC += J;
      }
      JC -= N;
      A[JC + 1] = Z;
      for (var J = 2; J <= N - 1; J++) {
        A[JC + J] = Y;
      }
    } else {
      // Set the lower triangle of A with a unit triangular matrix
      // of known condition number.

      for (var I = 2; I <= N - 1; I++) {
        A[I] = Y;
      }
      A[N] = Z;
      var JC = N + 1;
      for (var J = 2; J <= N - 1; J++) {
        A[JC + 1] = WORK[J - 1];
        if (J < N - 1) A[JC + 2] = WORK[N + J - 1];
        A[JC + N - J] = Y;
        JC += N - J + 1;
      }
    }

    // Fill in the zeros using Givens rotations

    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N - 1; J++) {
        final JCNEXT = JC + J;
        final RA = Box(A[JCNEXT + J - 1]);
        final RB = Box(TWO);
        final C = Box(0.0), S = Box(0.0);
        drotg(RA, RB, C, S);

        // Multiply by [ c  s; -s  c] on the left.

        if (N > J + 1) {
          var JX = JCNEXT + J;
          for (var I = J + 2; I <= N; I++) {
            final STEMP = C.value * A[JX + J] + S.value * A[JX + J + 1];
            A[JX + J + 1] = -S.value * A[JX + J] + C.value * A[JX + J + 1];
            A[JX + J] = STEMP;
            JX += I;
          }
        }

        // Multiply by [-c -s;  s -c] on the right.

        if (J > 1) drot(J - 1, A(JCNEXT), 1, A(JC), 1, -C.value, -S.value);

        // Negate A(J,J+1).

        A[JCNEXT + J - 1] = -A[JCNEXT + J - 1];
        JC = JCNEXT;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N - 1; J++) {
        final JCNEXT = JC + N - J + 1;
        final RA = Box(A[JC + 1]);
        final RB = Box(TWO);
        final C = Box(0.0), S = Box(0.0);
        drotg(RA, RB, C, S);

        // Multiply by [ c -s;  s  c] on the right.

        if (N > J + 1) {
          drot(N - J - 1, A(JCNEXT + 1), 1, A(JC + 2), 1, C.value, -S.value);
        }

        // Multiply by [-c  s; -s -c] on the left.

        if (J > 1) {
          var JX = 1;
          for (var I = 1; I <= J - 1; I++) {
            final STEMP =
                -C.value * A[JX + J - I] + S.value * A[JX + J - I + 1];
            A[JX + J - I + 1] =
                -S.value * A[JX + J - I] - C.value * A[JX + J - I + 1];
            A[JX + J - I] = STEMP;
            JX += N - I + 1;
          }
        }

        // Negate A(J+1,J).

        A[JC + 1] = -A[JC + 1];
        JC = JCNEXT;
      }
    }

    // IMAT > 10:  Pathological test cases.  These triangular matrices
    // are badly scaled or badly conditioned, so when used in solving a
    // triangular system they may cause overflow in the solution vector.
  } else if (IMAT == 11) {
    // Type 11:  Generate a triangular matrix with elements between
    // -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
    // Make the right hand side large so that it requires scaling.

    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J, A(JC));
        A[JC + J - 1] = sign(TWO, A[JC + J - 1]);
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J + 1, A(JC));
        A[JC] = sign(TWO, A[JC]);
        JC += N - J + 1;
      }
    }

    // Set the right hand side so that the largest value is BIGNUM.

    dlarnv(2, ISEED, N, B);
    final IY = idamax(N, B, 1);
    final BNORM = B[IY].abs();
    final BSCAL = BIGNUM / max(ONE, BNORM);
    dscal(N, BSCAL, B, 1);
  } else if (IMAT == 12) {
    // Type 12:  Make the first diagonal element in the solve small to
    // cause immediate overflow when dividing by T(j,j).
    // In type 12, the offdiagonal elements are small (CNORM(j) < 1).

    dlarnv(2, ISEED, N, B);
    final TSCAL = ONE / max(ONE, N - 1);
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J - 1, A(JC));
        dscal(J - 1, TSCAL, A(JC), 1);
        A[JC + J - 1] = sign(ONE, dlarnd(2, ISEED));
        JC += J;
      }
      A[N * (N + 1) ~/ 2] = SMLNUM;
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J, A(JC + 1));
        dscal(N - J, TSCAL, A(JC + 1), 1);
        A[JC] = sign(ONE, dlarnd(2, ISEED));
        JC += N - J + 1;
      }
      A[1] = SMLNUM;
    }
  } else if (IMAT == 13) {
    // Type 13:  Make the first diagonal element in the solve small to
    // cause immediate overflow when dividing by T(j,j).
    // In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).

    dlarnv(2, ISEED, N, B);
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J - 1, A(JC));
        A[JC + J - 1] = sign(ONE, dlarnd(2, ISEED));
        JC += J;
      }
      A[N * (N + 1) ~/ 2] = SMLNUM;
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J, A(JC + 1));
        A[JC] = sign(ONE, dlarnd(2, ISEED));
        JC += N - J + 1;
      }
      A[1] = SMLNUM;
    }
  } else if (IMAT == 14) {
    // Type 14:  T is diagonal with small numbers on the diagonal to
    // make the growth factor underflow, but a small right hand side
    // chosen so that the solution does not overflow.

    if (UPPER) {
      var JCOUNT = 1;
      var JC = (N - 1) * N ~/ 2 + 1;
      for (var J = N; J >= 1; J--) {
        for (var I = 1; I <= J - 1; I++) {
          A[JC + I - 1] = ZERO;
        }
        if (JCOUNT <= 2) {
          A[JC + J - 1] = SMLNUM;
        } else {
          A[JC + J - 1] = ONE;
        }
        JCOUNT++;
        if (JCOUNT > 4) JCOUNT = 1;
        JC -= J - 1;
      }
    } else {
      var JCOUNT = 1;
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = J + 1; I <= N; I++) {
          A[JC + I - J] = ZERO;
        }
        if (JCOUNT <= 2) {
          A[JC] = SMLNUM;
        } else {
          A[JC] = ONE;
        }
        JCOUNT++;
        if (JCOUNT > 4) JCOUNT = 1;
        JC += N - J + 1;
      }
    }

    // Set the right hand side alternately zero and small.

    if (UPPER) {
      B[1] = ZERO;
      for (var I = N; I >= 2; I -= 2) {
        B[I] = ZERO;
        B[I - 1] = SMLNUM;
      }
    } else {
      B[N] = ZERO;
      for (var I = 1; I <= N - 1; I += 2) {
        B[I] = ZERO;
        B[I + 1] = SMLNUM;
      }
    }
  } else if (IMAT == 15) {
    // Type 15:  Make the diagonal elements small to cause gradual
    // overflow when dividing by T(j,j).  To control the amount of
    // scaling needed, the matrix is bidiagonal.

    final TEXP = ONE / max(ONE, N - 1);
    final TSCAL = pow(SMLNUM, TEXP).toDouble();
    dlarnv(2, ISEED, N, B);
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J - 2; I++) {
          A[JC + I - 1] = ZERO;
        }
        if (J > 1) A[JC + J - 2] = -ONE;
        A[JC + J - 1] = TSCAL;
        JC += J;
      }
      B[N] = ONE;
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = J + 2; I <= N; I++) {
          A[JC + I - J] = ZERO;
        }
        if (J < N) A[JC + 1] = -ONE;
        A[JC] = TSCAL;
        JC += N - J + 1;
      }
      B[1] = ONE;
    }
  } else if (IMAT == 16) {
    // Type 16:  One zero diagonal element.

    final IY = N ~/ 2 + 1;
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J, A(JC));
        if (J != IY) {
          A[JC + J - 1] = sign(TWO, A[JC + J - 1]);
        } else {
          A[JC + J - 1] = ZERO;
        }
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J + 1, A(JC));
        if (J != IY) {
          A[JC] = sign(TWO, A[JC]);
        } else {
          A[JC] = ZERO;
        }
        JC += N - J + 1;
      }
    }
    dlarnv(2, ISEED, N, B);
    dscal(N, TWO, B, 1);
  } else if (IMAT == 17) {
    // Type 17:  Make the offdiagonal elements large to cause overflow
    // when adding a column of T.  In the non-transposed case, the
    // matrix is constructed to cause overflow when adding a column in
    // every other step.

    final TSCAL = (ONE - ULP) / (UNFL / ULP);
    for (var J = 1; J <= N * (N + 1) ~/ 2; J++) {
      A[J] = ZERO;
    }
    var TEXP = ONE;
    if (UPPER) {
      var JC = (N - 1) * N ~/ 2 + 1;
      for (var J = N; J >= 2; J -= 2) {
        A[JC] = -TSCAL / (N + 1);
        A[JC + J - 1] = ONE;
        B[J] = TEXP * (ONE - ULP);
        JC -= J - 1;
        A[JC] = -(TSCAL / (N + 1)) / (N + 2);
        A[JC + J - 2] = ONE;
        B[J - 1] = TEXP * (N * N + N - 1);
        TEXP *= TWO;
        JC -= J - 2;
      }
      B[1] = ((N + 1) / (N + 2)) * TSCAL;
    } else {
      var JC = 1;
      for (var J = 1; J <= N - 1; J += 2) {
        A[JC + N - J] = -TSCAL / (N + 1);
        A[JC] = ONE;
        B[J] = TEXP * (ONE - ULP);
        JC += N - J + 1;
        A[JC + N - J - 1] = -(TSCAL / (N + 1)) / (N + 2);
        A[JC] = ONE;
        B[J + 1] = TEXP * (N * N + N - 1);
        TEXP *= TWO;
        JC += N - J;
      }
      B[N] = ((N + 1) / (N + 2)) * TSCAL;
    }
  } else if (IMAT == 18) {
    // Type 18:  Generate a unit triangular matrix with elements
    // between -1 and 1, and make the right hand side large so that it
    // requires scaling.

    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J - 1, A(JC));
        A[JC + J - 1] = ZERO;
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        if (J < N) dlarnv(2, ISEED, N - J, A(JC + 1));
        A[JC] = ZERO;
        JC += N - J + 1;
      }
    }

    // Set the right hand side so that the largest value is BIGNUM.

    dlarnv(2, ISEED, N, B);
    final IY = idamax(N, B, 1);
    final BNORM = B[IY].abs();
    final BSCAL = BIGNUM / max(ONE, BNORM);
    dscal(N, BSCAL, B, 1);
  } else if (IMAT == 19) {
    // Type 19:  Generate a triangular matrix with elements between
    // BIGNUM/(n-1) and BIGNUM so that at least one of the column
    // norms will exceed BIGNUM.

    final TLEFT = BIGNUM / max(ONE, N - 1);
    final TSCAL = BIGNUM * ((N - 1) / max(ONE, N));
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J, A(JC));
        for (var I = 1; I <= J; I++) {
          A[JC + I - 1] = sign(TLEFT, A[JC + I - 1]) + TSCAL * A[JC + I - 1];
        }
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J + 1, A(JC));
        for (var I = J; I <= N; I++) {
          A[JC + I - J] = sign(TLEFT, A[JC + I - J]) + TSCAL * A[JC + I - J];
        }
        JC += N - J + 1;
      }
    }
    dlarnv(2, ISEED, N, B);
    dscal(N, TWO, B, 1);
  }

  // Flip the matrix across its counter-diagonal if the transpose will
  // be used.

  if (!lsame(TRANS, 'N')) {
    if (UPPER) {
      var JJ = 1;
      var JR = N * (N + 1) ~/ 2;
      for (var J = 1; J <= N ~/ 2; J++) {
        var JL = JJ;
        for (var I = J; I <= N - J; I++) {
          final T = A[JR - I + J];
          A[JR - I + J] = A[JL];
          A[JL] = T;
          JL += I;
        }
        JJ += J + 1;
        JR -= (N - J + 1);
      }
    } else {
      var JL = 1;
      var JJ = N * (N + 1) ~/ 2;
      for (var J = 1; J <= N ~/ 2; J++) {
        var JR = JJ;
        for (var I = J; I <= N - J; I++) {
          final T = A[JL + I - J];
          A[JL + I - J] = A[JR];
          A[JR] = T;
          JR -= I;
        }
        JL += N - J + 1;
        JJ -= J + 1;
      }
    }
  }
}
