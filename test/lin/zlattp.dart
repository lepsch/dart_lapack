// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zrotg.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlarnv.dart';
import 'package:lapack/src/zrot.dart';

import '../matgen/zlarnd.dart';
import '../matgen/zlatms.dart';
import 'zlatb4.dart';

void zlattp(
  final int IMAT,
  final String UPLO,
  final String TRANS,
  final Box<String> DIAG,
  final Array<int> ISEED_,
  final int N,
  final Array<Complex> AP_,
  final Array<Complex> B_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final B = B_.having();
  final ISEED = ISEED_.having(length: 4);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ONE = 1.0, TWO = 2.0, ZERO = 0.0;

  final PATH = '${'Zomplex precision'[0]}TP';
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

  // Call ZLATB4 to set parameters for CLATMS.

  final UPPER = lsame(UPLO, 'U');
  final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
      zlatb4(PATH, UPPER ? IMAT : -IMAT, N, N);
  final PACKIT = UPPER ? 'C' : 'R';

  // IMAT <= 6:  Non-unit triangular matrix

  if (IMAT <= 6) {
    zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT,
        AP.asMatrix(), N, WORK, INFO);

    // IMAT > 6:  Unit triangular matrix
    // The diagonal is deliberately set to something other than 1.

    // IMAT = 7:  Matrix is the identity
  } else if (IMAT == 7) {
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J - 1; I++) {
          AP[JC + I - 1] = Complex.zero;
        }
        AP[JC + J - 1] = J.toComplex();
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        AP[JC] = J.toComplex();
        for (var I = J + 1; I <= N; I++) {
          AP[JC + I - J] = Complex.zero;
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
          AP[JC + I] = Complex.zero;
        }
        AP[JC + J] = J.toComplex();
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        AP[JC] = J.toComplex();
        for (var I = J + 1; I <= N; I++) {
          AP[JC + I - J] = Complex.zero;
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

    //       [  c   s ] [  1   w ] [ -c  -s ] =  [  1  -w ]
    //       [ -s   c ] [  0   1 ] [  s  -c ]    [  0   1 ]
    // and
    //       [ -c  -s ] [  1   0 ] [  c   s ] =  [  1   0 ]
    //       [  s  -c ] [  w   1 ] [ -s   c ]    [ -w   1 ]

    // where c = w / sqrt(w**2+4) and s = 2 / sqrt(w**2+4).

    var STAR1 = 0.25.toComplex() * zlarnd(5, ISEED);
    final SFAC = 0.5;
    var PLUS1 = SFAC.toComplex() * zlarnd(5, ISEED);
    for (var J = 1; J <= N; J += 2) {
      final PLUS2 = STAR1 / PLUS1;
      WORK[J] = PLUS1;
      WORK[N + J] = STAR1;
      if (J + 1 <= N) {
        WORK[J + 1] = PLUS2;
        WORK[N + J + 1] = Complex.zero;
        PLUS1 = STAR1 / PLUS2;
        final REXP = zlarnd(2, ISEED).real;
        if (REXP < ZERO) {
          STAR1 = -pow(SFAC, ONE - REXP).toComplex() * zlarnd(5, ISEED);
        } else {
          STAR1 = pow(SFAC, ONE + REXP).toComplex() * zlarnd(5, ISEED);
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
        AP[JC + 1] = Y.toComplex();
        if (J > 2) AP[JC + J - 1] = WORK[J - 2];
        if (J > 3) AP[JC + J - 2] = WORK[N + J - 3];
        JC += J;
      }
      JC -= N;
      AP[JC + 1] = Z.toComplex();
      for (var J = 2; J <= N - 1; J++) {
        AP[JC + J] = Y.toComplex();
      }
    } else {
      // Set the lower triangle of A with a unit triangular matrix
      // of known condition number.

      for (var I = 2; I <= N - 1; I++) {
        AP[I] = Y.toComplex();
      }
      AP[N] = Z.toComplex();
      var JC = N + 1;
      for (var J = 2; J <= N - 1; J++) {
        AP[JC + 1] = WORK[J - 1];
        if (J < N - 1) AP[JC + 2] = WORK[N + J - 1];
        AP[JC + N - J] = Y.toComplex();
        JC += N - J + 1;
      }
    }

    // Fill in the zeros using Givens rotations

    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N - 1; J++) {
        final JCNEXT = JC + J;
        final RA = Box(AP[JCNEXT + J - 1]);
        final RB = TWO.toComplex();
        final C = Box(ZERO);
        final S = Box(Complex.zero);
        zrotg(RA, RB, C, S);

        // Multiply by [ c  s; -conjg(s)  c] on the left.

        if (N > J + 1) {
          var JX = JCNEXT + J;
          for (var I = J + 2; I <= N; I++) {
            final CTEMP =
                C.value.toComplex() * AP[JX + J] + S.value * AP[JX + J + 1];
            AP[JX + J + 1] = -S.value.conjugate() * AP[JX + J] +
                C.value.toComplex() * AP[JX + J + 1];
            AP[JX + J] = CTEMP;
            JX += I;
          }
        }

        // Multiply by [-c -s;  conjg(s) -c] on the right.

        if (J > 1) zrot(J - 1, AP(JCNEXT), 1, AP(JC), 1, -C.value, -S.value);

        // Negate A(J,J+1).

        AP[JCNEXT + J - 1] = -AP[JCNEXT + J - 1];
        JC = JCNEXT;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N - 1; J++) {
        final JCNEXT = JC + N - J + 1;
        final RA = Box(AP[JC + 1]);
        final RB = TWO.toComplex();
        final C = Box(ZERO);
        final S = Box(Complex.zero);
        zrotg(RA, RB, C, S);
        S.value = S.value.conjugate();

        // Multiply by [ c -s;  conjg(s) c] on the right.

        if (N > J + 1) {
          zrot(N - J - 1, AP(JCNEXT + 1), 1, AP(JC + 2), 1, C.value, -S.value);
        }

        // Multiply by [-c  s; -conjg(s) -c] on the left.

        if (J > 1) {
          var JX = 1;
          for (var I = 1; I <= J - 1; I++) {
            final CTEMP = -C.value.toComplex() * AP[JX + J - I] +
                S.value * AP[JX + J - I + 1];
            AP[JX + J - I + 1] = -S.value.conjugate() * AP[JX + J - I] -
                C.value.toComplex() * AP[JX + J - I + 1];
            AP[JX + J - I] = CTEMP;
            JX += N - I + 1;
          }
        }

        // Negate A(J+1,J).

        AP[JC + 1] = -AP[JC + 1];
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
        zlarnv(4, ISEED, J - 1, AP(JC));
        AP[JC + J - 1] = zlarnd(5, ISEED) * TWO.toComplex();
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        if (J < N) zlarnv(4, ISEED, N - J, AP(JC + 1));
        AP[JC] = zlarnd(5, ISEED) * TWO.toComplex();
        JC += N - J + 1;
      }
    }

    // Set the right hand side so that the largest value is BIGNUM.

    zlarnv(2, ISEED, N, B);
    final IY = izamax(N, B, 1);
    final BNORM = B[IY].abs();
    final BSCAL = BIGNUM / max(ONE, BNORM);
    zdscal(N, BSCAL, B, 1);
  } else if (IMAT == 12) {
    // Type 12:  Make the first diagonal element in the solve small to
    // cause immediate overflow when dividing by T(j,j).
    // In type 12, the offdiagonal elements are small (CNORM(j) < 1).

    zlarnv(2, ISEED, N, B);
    final TSCAL = ONE / max(ONE, N - 1);
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        zlarnv(4, ISEED, J - 1, AP(JC));
        zdscal(J - 1, TSCAL, AP(JC), 1);
        AP[JC + J - 1] = zlarnd(5, ISEED);
        JC += J;
      }
      AP[N * (N + 1) ~/ 2] = SMLNUM.toComplex() * AP[N * (N + 1) ~/ 2];
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        zlarnv(2, ISEED, N - J, AP(JC + 1));
        zdscal(N - J, TSCAL, AP(JC + 1), 1);
        AP[JC] = zlarnd(5, ISEED);
        JC += N - J + 1;
      }
      AP[1] = SMLNUM.toComplex() * AP[1];
    }
  } else if (IMAT == 13) {
    // Type 13:  Make the first diagonal element in the solve small to
    // cause immediate overflow when dividing by T(j,j).
    // In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).

    zlarnv(2, ISEED, N, B);
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        zlarnv(4, ISEED, J - 1, AP(JC));
        AP[JC + J - 1] = zlarnd(5, ISEED);
        JC += J;
      }
      AP[N * (N + 1) ~/ 2] = SMLNUM.toComplex() * AP[N * (N + 1) ~/ 2];
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        zlarnv(4, ISEED, N - J, AP(JC + 1));
        AP[JC] = zlarnd(5, ISEED);
        JC += N - J + 1;
      }
      AP[1] = SMLNUM.toComplex() * AP[1];
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
          AP[JC + I - 1] = Complex.zero;
        }
        if (JCOUNT <= 2) {
          AP[JC + J - 1] = SMLNUM.toComplex() * zlarnd(5, ISEED);
        } else {
          AP[JC + J - 1] = zlarnd(5, ISEED);
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
          AP[JC + I - J] = Complex.zero;
        }
        if (JCOUNT <= 2) {
          AP[JC] = SMLNUM.toComplex() * zlarnd(5, ISEED);
        } else {
          AP[JC] = zlarnd(5, ISEED);
        }
        JCOUNT++;
        if (JCOUNT > 4) JCOUNT = 1;
        JC += N - J + 1;
      }
    }

    // Set the right hand side alternately zero and small.

    if (UPPER) {
      B[1] = Complex.zero;
      for (var I = N; I >= 2; I -= 2) {
        B[I] = Complex.zero;
        B[I - 1] = SMLNUM.toComplex() * zlarnd(5, ISEED);
      }
    } else {
      B[N] = Complex.zero;
      for (var I = 1; I <= N - 1; I += 2) {
        B[I] = Complex.zero;
        B[I + 1] = SMLNUM.toComplex() * zlarnd(5, ISEED);
      }
    }
  } else if (IMAT == 15) {
    // Type 15:  Make the diagonal elements small to cause gradual
    // overflow when dividing by T(j,j).  To control the amount of
    // scaling needed, the matrix is bidiagonal.

    final TEXP = ONE / max(ONE, N - 1);
    final TSCAL = pow(SMLNUM, TEXP).toDouble();
    zlarnv(4, ISEED, N, B);
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J - 2; I++) {
          AP[JC + I - 1] = Complex.zero;
        }
        if (J > 1) AP[JC + J - 2] = Complex(-ONE, -ONE);
        AP[JC + J - 1] = TSCAL.toComplex() * zlarnd(5, ISEED);
        JC += J;
      }
      B[N] = Complex(ONE, ONE);
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = J + 2; I <= N; I++) {
          AP[JC + I - J] = Complex.zero;
        }
        if (J < N) AP[JC + 1] = Complex(-ONE, -ONE);
        AP[JC] = TSCAL.toComplex() * zlarnd(5, ISEED);
        JC += N - J + 1;
      }
      B[1] = Complex(ONE, ONE);
    }
  } else if (IMAT == 16) {
    // Type 16:  One zero diagonal element.

    final IY = N ~/ 2 + 1;
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        zlarnv(4, ISEED, J, AP(JC));
        if (J != IY) {
          AP[JC + J - 1] = zlarnd(5, ISEED) * TWO.toComplex();
        } else {
          AP[JC + J - 1] = Complex.zero;
        }
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        zlarnv(4, ISEED, N - J + 1, AP(JC));
        if (J != IY) {
          AP[JC] = zlarnd(5, ISEED) * TWO.toComplex();
        } else {
          AP[JC] = Complex.zero;
        }
        JC += N - J + 1;
      }
    }
    zlarnv(2, ISEED, N, B);
    zdscal(N, TWO, B, 1);
  } else if (IMAT == 17) {
    // Type 17:  Make the offdiagonal elements large to cause overflow
    // when adding a column of T.  In the non-transposed case, the
    // matrix is constructed to cause overflow when adding a column in
    // every other step.

    final TSCAL = (ONE - ULP) / (UNFL / ULP);
    for (var J = 1; J <= N * (N + 1) ~/ 2; J++) {
      AP[J] = Complex.zero;
    }
    var TEXP = ONE;
    if (UPPER) {
      var JC = (N - 1) * N ~/ 2 + 1;
      for (var J = N; J >= 2; J -= 2) {
        AP[JC] = (-TSCAL / (N + 1)).toComplex();
        AP[JC + J - 1] = Complex.one;
        B[J] = (TEXP * (ONE - ULP)).toComplex();
        JC -= J - 1;
        AP[JC] = (-(TSCAL / (N + 1)) / (N + 2)).toComplex();
        AP[JC + J - 2] = Complex.one;
        B[J - 1] = (TEXP * (N * N + N - 1)).toComplex();
        TEXP *= TWO;
        JC -= J - 2;
      }
      B[1] = (((N + 1) / (N + 2)) * TSCAL).toComplex();
    } else {
      var JC = 1;
      for (var J = 1; J <= N - 1; J += 2) {
        AP[JC + N - J] = (-TSCAL / (N + 1)).toComplex();
        AP[JC] = Complex.one;
        B[J] = (TEXP * (ONE - ULP)).toComplex();
        JC += N - J + 1;
        AP[JC + N - J - 1] = (-(TSCAL / (N + 1)) / (N + 2)).toComplex();
        AP[JC] = Complex.one;
        B[J + 1] = (TEXP * (N * N + N - 1)).toComplex();
        TEXP *= TWO;
        JC += N - J;
      }
      B[N] = (((N + 1) / (N + 2)) * TSCAL).toComplex();
    }
  } else if (IMAT == 18) {
    // Type 18:  Generate a unit triangular matrix with elements
    // between -1 and 1, and make the right hand side large so that it
    // requires scaling.

    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        zlarnv(4, ISEED, J - 1, AP(JC));
        AP[JC + J - 1] = Complex.zero;
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        if (J < N) zlarnv(4, ISEED, N - J, AP(JC + 1));
        AP[JC] = Complex.zero;
        JC += N - J + 1;
      }
    }

    // Set the right hand side so that the largest value is BIGNUM.

    zlarnv(2, ISEED, N, B);
    final IY = izamax(N, B, 1);
    final BNORM = B[IY].abs();
    final BSCAL = BIGNUM / max(ONE, BNORM);
    zdscal(N, BSCAL, B, 1);
  } else if (IMAT == 19) {
    // Type 19:  Generate a triangular matrix with elements between
    // BIGNUM/(n-1) and BIGNUM so that at least one of the column
    // norms will exceed BIGNUM.
    // 1/3/91:  ZLATPS no longer can handle this case

    final TLEFT = BIGNUM / max(ONE, N - 1);
    final TSCAL = BIGNUM * ((N - 1) / max(ONE, N));
    if (UPPER) {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        zlarnv(5, ISEED, J, AP(JC));
        dlarnv(1, ISEED, J, RWORK);
        for (var I = 1; I <= J; I++) {
          AP[JC + I - 1] =
              AP[JC + I - 1] * (TLEFT + RWORK[I] * TSCAL).toComplex();
        }
        JC += J;
      }
    } else {
      var JC = 1;
      for (var J = 1; J <= N; J++) {
        zlarnv(5, ISEED, N - J + 1, AP(JC));
        dlarnv(1, ISEED, N - J + 1, RWORK);
        for (var I = J; I <= N; I++) {
          AP[JC + I - J] =
              AP[JC + I - J] * (TLEFT + RWORK[I - J + 1] * TSCAL).toComplex();
        }
        JC += N - J + 1;
      }
    }
    zlarnv(2, ISEED, N, B);
    zdscal(N, TWO, B, 1);
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
          final T = AP[JR - I + J].real;
          AP[JR - I + J] = AP[JL];
          AP[JL] = T.toComplex();
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
          final T = AP[JL + I - J].real;
          AP[JL + I - J] = AP[JR];
          AP[JR] = T.toComplex();
          JR -= I;
        }
        JL += N - J + 1;
        JJ -= J + 1;
      }
    }
  }
}
