import 'dart:math' hide pow;

import 'package:lapack/lapack.dart';

import '../matgen/dlarnd.dart';
import '../matgen/dlatms.dart';
import 'dlatb4.dart';

void dlattr(
  final int IMAT,
  final String UPLO,
  final String TRANS,
  final Box<String> DIAG,
  final Array<int> ISEED_,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> B_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final B = B_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, TWO = 2.0, ZERO = 0.0;

  final PATH = '${'Double precision'[0]}TR';
  final UNFL = dlamch('Safe minimum');
  final ULP = dlamch('Epsilon') * dlamch('Base');
  final SMLNUM = UNFL;
  final BIGNUM = (ONE - ULP) / SMLNUM;
  DIAG.value = (IMAT >= 7 && IMAT <= 10) || IMAT == 18 ? 'U' : 'N';
  INFO.value = 0;

  // Quick return if N <= 0.
  if (N <= 0) return;

  // Call DLATB4 to set parameters for DLATMS.
  final UPPER = lsame(UPLO, 'U');
  final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
      dlatb4(PATH, UPPER ? IMAT : -IMAT, N, N);

  // IMAT <= 6:  Non-unit triangular matrix
  if (IMAT <= 6) {
    dlatms(N, N, DIST, ISEED, TYPE, B, MODE, CNDNUM, ANORM, KL, KU,
        'No packing', A, LDA, WORK, INFO);

    // IMAT > 6:  Unit triangular matrix
    // The diagonal is deliberately set to something other than 1.

    // IMAT = 7:  Matrix is the identity
  } else if (IMAT == 7) {
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J - 1; I++) {
          A[I][J] = ZERO;
        }
        A[J][J] = J.toDouble();
      }
    } else {
      for (var J = 1; J <= N; J++) {
        A[J][J] = J.toDouble();
        for (var I = J + 1; I <= N; I++) {
          A[I][J] = ZERO;
        }
      }
    }

    // IMAT > 7:  Non-trivial unit triangular matrix

    // Generate a unit triangular matrix T with condition CNDNUM by
    // forming a triangular matrix with known singular values and
    // filling in the zero entries with Givens rotations.
  } else if (IMAT <= 10) {
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J - 1; I++) {
          A[I][J] = ZERO;
        }
        A[J][J] = J.toDouble();
      }
    } else {
      for (var J = 1; J <= N; J++) {
        A[J][J] = J.toDouble();
        for (var I = J + 1; I <= N; I++) {
          A[I][J] = ZERO;
        }
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

    final X = sqrt(CNDNUM) - 1 / sqrt(CNDNUM);
    final Y = N > 2 ? sqrt(2 / (N - 2)) * X : ZERO;
    final Z = X * X;

    if (UPPER) {
      if (N > 3) {
        dcopy(N - 3, WORK, 1, A(2, 3).asArray(), LDA + 1);
        if (N > 4) dcopy(N - 4, WORK(N + 1), 1, A(2, 4).asArray(), LDA + 1);
      }
      for (var J = 2; J <= N - 1; J++) {
        A[1][J] = Y;
        A[J][N] = Y;
      }
      A[1][N] = Z;
    } else {
      if (N > 3) {
        dcopy(N - 3, WORK, 1, A(3, 2).asArray(), LDA + 1);
        if (N > 4) dcopy(N - 4, WORK(N + 1), 1, A(4, 2).asArray(), LDA + 1);
      }
      for (var J = 2; J <= N - 1; J++) {
        A[J][1] = Y;
        A[N][J] = Y;
      }
      A[N][1] = Z;
    }

    // Fill in the zeros using Givens rotations.
    if (UPPER) {
      for (var J = 1; J <= N - 1; J++) {
        final RA = Box(A[J][J + 1]);
        final RB = Box(2.0), C = Box(ZERO), S = Box(ZERO);
        drotg(RA, RB, C, S);

        // Multiply by [ c  s; -s  c] on the left.
        if (N > J + 1) {
          drot(N - J - 1, A(J, J + 2).asArray(), LDA, A(J + 1, J + 2).asArray(),
              LDA, C.value, S.value);
        }

        // Multiply by [-c -s;  s -c] on the right.
        if (J > 1) {
          drot(J - 1, A(1, J + 1).asArray(), 1, A(1, J).asArray(), 1, -C.value,
              -S.value);
        }

        // Negate A(J,J+1).
        A[J][J + 1] = -A[J][J + 1];
      }
    } else {
      for (var J = 1; J <= N - 1; J++) {
        final RA = Box(A[J + 1][J]);
        final RB = Box(2.0);
        final C = Box(ZERO), S = Box(ZERO);
        drotg(RA, RB, C, S);

        // Multiply by [ c -s;  s  c] on the right.
        if (N > J + 1) {
          drot(N - J - 1, A(J + 2, J + 1).asArray(), 1, A(J + 2, J).asArray(),
              1, C.value, -S.value);
        }

        // Multiply by [-c  s; -s -c] on the left.
        if (J > 1) {
          drot(J - 1, A(J, 1).asArray(), LDA, A(J + 1, 1).asArray(), LDA,
              -C.value, S.value);
        }

        // Negate A(J+1,J).
        A[J + 1][J] = -A[J + 1][J];
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
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J, A(1, J).asArray());
        A[J][J] = sign(TWO, A[J][J]);
      }
    } else {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J + 1, A(J, J).asArray());
        A[J][J] = sign(TWO, A[J][J]);
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
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J, A(1, J).asArray());
        dscal(J - 1, TSCAL, A(1, J).asArray(), 1);
        A[J][J] = sign(ONE, A[J][J]);
      }
      A[N][N] = SMLNUM * A[N][N];
    } else {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J + 1, A(J, J).asArray());
        if (N > J) dscal(N - J, TSCAL, A(J + 1, J).asArray(), 1);
        A[J][J] = sign(ONE, A[J][J]);
      }
      A[1][1] = SMLNUM * A[1][1];
    }
  } else if (IMAT == 13) {
    // Type 13:  Make the first diagonal element in the solve small to
    // cause immediate overflow when dividing by T(j,j).
    // In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).

    dlarnv(2, ISEED, N, B);
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J, A(1, J).asArray());
        A[J][J] = sign(ONE, A[J][J]);
      }
      A[N][N] = SMLNUM * A[N][N];
    } else {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J + 1, A(J, J).asArray());
        A[J][J] = sign(ONE, A[J][J]);
      }
      A[1][1] = SMLNUM * A[1][1];
    }
  } else if (IMAT == 14) {
    // Type 14:  T is diagonal with small numbers on the diagonal to
    // make the growth factor underflow, but a small right hand side
    // chosen so that the solution does not overflow.

    if (UPPER) {
      var JCOUNT = 1;
      for (var J = N; J >= 1; J--) {
        for (var I = 1; I <= J - 1; I++) {
          A[I][J] = ZERO;
        }
        if (JCOUNT <= 2) {
          A[J][J] = SMLNUM;
        } else {
          A[J][J] = ONE;
        }
        JCOUNT++;
        if (JCOUNT > 4) JCOUNT = 1;
      }
    } else {
      var JCOUNT = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = J + 1; I <= N; I++) {
          A[I][J] = ZERO;
        }
        if (JCOUNT <= 2) {
          A[J][J] = SMLNUM;
        } else {
          A[J][J] = ONE;
        }
        JCOUNT++;
        if (JCOUNT > 4) JCOUNT = 1;
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
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J - 2; I++) {
          A[I][J] = 0.0;
        }
        if (J > 1) A[J - 1][J] = -ONE;
        A[J][J] = TSCAL;
      }
      B[N] = ONE;
    } else {
      for (var J = 1; J <= N; J++) {
        for (var I = J + 2; I <= N; I++) {
          A[I][J] = 0.0;
        }
        if (J < N) A[J + 1][J] = -ONE;
        A[J][J] = TSCAL;
      }
      B[1] = ONE;
    }
  } else if (IMAT == 16) {
    // Type 16:  One zero diagonal element.

    final IY = N ~/ 2 + 1;
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J, A(1, J).asArray());
        if (J != IY) {
          A[J][J] = sign(TWO, A[J][J]);
        } else {
          A[J][J] = ZERO;
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J + 1, A(J, J).asArray());
        if (J != IY) {
          A[J][J] = sign(TWO, A[J][J]);
        } else {
          A[J][J] = ZERO;
        }
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
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= N; I++) {
        A[I][J] = 0.0;
      }
    }
    var TEXP = ONE;
    if (UPPER) {
      for (var J = N; J >= 2; J -= 2) {
        A[1][J] = -TSCAL / (N + 1);
        A[J][J] = ONE;
        B[J] = TEXP * (ONE - ULP);
        A[1][J - 1] = -(TSCAL / (N + 1)) / (N + 2);
        A[J - 1][J - 1] = ONE;
        B[J - 1] = TEXP * (N * N + N - 1);
        TEXP *= 2.0;
      }
      B[1] = ((N + 1) / (N + 2)) * TSCAL;
    } else {
      for (var J = 1; J <= N - 1; J += 2) {
        A[N][J] = -TSCAL / (N + 1);
        A[J][J] = ONE;
        B[J] = TEXP * (ONE - ULP);
        A[N][J + 1] = -(TSCAL / (N + 1)) / (N + 2);
        A[J + 1][J + 1] = ONE;
        B[J + 1] = TEXP * (N * N + N - 1);
        TEXP *= 2.0;
      }
      B[N] = ((N + 1) / (N + 2)) * TSCAL;
    }
  } else if (IMAT == 18) {
    // Type 18:  Generate a unit triangular matrix with elements
    // between -1 and 1, and make the right hand side large so that it
    // requires scaling.

    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J - 1, A(1, J).asArray());
        A[J][J] = ZERO;
      }
    } else {
      for (var J = 1; J <= N; J++) {
        if (J < N) dlarnv(2, ISEED, N - J, A(J + 1, J).asArray());
        A[J][J] = ZERO;
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
    // 1/3/91:  DLATRS no longer can handle this case

    final TLEFT = BIGNUM / max(ONE, N - 1);
    final TSCAL = BIGNUM * ((N - 1) / max(ONE, N));
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, J, A(1, J).asArray());
        for (var I = 1; I <= J; I++) {
          A[I][J] = sign(TLEFT, A[I][J]) + TSCAL * A[I][J];
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, N - J + 1, A(J, J).asArray());
        for (var I = J; I <= N; I++) {
          A[I][J] = sign(TLEFT, A[I][J]) + TSCAL * A[I][J];
        }
      }
    }
    dlarnv(2, ISEED, N, B);
    dscal(N, TWO, B, 1);
  }

  // Flip the matrix if the transpose will be used.

  if (!lsame(TRANS, 'N')) {
    if (UPPER) {
      for (var J = 1; J <= N ~/ 2; J++) {
        dswap(N - 2 * J + 1, A(J, J).asArray(), LDA,
            A(J + 1, N - J + 1).asArray(), -1);
      }
    } else {
      for (var J = 1; J <= N ~/ 2; J++) {
        dswap(N - 2 * J + 1, A(J, J).asArray(), 1,
            A(N - J + 1, J + 1).asArray(), -LDA);
      }
    }
  }
}
