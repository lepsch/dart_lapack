import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlarnv.dart';

import '../matgen/dlarnd.dart';
import '../matgen/zlarnd.dart';
import '../matgen/zlatms.dart';
import 'zlatb4.dart';

void zlattb(
  final int IMAT,
  final String UPLO,
  final String TRANS,
  final Box<String> DIAG,
  final Array<int> ISEED_,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<Complex> B_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having(length: 4);
  final B = B_.having();
  final AB = AB_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ONE = 1.0, TWO = 2.0, ZERO = 0.0;

  final PATH = '${'Zomplex precision'[0]}TB';
  final UNFL = dlamch('Safe minimum');
  final ULP = dlamch('Epsilon') * dlamch('Base');
  final SMLNUM = UNFL;
  final BIGNUM = (ONE - ULP) / SMLNUM;
  if ((IMAT >= 6 && IMAT <= 9) || IMAT == 17) {
    DIAG.value = 'U';
  } else {
    DIAG.value = 'N';
  }
  INFO.value = 0;

  // Quick return if N <= 0.

  if (N <= 0) return;

  // Call ZLATB4 to set parameters for ZLATMS.

  final UPPER = lsame(UPLO, 'U');
  final (:TYPE, KL: _, KU: _, :ANORM, :MODE, :CNDNUM, :DIST) =
      zlatb4(PATH, UPPER ? IMAT : -IMAT, N, N);
  final (KU, KL, IOFF, PACKIT) =
      UPPER ? (KD, 0, 1 + max(0, KD - N + 1).toInt(), 'Q') : (0, KD, 1, 'B');

  // IMAT <= 5:  Non-unit triangular matrix

  if (IMAT <= 5) {
    zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT,
        AB(IOFF, 1), LDAB, WORK, INFO);

    // IMAT > 5:  Unit triangular matrix
    // The diagonal is deliberately set to something other than 1.

    // IMAT = 6:  Matrix is the identity
  } else if (IMAT == 6) {
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = max(1, KD + 2 - J); I <= KD; I++) {
          AB[I][J] = Complex.zero;
        }
        AB[KD + 1][J] = J.toComplex();
      }
    } else {
      for (var J = 1; J <= N; J++) {
        AB[1][J] = J.toComplex();
        for (var I = 2; I <= min(KD + 1, N - J + 1); I++) {
          AB[I][J] = Complex.zero;
        }
      }
    }

    // IMAT > 6:  Non-trivial unit triangular matrix

    // A unit triangular matrix T with condition CNDNUM is formed.
    // In this version, T only has bandwidth 2, the rest of it is zero.
  } else if (IMAT <= 9) {
    final TNORM = sqrt(CNDNUM);

    // Initialize AB to zero.

    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = max(1, KD + 2 - J); I <= KD; I++) {
          AB[I][J] = Complex.zero;
        }
        AB[KD + 1][J] = J.toComplex();
      }
    } else {
      for (var J = 1; J <= N; J++) {
        for (var I = 2; I <= min(KD + 1, N - J + 1); I++) {
          AB[I][J] = Complex.zero;
        }
        AB[1][J] = J.toComplex();
      }
    }

    // Special case:  T is tridiagonal.  Set every other offdiagonal
    // so that the matrix has norm TNORM+1.

    if (KD == 1) {
      if (UPPER) {
        AB[1][2] = TNORM.toComplex() * zlarnd(5, ISEED);
        final LENJ = (N - 3) ~/ 2;
        zlarnv(2, ISEED, LENJ, WORK);
        for (var J = 1; J <= LENJ; J++) {
          AB[1][2 * (J + 1)] = TNORM.toComplex() * WORK[J];
        }
      } else {
        AB[2][1] = TNORM.toComplex() * zlarnd(5, ISEED);
        final LENJ = (N - 3) ~/ 2;
        zlarnv(2, ISEED, LENJ, WORK);
        for (var J = 1; J <= LENJ; J++) {
          AB[2][2 * J + 1] = TNORM.toComplex() * WORK[J];
        }
      }
    } else if (KD > 1) {
      // Form a unit triangular matrix T with condition CNDNUM.  T is
      // given by
      //         | 1   +   *                      |
      //         |     1   +                      |
      //     T = |         1   +   *              |
      //         |             1   +              |
      //         |                 1   +   *      |
      //         |                     1   +      |
      //         |                          . . . |
      // Each element marked with a '*' is formed by taking the product
      // of the adjacent elements marked with '+'.  The '*'s can be
      // chosen freely, and the '+'s are chosen so that the inverse of
      // T will have elements of the same magnitude as T.

      // The two offdiagonals of T are stored in WORK.

      var STAR1 = TNORM.toComplex() * zlarnd(5, ISEED);
      final SFAC = sqrt(TNORM);
      var PLUS1 = SFAC.toComplex() * zlarnd(5, ISEED);
      for (var J = 1; J <= N; J += 2) {
        final PLUS2 = STAR1 / PLUS1;
        WORK[J] = PLUS1;
        WORK[N + J] = STAR1;
        if (J + 1 <= N) {
          WORK[J + 1] = PLUS2;
          WORK[N + J + 1] = Complex.zero;
          PLUS1 = STAR1 / PLUS2;

          // Generate a new *-value with norm between sqrt(TNORM)
          // and TNORM.

          final REXP = dlarnd(2, ISEED);
          if (REXP < ZERO) {
            STAR1 = pow(-SFAC, ONE - REXP).toComplex() * zlarnd(5, ISEED);
          } else {
            STAR1 = pow(SFAC, ONE + REXP).toComplex() * zlarnd(5, ISEED);
          }
        }
      }

      // Copy the tridiagonal T to AB.

      if (UPPER) {
        zcopy(N - 1, WORK, 1, AB(KD, 2).asArray(), LDAB);
        zcopy(N - 2, WORK(N + 1), 1, AB(KD - 1, 3).asArray(), LDAB);
      } else {
        zcopy(N - 1, WORK, 1, AB(2, 1).asArray(), LDAB);
        zcopy(N - 2, WORK(N + 1), 1, AB(3, 1).asArray(), LDAB);
      }
    }

    // IMAT > 9:  Pathological test cases.  These triangular matrices
    // are badly scaled or badly conditioned, so when used in solving a
    // triangular system they may cause overflow in the solution vector.
  } else if (IMAT == 10) {
    // Type 10:  Generate a triangular matrix with elements between
    // -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
    // Make the right hand side large so that it requires scaling.

    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J - 1, KD);
        zlarnv(4, ISEED, LENJ, AB(KD + 1 - LENJ, J).asArray());
        AB[KD + 1][J] = zlarnd(5, ISEED) * TWO.toComplex();
      }
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J, KD);
        if (LENJ > 0) zlarnv(4, ISEED, LENJ, AB(2, J).asArray());
        AB[1][J] = zlarnd(5, ISEED) * TWO.toComplex();
      }
    }

    // Set the right hand side so that the largest value is BIGNUM.

    zlarnv(2, ISEED, N, B);
    final IY = izamax(N, B, 1);
    final BNORM = B[IY].abs();
    final BSCAL = BIGNUM / max(ONE, BNORM);
    zdscal(N, BSCAL, B, 1);
  } else if (IMAT == 11) {
    // Type 11:  Make the first diagonal element in the solve small to
    // cause immediate overflow when dividing by T(j,j).
    // In type 11, the offdiagonal elements are small (CNORM(j) < 1).

    zlarnv(2, ISEED, N, B);
    final TSCAL = ONE / (KD + 1).toDouble();
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J - 1, KD);
        if (LENJ > 0) {
          zlarnv(4, ISEED, LENJ, AB(KD + 2 - LENJ, J).asArray());
          zdscal(LENJ, TSCAL, AB(KD + 2 - LENJ, J).asArray(), 1);
        }
        AB[KD + 1][J] = zlarnd(5, ISEED);
      }
      AB[KD + 1][N] = SMLNUM.toComplex() * AB[KD + 1][N];
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J, KD);
        if (LENJ > 0) {
          zlarnv(4, ISEED, LENJ, AB(2, J).asArray());
          zdscal(LENJ, TSCAL, AB(2, J).asArray(), 1);
        }
        AB[1][J] = zlarnd(5, ISEED);
      }
      AB[1][1] = SMLNUM.toComplex() * AB[1][1];
    }
  } else if (IMAT == 12) {
    // Type 12:  Make the first diagonal element in the solve small to
    // cause immediate overflow when dividing by T(j,j).
    // In type 12, the offdiagonal elements are O(1) (CNORM(j) > 1).

    zlarnv(2, ISEED, N, B);
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J - 1, KD);
        if (LENJ > 0) zlarnv(4, ISEED, LENJ, AB(KD + 2 - LENJ, J).asArray());
        AB[KD + 1][J] = zlarnd(5, ISEED);
      }
      AB[KD + 1][N] = SMLNUM.toComplex() * AB[KD + 1][N];
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J, KD);
        if (LENJ > 0) zlarnv(4, ISEED, LENJ, AB(2, J).asArray());
        AB[1][J] = zlarnd(5, ISEED);
      }
      AB[1][1] = SMLNUM.toComplex() * AB[1][1];
    }
  } else if (IMAT == 13) {
    // Type 13:  T is diagonal with small numbers on the diagonal to
    // make the growth factor underflow, but a small right hand side
    // chosen so that the solution does not overflow.

    if (UPPER) {
      var JCOUNT = 1;
      for (var J = N; J >= 1; J--) {
        for (var I = max(1, KD + 1 - (J - 1)); I <= KD; I++) {
          AB[I][J] = Complex.zero;
        }
        if (JCOUNT <= 2) {
          AB[KD + 1][J] = SMLNUM.toComplex() * zlarnd(5, ISEED);
        } else {
          AB[KD + 1][J] = zlarnd(5, ISEED);
        }
        JCOUNT++;
        if (JCOUNT > 4) JCOUNT = 1;
      }
    } else {
      var JCOUNT = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = 2; I <= min(N - J + 1, KD + 1); I++) {
          AB[I][J] = Complex.zero;
        }
        if (JCOUNT <= 2) {
          AB[1][J] = SMLNUM.toComplex() * zlarnd(5, ISEED);
        } else {
          AB[1][J] = zlarnd(5, ISEED);
        }
        JCOUNT++;
        if (JCOUNT > 4) JCOUNT = 1;
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
  } else if (IMAT == 14) {
    // Type 14:  Make the diagonal elements small to cause gradual
    // overflow when dividing by T(j,j).  To control the amount of
    // scaling needed, the matrix is bidiagonal.

    final TEXP = ONE / (KD + 1).toDouble();
    final TSCAL = pow(SMLNUM, TEXP).toInt();
    zlarnv(4, ISEED, N, B);
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = max(1, KD + 2 - J); I <= KD; I++) {
          AB[I][J] = Complex.zero;
        }
        if (J > 1 && KD > 0) AB[KD][J] = Complex(-ONE, -ONE);
        AB[KD + 1][J] = TSCAL.toComplex() * zlarnd(5, ISEED);
      }
      B[N] = Complex(ONE, ONE);
    } else {
      for (var J = 1; J <= N; J++) {
        for (var I = 3; I <= min(N - J + 1, KD + 1); I++) {
          AB[I][J] = Complex.zero;
        }
        if (J < N && KD > 0) AB[2][J] = Complex(-ONE, -ONE);
        AB[1][J] = TSCAL.toComplex() * zlarnd(5, ISEED);
      }
      B[1] = Complex(ONE, ONE);
    }
  } else if (IMAT == 15) {
    // Type 15:  One zero diagonal element.

    final IY = N ~/ 2 + 1;
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J, KD + 1);
        zlarnv(4, ISEED, LENJ, AB(KD + 2 - LENJ, J).asArray());
        if (J != IY) {
          AB[KD + 1][J] = zlarnd(5, ISEED) * TWO.toComplex();
        } else {
          AB[KD + 1][J] = Complex.zero;
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J + 1, KD + 1);
        zlarnv(4, ISEED, LENJ, AB(1, J).asArray());
        if (J != IY) {
          AB[1][J] = zlarnd(5, ISEED) * TWO.toComplex();
        } else {
          AB[1][J] = Complex.zero;
        }
      }
    }
    zlarnv(2, ISEED, N, B);
    zdscal(N, TWO, B, 1);
  } else if (IMAT == 16) {
    // Type 16:  Make the offdiagonal elements large to cause overflow
    // when adding a column of T.  In the non-transposed case, the
    // matrix is constructed to cause overflow when adding a column in
    // every other step.

    final TSCAL = (ONE - ULP) / (UNFL / ULP);
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= KD + 1; I++) {
        AB[I][J] = Complex.zero;
      }
    }
    var TEXP = ONE;
    if (KD > 0) {
      if (UPPER) {
        for (var J = N; J >= 1; J -= KD) {
          for (var I = J; I >= max(1, J - KD + 1); I -= 2) {
            AB[1 + (J - I)][I] = (-TSCAL / (KD + 2)).toComplex();
            AB[KD + 1][I] = Complex.one;
            B[I] = (TEXP * (ONE - ULP)).toComplex();
            if (I > max(1, J - KD + 1)) {
              AB[2 + (J - I)][I - 1] =
                  (-(TSCAL / (KD + 2)) / (KD + 3)).toComplex();
              AB[KD + 1][I - 1] = Complex.one;
              B[I - 1] = (TEXP * ((KD + 1) * (KD + 1) + KD)).toComplex();
            }
            TEXP *= TWO;
          }
          B[max(1, J - KD + 1)] = (((KD + 2) / (KD + 3)) * TSCAL).toComplex();
        }
      } else {
        for (var J = 1; J <= N; J += KD) {
          var TEXP = ONE;
          final LENJ = min(KD + 1, N - J + 1);
          for (var I = J; I <= min(N, J + KD - 1); I += 2) {
            AB[LENJ - (I - J)][J] = (-TSCAL / (KD + 2)).toComplex();
            AB[1][J] = Complex.one;
            B[J] = (TEXP * (ONE - ULP)).toComplex();
            if (I < min(N, J + KD - 1)) {
              AB[LENJ - (I - J + 1)][I + 1] =
                  (-(TSCAL / (KD + 2)) / (KD + 3)).toComplex();
              AB[1][I + 1] = Complex.one;
              B[I + 1] = (TEXP * ((KD + 1) * (KD + 1) + KD)).toComplex();
            }
            TEXP *= TWO;
          }
          B[min(N, J + KD - 1)] = (((KD + 2) / (KD + 3)) * TSCAL).toComplex();
        }
      }
    }
  } else if (IMAT == 17) {
    // Type 17:  Generate a unit triangular matrix with elements
    // between -1 and 1, and make the right hand side large so that it
    // requires scaling.

    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J - 1, KD);
        zlarnv(4, ISEED, LENJ, AB(KD + 1 - LENJ, J).asArray());
        AB[KD + 1][J] = J.toComplex();
      }
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J, KD);
        if (LENJ > 0) zlarnv(4, ISEED, LENJ, AB(2, J).asArray());
        AB[1][J] = J.toComplex();
      }
    }

    // Set the right hand side so that the largest value is BIGNUM.

    zlarnv(2, ISEED, N, B);
    final IY = izamax(N, B, 1);
    final BNORM = B[IY].abs();
    final BSCAL = BIGNUM / max(ONE, BNORM);
    zdscal(N, BSCAL, B, 1);
  } else if (IMAT == 18) {
    // Type 18:  Generate a triangular matrix with elements between
    // BIGNUM/(KD+1) and BIGNUM so that at least one of the column
    // norms will exceed BIGNUM.
    // 1/3/91:  ZLATBS no longer can handle this case

    final TLEFT = BIGNUM / (KD + 1);
    final TSCAL = BIGNUM * (KD + 1) / (KD + 2);
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J, KD + 1);
        zlarnv(5, ISEED, LENJ, AB(KD + 2 - LENJ, J).asArray());
        dlarnv(1, ISEED, LENJ, RWORK(KD + 2 - LENJ));
        for (var I = KD + 2 - LENJ; I <= KD + 1; I++) {
          AB[I][J] *= (TLEFT + RWORK[I] * TSCAL).toComplex();
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J + 1, KD + 1);
        zlarnv(5, ISEED, LENJ, AB(1, J).asArray());
        dlarnv(1, ISEED, LENJ, RWORK);
        for (var I = 1; I <= LENJ; I++) {
          AB[I][J] *= (TLEFT + RWORK[I] * TSCAL).toComplex();
        }
      }
    }
    zlarnv(2, ISEED, N, B);
    zdscal(N, TWO, B, 1);
  }

  // Flip the matrix if the transpose will be used.

  if (!lsame(TRANS, 'N')) {
    if (UPPER) {
      for (var J = 1; J <= N ~/ 2; J++) {
        final LENJ = min(N - 2 * J + 1, KD + 1);
        zswap(LENJ, AB(KD + 1, J).asArray(), LDAB - 1,
            AB(KD + 2 - LENJ, N - J + 1).asArray(), -1);
      }
    } else {
      for (var J = 1; J <= N ~/ 2; J++) {
        final LENJ = min(N - 2 * J + 1, KD + 1);
        zswap(LENJ, AB(1, J).asArray(), 1, AB(LENJ, N - J + 2 - LENJ).asArray(),
            -LDAB + 1);
      }
    }
  }
}
