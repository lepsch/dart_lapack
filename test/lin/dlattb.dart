import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';

import '../matgen/dlarnd.dart';
import '../matgen/dlatms.dart';
import 'dlatb4.dart';

void dlattb(
  final int IMAT,
  final String UPLO,
  final String TRANS,
  final Box<String> DIAG,
  final Array<int> ISEED_,
  final int N,
  final int KD,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<double> B_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having(length: 4);
  final AB = AB_.having(ld: LDAB);
  final B = B_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, TWO = 2.0, ZERO = 0.0;

  final PATH = '${'Double precision'[0]}TB';
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

  // Call DLATB4 to set parameters for DLATMS.

  final UPPER = lsame(UPLO, 'U');
  final (:TYPE, KL: _, KU: _, :ANORM, :MODE, COND: CNDNUM, :DIST) =
      dlatb4(PATH, UPPER ? IMAT : -IMAT, N, N);
  final int IOFF, KL, KU;
  final String PACKIT;
  if (UPPER) {
    KU = KD;
    IOFF = 1 + max(0, KD - N + 1);
    KL = 0;
    PACKIT = 'Q';
  } else {
    KL = KD;
    IOFF = 1;
    KU = 0;
    PACKIT = 'B';
  }

  // IMAT <= 5:  Non-unit triangular matrix

  if (IMAT <= 5) {
    dlatms(N, N, DIST, ISEED, TYPE, B, MODE, CNDNUM, ANORM, KL, KU, PACKIT,
        AB(IOFF, 1), LDAB, WORK, INFO);

    // IMAT > 5:  Unit triangular matrix
    // The diagonal is deliberately set to something other than 1.

    // IMAT = 6:  Matrix is the identity
  } else if (IMAT == 6) {
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = max(1, KD + 2 - J); I <= KD; I++) {
          AB[I][J] = ZERO;
        }
        AB[KD + 1][J] = J.toDouble();
      }
    } else {
      for (var J = 1; J <= N; J++) {
        AB[1][J] = J.toDouble();
        for (var I = 2; I <= min(KD + 1, N - J + 1); I++) {
          AB[I][J] = ZERO;
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
          AB[I][J] = ZERO;
        }
        AB[KD + 1][J] = J.toDouble();
      }
    } else {
      for (var J = 1; J <= N; J++) {
        for (var I = 2; I <= min(KD + 1, N - J + 1); I++) {
          AB[I][J] = ZERO;
        }
        AB[1][J] = J.toDouble();
      }
    }

    // Special case:  T is tridiagonal.  Set every other offdiagonal
    // so that the matrix has norm TNORM+1.

    if (KD == 1) {
      if (UPPER) {
        AB[1][2] = sign(TNORM, dlarnd(2, ISEED));
        final LENJ = (N - 3) ~/ 2;
        dlarnv(2, ISEED, LENJ, WORK);
        for (var J = 1; J <= LENJ; J++) {
          AB[1][2 * (J + 1)] = TNORM * WORK[J];
        }
      } else {
        AB[2][1] = sign(TNORM, dlarnd(2, ISEED));
        final LENJ = (N - 3) ~/ 2;
        dlarnv(2, ISEED, LENJ, WORK);
        for (var J = 1; J <= LENJ; J++) {
          AB[2][2 * J + 1] = TNORM * WORK[J];
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

      var STAR1 = sign(TNORM, dlarnd(2, ISEED));
      final SFAC = sqrt(TNORM);
      var PLUS1 = sign(SFAC, dlarnd(2, ISEED));
      for (var J = 1; J <= N; J += 2) {
        final PLUS2 = STAR1 / PLUS1;
        WORK[J] = PLUS1;
        WORK[N + J] = STAR1;
        if (J + 1 <= N) {
          WORK[J + 1] = PLUS2;
          WORK[N + J + 1] = ZERO;
          PLUS1 = STAR1 / PLUS2;

          // Generate a new *-value with norm between sqrt(TNORM)
          // and TNORM.

          final REXP = dlarnd(2, ISEED);
          if (REXP < ZERO) {
            STAR1 = pow(-SFAC, ONE - REXP).toDouble();
          } else {
            STAR1 = pow(SFAC, ONE + REXP).toDouble();
          }
        }
      }

      // Copy the tridiagonal T to AB.

      if (UPPER) {
        dcopy(N - 1, WORK, 1, AB(KD, 2).asArray(), LDAB);
        dcopy(N - 2, WORK(N + 1), 1, AB(KD - 1, 3).asArray(), LDAB);
      } else {
        dcopy(N - 1, WORK, 1, AB(2, 1).asArray(), LDAB);
        dcopy(N - 2, WORK(N + 1), 1, AB(3, 1).asArray(), LDAB);
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
        final LENJ = min(J, KD + 1);
        dlarnv(2, ISEED, LENJ, AB(KD + 2 - LENJ, J).asArray());
        AB[KD + 1][J] = sign(TWO, AB[KD + 1][J]);
      }
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J + 1, KD + 1);
        if (LENJ > 0) dlarnv(2, ISEED, LENJ, AB(1, J).asArray());
        AB[1][J] = sign(TWO, AB[1][J]);
      }
    }

    // Set the right hand side so that the largest value is BIGNUM.

    dlarnv(2, ISEED, N, B);
    final IY = idamax(N, B, 1);
    final BNORM = B[IY].abs();
    final BSCAL = BIGNUM / max(ONE, BNORM);
    dscal(N, BSCAL, B, 1);
  } else if (IMAT == 11) {
    // Type 11:  Make the first diagonal element in the solve small to
    // cause immediate overflow when dividing by T(j,j).
    // In type 11, the offdiagonal elements are small (CNORM(j) < 1).

    dlarnv(2, ISEED, N, B);
    final TSCAL = ONE / (KD + 1);
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J, KD + 1);
        dlarnv(2, ISEED, LENJ, AB(KD + 2 - LENJ, J).asArray());
        dscal(LENJ - 1, TSCAL, AB(KD + 2 - LENJ, J).asArray(), 1);
        AB[KD + 1][J] = sign(ONE, AB[KD + 1][J]);
      }
      AB[KD + 1][N] = SMLNUM * AB[KD + 1][N];
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J + 1, KD + 1);
        dlarnv(2, ISEED, LENJ, AB(1, J).asArray());
        if (LENJ > 1) dscal(LENJ - 1, TSCAL, AB(2, J).asArray(), 1);
        AB[1][J] = sign(ONE, AB[1][J]);
      }
      AB[1][1] = SMLNUM * AB[1][1];
    }
  } else if (IMAT == 12) {
    // Type 12:  Make the first diagonal element in the solve small to
    // cause immediate overflow when dividing by T(j,j).
    // In type 12, the offdiagonal elements are O(1) (CNORM(j) > 1).

    dlarnv(2, ISEED, N, B);
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J, KD + 1);
        dlarnv(2, ISEED, LENJ, AB(KD + 2 - LENJ, J).asArray());
        AB[KD + 1][J] = sign(ONE, AB[KD + 1][J]);
      }
      AB[KD + 1][N] = SMLNUM * AB[KD + 1][N];
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J + 1, KD + 1);
        dlarnv(2, ISEED, LENJ, AB(1, J).asArray());
        AB[1][J] = sign(ONE, AB[1][J]);
      }
      AB[1][1] = SMLNUM * AB[1][1];
    }
  } else if (IMAT == 13) {
    // Type 13:  T is diagonal with small numbers on the diagonal to
    // make the growth factor underflow, but a small right hand side
    // chosen so that the solution does not overflow.

    if (UPPER) {
      var JCOUNT = 1;
      for (var J = N; J >= 1; J--) {
        for (var I = max(1, KD + 1 - (J - 1)); I <= KD; I++) {
          AB[I][J] = ZERO;
        }
        if (JCOUNT <= 2) {
          AB[KD + 1][J] = SMLNUM;
        } else {
          AB[KD + 1][J] = ONE;
        }
        JCOUNT++;
        if (JCOUNT > 4) JCOUNT = 1;
      }
    } else {
      var JCOUNT = 1;
      for (var J = 1; J <= N; J++) {
        for (var I = 2; I <= min(N - J + 1, KD + 1); I++) {
          AB[I][J] = ZERO;
        }
        if (JCOUNT <= 2) {
          AB[1][J] = SMLNUM;
        } else {
          AB[1][J] = ONE;
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
  } else if (IMAT == 14) {
    // Type 14:  Make the diagonal elements small to cause gradual
    // overflow when dividing by T(j,j).  To control the amount of
    // scaling needed, the matrix is bidiagonal.

    final TEXP = ONE / (KD + 1);
    final TSCAL = pow(SMLNUM, TEXP).toDouble();
    dlarnv(2, ISEED, N, B);
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = max(1, KD + 2 - J); I <= KD; I++) {
          AB[I][J] = ZERO;
        }
        if (J > 1 && KD > 0) AB[KD][J] = -ONE;
        AB[KD + 1][J] = TSCAL;
      }
      B[N] = ONE;
    } else {
      for (var J = 1; J <= N; J++) {
        for (var I = 3; I <= min(N - J + 1, KD + 1); I++) {
          AB[I][J] = ZERO;
        }
        if (J < N && KD > 0) AB[2][J] = -ONE;
        AB[1][J] = TSCAL;
      }
      B[1] = ONE;
    }
  } else if (IMAT == 15) {
    // Type 15:  One zero diagonal element.

    final IY = N ~/ 2 + 1;
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J, KD + 1);
        dlarnv(2, ISEED, LENJ, AB(KD + 2 - LENJ, J).asArray());
        if (J != IY) {
          AB[KD + 1][J] = sign(TWO, AB[KD + 1][J]);
        } else {
          AB[KD + 1][J] = ZERO;
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J + 1, KD + 1);
        dlarnv(2, ISEED, LENJ, AB(1, J).asArray());
        if (J != IY) {
          AB[1][J] = sign(TWO, AB[1][J]);
        } else {
          AB[1][J] = ZERO;
        }
      }
    }
    dlarnv(2, ISEED, N, B);
    dscal(N, TWO, B, 1);
  } else if (IMAT == 16) {
    // Type 16:  Make the offdiagonal elements large to cause overflow
    // when adding a column of T.  In the non-transposed case, the
    // matrix is constructed to cause overflow when adding a column in
    // every other step.

    final TSCAL = (ONE - ULP) / (UNFL / ULP);
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= KD + 1; I++) {
        AB[I][J] = ZERO;
      }
    }
    var TEXP = ONE;
    if (KD > 0) {
      if (UPPER) {
        for (var J = N; J >= 1; J -= KD) {
          for (var I = J; I >= max(1, J - KD + 1); I -= 2) {
            AB[1 + (J - I)][I] = -TSCAL / (KD + 2);
            AB[KD + 1][I] = ONE;
            B[I] = TEXP * (ONE - ULP);
            if (I > max(1, J - KD + 1)) {
              AB[2 + (J - I)][I - 1] = -(TSCAL / (KD + 2)) / (KD + 3);
              AB[KD + 1][I - 1] = ONE;
              B[I - 1] = TEXP * (KD + 1) * (KD + 1) + KD;
            }
            TEXP *= TWO;
          }
          B[max(1, J - KD + 1)] = ((KD + 2) / (KD + 3)) * TSCAL;
        }
      } else {
        for (var J = 1; J <= N; J += KD) {
          TEXP = ONE;
          final LENJ = min(KD + 1, N - J + 1);
          for (var I = J; I <= min(N, J + KD - 1); I += 2) {
            AB[LENJ - (I - J)][J] = -TSCAL / (KD + 2);
            AB[1][J] = ONE;
            B[J] = TEXP * (ONE - ULP);
            if (I < min(N, J + KD - 1)) {
              AB[LENJ - (I - J + 1)][I + 1] = -(TSCAL / (KD + 2)) / (KD + 3);
              AB[1][I + 1] = ONE;
              B[I + 1] = TEXP * (KD + 1) * (KD + 1) + KD;
            }
            TEXP *= TWO;
          }
          B[min(N, J + KD - 1)] = ((KD + 2) / (KD + 3)) * TSCAL;
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        AB[1][J] = ONE;
        B[J] = J.toDouble();
      }
    }
  } else if (IMAT == 17) {
    // Type 17:  Generate a unit triangular matrix with elements
    // between -1 and 1, and make the right hand side large so that it
    // requires scaling.

    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J - 1, KD);
        dlarnv(2, ISEED, LENJ, AB(KD + 1 - LENJ, J).asArray());
        AB[KD + 1][J] = J.toDouble();
      }
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J, KD);
        if (LENJ > 0) dlarnv(2, ISEED, LENJ, AB(2, J).asArray());
        AB[1][J] = J.toDouble();
      }
    }

    // Set the right hand side so that the largest value is BIGNUM.

    dlarnv(2, ISEED, N, B);
    final IY = idamax(N, B, 1);
    final BNORM = B[IY].abs();
    final BSCAL = BIGNUM / max(ONE, BNORM);
    dscal(N, BSCAL, B, 1);
  } else if (IMAT == 18) {
    // Type 18:  Generate a triangular matrix with elements between
    // BIGNUM/KD and BIGNUM so that at least one of the column
    // norms will exceed BIGNUM.

    final TLEFT = BIGNUM / max(ONE, KD);
    final TSCAL = BIGNUM * (KD / (KD + 1));
    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(J, KD + 1);
        dlarnv(2, ISEED, LENJ, AB(KD + 2 - LENJ, J).asArray());
        for (var I = KD + 2 - LENJ; I <= KD + 1; I++) {
          AB[I][J] = sign(TLEFT, AB[I][J]) + TSCAL * AB[I][J];
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        final LENJ = min(N - J + 1, KD + 1);
        dlarnv(2, ISEED, LENJ, AB(1, J).asArray());
        for (var I = 1; I <= LENJ; I++) {
          AB[I][J] = sign(TLEFT, AB[I][J]) + TSCAL * AB[I][J];
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
        final LENJ = min(N - 2 * J + 1, KD + 1);
        dswap(LENJ, AB(KD + 1, J).asArray(), LDAB - 1,
            AB(KD + 2 - LENJ, N - J + 1).asArray(), -1);
      }
    } else {
      for (var J = 1; J <= N ~/ 2; J++) {
        final LENJ = min(N - 2 * J + 1, KD + 1);
        dswap(LENJ, AB(1, J).asArray(), 1, AB(LENJ, N - J + 2 - LENJ).asArray(),
            -LDAB + 1);
      }
    }
  }
}
