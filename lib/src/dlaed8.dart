import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlamrg.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlaed8(
  final int ICOMPQ,
  final Box<int> K,
  final int N,
  final int QSIZ,
  final Array<double> D,
  final Matrix<double> Q,
  final int LDQ,
  final Array<int> INDXQ,
  final Box<double> RHO,
  final int CUTPNT,
  final Array<double> Z,
  final Array<double> DLAMBDA,
  final Matrix<double> Q2,
  final int LDQ2,
  final Array<double> W,
  final Array<int> PERM,
  final Box<int> GIVPTR,
  final Matrix<int> GIVCOL,
  final Matrix<double> GIVNUM,
  final Array<int> INDXP,
  final Array<int> INDX,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const MONE = -1.0, ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0;

  int I, IMAX, J, JLAM = 0, JMAX, JP, K2, N1, N1P1, N2;
  double C, EPS, S, T, TAU, TOL;

  // Test the input parameters.

  INFO.value = 0;

  if (ICOMPQ < 0 || ICOMPQ > 1) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (ICOMPQ == 1 && QSIZ < N) {
    INFO.value = -4;
  } else if (LDQ < max(1, N)) {
    INFO.value = -7;
  } else if (CUTPNT < min(1, N) || CUTPNT > N) {
    INFO.value = -10;
  } else if (LDQ2 < max(1, N)) {
    INFO.value = -14;
  }
  if (INFO.value != 0) {
    xerbla('DLAED8', -INFO.value);
    return;
  }

  // Need to initialize GIVPTR.value to O here in case of quick exit
  // to prevent an unspecified code behavior (usually sigfault)
  // when IWORK array on entry to *stedc is not zeroed
  // (or at least some IWORK entries which used in *laed7 for GIVPTR.value).

  GIVPTR.value = 0;

  // Quick return if possible

  if (N == 0) return;

  N1 = CUTPNT;
  N2 = N - N1;
  N1P1 = N1 + 1;

  if (RHO.value < ZERO) {
    dscal(N2, MONE, Z(N1P1), 1);
  }

  // Normalize z so that norm(z) = 1

  T = ONE / sqrt(TWO);
  for (J = 1; J <= N; J++) {
    INDX[J] = J;
  }
  dscal(N, T, Z, 1);
  RHO.value = (TWO * RHO.value).abs();

  // Sort the eigenvalues into increasing order

  for (I = CUTPNT + 1; I <= N; I++) {
    INDXQ[I] = INDXQ[I] + CUTPNT;
  }
  for (I = 1; I <= N; I++) {
    DLAMBDA[I] = D[INDXQ[I]];
    W[I] = Z[INDXQ[I]];
  }
  I = 1;
  J = CUTPNT + 1;
  dlamrg(N1, N2, DLAMBDA, 1, 1, INDX);
  for (I = 1; I <= N; I++) {
    D[I] = DLAMBDA[INDX[I]];
    Z[I] = W[INDX[I]];
  }

  // Calculate the allowable deflation tolerance

  IMAX = idamax(N, Z, 1);
  JMAX = idamax(N, D, 1);
  EPS = dlamch('Epsilon');
  TOL = EIGHT * EPS * (D[JMAX]).abs();

  // If the rank-1 modifier is small enough, no more needs to be done
  // except to reorganize Q so that its columns correspond with the
  // elements in D.

  if (RHO.value * (Z[IMAX]).abs() <= TOL) {
    K.value = 0;
    if (ICOMPQ == 0) {
      for (J = 1; J <= N; J++) {
        PERM[J] = INDXQ[INDX[J]];
      }
    } else {
      for (J = 1; J <= N; J++) {
        PERM[J] = INDXQ[INDX[J]];
        dcopy(QSIZ, Q(1, PERM[J]).asArray(), 1, Q2(1, J).asArray(), 1);
      }
      dlacpy('A', QSIZ, N, Q2, LDQ2, Q, LDQ);
    }
    return;
  }

  // If there are multiple eigenvalues then the problem deflates.  Here
  // the number of equal eigenvalues are found.  As each equal
  // eigenvalue is found, an elementary reflector is computed to rotate
  // the corresponding eigensubspace so that the corresponding
  // components of Z are zero in this new basis.

  K.value = 0;
  K2 = N + 1;
  var isLastItem = false;
  for (J = 1; J <= N; J++) {
    if (RHO.value * (Z[J]).abs() <= TOL) {
      // Deflate due to small z component.

      K2 = K2 - 1;
      INDXP[K2] = J;
      if (J == N) {
        isLastItem = true;
        break;
      }
    } else {
      JLAM = J;
      break;
    }
  }

  while (!isLastItem) {
    J = J + 1;
    if (J > N) break;
    if (RHO.value * (Z[J]).abs() <= TOL) {
      // Deflate due to small z component.

      K2 = K2 - 1;
      INDXP[K2] = J;
    } else {
      // Check if eigenvalues are close enough to allow deflation.

      S = Z[JLAM];
      C = Z[J];

      // Find sqrt(a**2+b**2) without overflow or
      // destructive underflow.

      TAU = dlapy2(C, S);
      T = D[J] - D[JLAM];
      C = C / TAU;
      S = -S / TAU;
      if ((T * C * S).abs() <= TOL) {
        // Deflation is possible.

        Z[J] = TAU;
        Z[JLAM] = ZERO;

        // Record the appropriate Givens rotation

        GIVPTR.value = GIVPTR.value + 1;
        GIVCOL[1][GIVPTR.value] = INDXQ[INDX[JLAM]];
        GIVCOL[2][GIVPTR.value] = INDXQ[INDX[J]];
        GIVNUM[1][GIVPTR.value] = C;
        GIVNUM[2][GIVPTR.value] = S;
        if (ICOMPQ == 1) {
          drot(QSIZ, Q(1, INDXQ[INDX[JLAM]]).asArray(), 1,
              Q(1, INDXQ[INDX[J]]).asArray(), 1, C, S);
        }
        T = D[JLAM] * C * C + D[J] * S * S;
        D[J] = D[JLAM] * S * S + D[J] * C * C;
        D[JLAM] = T;
        K2 = K2 - 1;
        I = 1;

        while (K2 + I <= N && D[JLAM] < D[INDXP[K2 + I]]) {
          INDXP[K2 + I - 1] = INDXP[K2 + I];
          INDXP[K2 + I] = JLAM;
          I = I + 1;
        }
        INDXP[K2 + I - 1] = JLAM;

        JLAM = J;
      } else {
        K.value = K.value + 1;
        W[K.value] = Z[JLAM];
        DLAMBDA[K.value] = D[JLAM];
        INDXP[K.value] = JLAM;
        JLAM = J;
      }
    }
  }

  // Record the last eigenvalue.
  if (!isLastItem) {
    K.value = K.value + 1;
    W[K.value] = Z[JLAM];
    DLAMBDA[K.value] = D[JLAM];
    INDXP[K.value] = JLAM;
  }

  // Sort the eigenvalues and corresponding eigenvectors into DLAMBDA
  // and Q2 respectively.  The eigenvalues/vectors which were not
  // deflated go into the first K.value slots of DLAMBDA and Q2 respectively,
  // while those which were deflated go into the last N - K.value slots.

  if (ICOMPQ == 0) {
    for (J = 1; J <= N; J++) {
      JP = INDXP[J];
      DLAMBDA[J] = D[JP];
      PERM[J] = INDXQ[INDX[JP]];
    }
  } else {
    for (J = 1; J <= N; J++) {
      JP = INDXP[J];
      DLAMBDA[J] = D[JP];
      PERM[J] = INDXQ[INDX[JP]];
      dcopy(QSIZ, Q(1, PERM[J]).asArray(), 1, Q2(1, J).asArray(), 1);
    }
  }

  // The deflated eigenvalues and their corresponding vectors go back
  // into the last N - K.value slots of D and Q respectively.

  if (K.value < N) {
    if (ICOMPQ == 0) {
      dcopy(N - K.value, DLAMBDA(K.value + 1), 1, D(K.value + 1), 1);
    } else {
      dcopy(N - K.value, DLAMBDA(K.value + 1), 1, D(K.value + 1), 1);
      dlacpy('A', QSIZ, N - K.value, Q2(1, K.value + 1), LDQ2,
          Q(1, K.value + 1), LDQ);
    }
  }
}
