// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/drot.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlamrg.dart';
import 'package:dart_lapack/src/dlapy2.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlaed2(
  final Box<int> K,
  final int N,
  final int N1,
  final Array<double> D_,
  final Matrix<double> Q_,
  final int LDQ,
  final Array<int> INDXQ_,
  final Box<double> RHO,
  final Array<double> Z_,
  final Array<double> DLAMBDA_,
  final Array<double> W_,
  final Array<double> Q2_,
  final Array<int> INDX_,
  final Array<int> INDXC_,
  final Array<int> INDXP_,
  final Array<int> COLTYP_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final Q = Q_.having(ld: LDQ);
  final INDXQ = INDXQ_.having();
  final Z = Z_.having();
  final DLAMBDA = DLAMBDA_.having();
  final W = W_.having();
  final Q2 = Q2_.having();
  final INDX = INDX_.having();
  final INDXC = INDXC_.having();
  final INDXP = INDXP_.having();
  final COLTYP = COLTYP_.having();
  const MONE = -1.0, ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0;
  final CTOT = Array<int>(4), PSM = Array<int>(4);
  int CT, I, IMAX, IQ1, IQ2, J, JMAX, JS, K2, N1P1, N2, NJ, PJ = 0;
  double C, EPS, S, T, TAU, TOL;

  // Test the input parameters.

  INFO.value = 0;

  if (N < 0) {
    INFO.value = -2;
  } else if (LDQ < max(1, N)) {
    INFO.value = -6;
  } else if (min(1, N ~/ 2) > N1 || (N ~/ 2) < N1) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('DLAED2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  N2 = N - N1;
  N1P1 = N1 + 1;

  if (RHO.value < ZERO) {
    dscal(N2, MONE, Z(N1P1), 1);
  }

  // Normalize z so that norm(z) = 1.  Since z is the concatenation of
  // two normalized vectors, norm2(z) = sqrt(2).

  T = ONE / sqrt(TWO);
  dscal(N, T, Z, 1);

  // RHO = ABS( norm(z)**2 * RHO )

  RHO.value = (TWO * RHO.value).abs();

  // Sort the eigenvalues into increasing order

  for (I = N1P1; I <= N; I++) {
    INDXQ[I] += N1;
  }

  // re-integrate the deflated parts from the last pass

  for (I = 1; I <= N; I++) {
    DLAMBDA[I] = D[INDXQ[I]];
  }
  dlamrg(N1, N2, DLAMBDA, 1, 1, INDXC);
  for (I = 1; I <= N; I++) {
    INDX[I] = INDXQ[INDXC[I]];
  }

  // Calculate the allowable deflation tolerance

  IMAX = idamax(N, Z, 1);
  JMAX = idamax(N, D, 1);
  EPS = dlamch('Epsilon');
  TOL = EIGHT * EPS * max(D[JMAX].abs(), Z[IMAX].abs());

  // If the rank-1 modifier is small enough, no more needs to be done
  // except to reorganize Q so that its columns correspond with the
  // elements in D.

  if (RHO.value * Z[IMAX].abs() <= TOL) {
    K.value = 0;
    IQ2 = 1;
    for (J = 1; J <= N; J++) {
      I = INDX[J];
      dcopy(N, Q(1, I).asArray(), 1, Q2(IQ2), 1);
      DLAMBDA[J] = D[I];
      IQ2 += N;
    }
    dlacpy('A', N, N, Q2.asMatrix(N), N, Q, LDQ);
    dcopy(N, DLAMBDA, 1, D, 1);
    return;
  }

  // If there are multiple eigenvalues then the problem deflates.  Here
  // the number of equal eigenvalues are found.  As each equal
  // eigenvalue is found, an elementary reflector is computed to rotate
  // the corresponding eigensubspace so that the corresponding
  // components of Z are zero in this new basis.

  for (I = 1; I <= N1; I++) {
    COLTYP[I] = 1;
  }
  for (I = N1P1; I <= N; I++) {
    COLTYP[I] = 3;
  }

  K.value = 0;
  K2 = N + 1;
  var isLastItem = false;
  for (J = 1; J <= N; J++) {
    NJ = INDX[J];
    if (RHO.value * Z[NJ].abs() <= TOL) {
      // Deflate due to small z component.

      K2--;
      COLTYP[NJ] = 4;
      INDXP[K2] = NJ;
      if (J == N) {
        isLastItem = true;
        break;
      }
    } else {
      PJ = NJ;
      break;
    }
  }

  while (!isLastItem) {
    J++;
    NJ = INDX[J];
    if (J > N) break;
    if (RHO.value * Z[NJ].abs() <= TOL) {
      // Deflate due to small z component.

      K2--;
      COLTYP[NJ] = 4;
      INDXP[K2] = NJ;
    } else {
      // Check if eigenvalues are close enough to allow deflation.

      S = Z[PJ];
      C = Z[NJ];

      // Find sqrt(a**2+b**2) without overflow or
      // destructive underflow.

      TAU = dlapy2(C, S);
      T = D[NJ] - D[PJ];
      C /= TAU;
      S = -S / TAU;
      if ((T * C * S).abs() <= TOL) {
        // Deflation is possible.

        Z[NJ] = TAU;
        Z[PJ] = ZERO;
        if (COLTYP[NJ] != COLTYP[PJ]) COLTYP[NJ] = 2;
        COLTYP[PJ] = 4;
        drot(N, Q(1, PJ).asArray(), 1, Q(1, NJ).asArray(), 1, C, S);
        T = D[PJ] * pow(C, 2) + D[NJ] * pow(S, 2);
        D[NJ] = D[PJ] * pow(S, 2) + D[NJ] * pow(C, 2);
        D[PJ] = T;
        K2--;
        I = 1;

        while (K2 + I <= N && D[PJ] >= D[INDXP[K2 + I]]) {
          INDXP[K2 + I - 1] = INDXP[K2 + I];
          INDXP[K2 + I] = PJ;
          I++;
        }
        INDXP[K2 + I - 1] = PJ;

        PJ = NJ;
      } else {
        K.value++;
        DLAMBDA[K.value] = D[PJ];
        W[K.value] = Z[PJ];
        INDXP[K.value] = PJ;
        PJ = NJ;
      }
    }
  }

  // Record the last eigenvalue.

  K.value++;
  DLAMBDA[K.value] = D[PJ];
  W[K.value] = Z[PJ];
  INDXP[K.value] = PJ;

  // Count up the total number of the various types of columns, then
  // form a permutation which positions the four column types into
  // four uniform groups (although one or more of these groups may be
  // empty).

  for (J = 1; J <= 4; J++) {
    CTOT[J] = 0;
  }
  for (J = 1; J <= N; J++) {
    CT = COLTYP[J];
    CTOT[CT]++;
  }

  // PSM[*] = Position in SubMatrix (of types 1 through 4)

  PSM[1] = 1;
  PSM[2] = 1 + CTOT[1];
  PSM[3] = PSM[2] + CTOT[2];
  PSM[4] = PSM[3] + CTOT[3];
  K.value = N - CTOT[4];

  // Fill out the INDXC array so that the permutation which it induces
  // will place all type-1 columns first, all type-2 columns next,
  // then all type-3's, and finally all type-4's.

  for (J = 1; J <= N; J++) {
    JS = INDXP[J];
    CT = COLTYP[JS];
    INDX[PSM[CT]] = JS;
    INDXC[PSM[CT]] = J;
    PSM[CT]++;
  }

  // Sort the eigenvalues and corresponding eigenvectors into DLAMBDA
  // and Q2 respectively.  The eigenvalues/vectors which were not
  // deflated go into the first K slots of DLAMBDA and Q2 respectively,
  // while those which were deflated go into the last N - K slots.

  I = 1;
  IQ1 = 1;
  IQ2 = 1 + (CTOT[1] + CTOT[2]) * N1;
  for (J = 1; J <= CTOT[1]; J++) {
    JS = INDX[I];
    dcopy(N1, Q(1, JS).asArray(), 1, Q2(IQ1), 1);
    Z[I] = D[JS];
    I++;
    IQ1 += N1;
  }

  for (J = 1; J <= CTOT[2]; J++) {
    JS = INDX[I];
    dcopy(N1, Q(1, JS).asArray(), 1, Q2(IQ1), 1);
    dcopy(N2, Q(N1 + 1, JS).asArray(), 1, Q2(IQ2), 1);
    Z[I] = D[JS];
    I++;
    IQ1 += N1;
    IQ2 += N2;
  }

  for (J = 1; J <= CTOT[3]; J++) {
    JS = INDX[I];
    dcopy(N2, Q(N1 + 1, JS).asArray(), 1, Q2(IQ2), 1);
    Z[I] = D[JS];
    I++;
    IQ2 += N2;
  }

  IQ1 = IQ2;
  for (J = 1; J <= CTOT[4]; J++) {
    JS = INDX[I];
    dcopy(N, Q(1, JS).asArray(), 1, Q2(IQ2), 1);
    IQ2 += N;
    Z[I] = D[JS];
    I++;
  }

  // The deflated eigenvalues and their corresponding vectors go back
  // into the last N - K slots of D and Q respectively.

  if (K.value < N) {
    dlacpy('A', N, CTOT[4], Q2(IQ1).asMatrix(N), N, Q(1, K.value + 1), LDQ);
    dcopy(N - K.value, Z(K.value + 1), 1, D(K.value + 1), 1);
  }

  // Copy CTOT into COLTYP for referencing in DLAED3.

  for (J = 1; J <= 4; J++) {
    COLTYP[J] = CTOT[J];
  }
}
