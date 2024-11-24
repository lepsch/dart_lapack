// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/drot.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlanv2.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/dlarfx.dart';
import 'package:dart_lapack/src/dlartg.dart';
import 'package:dart_lapack/src/dlasy2.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlaexc(
  final bool WANTQ,
  final int N,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> Q_,
  final int LDQ,
  final int J1,
  final int N1,
  final int N2,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final T = T_.having(ld: LDT);
  final Q = Q_.having(ld: LDQ);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const TEN = 1.0e+1;
  const LDD = 4, LDX = 2;
  int J2, J3, J4, K, ND;
  double DNORM, EPS, SMLNUM, T11, T22, T33, THRESH;
  final D = Matrix<double>(LDD, 4), X = Matrix<double>(LDX, 2);
  final U = Array<double>(3), U1 = Array<double>(3), U2 = Array<double>(3);
  final CS = Box(0.0),
      SN = Box(0.0),
      SCALE = Box(0.0),
      XNORM = Box(0.0),
      TEMP = Box(0.0),
      TAU = Box(0.0),
      TAU1 = Box(0.0),
      TAU2 = Box(0.0),
      WI1 = Box(0.0),
      WI2 = Box(0.0),
      WR1 = Box(0.0),
      WR2 = Box(0.0);
  final IERR = Box(0);

  INFO.value = 0;

  // Quick return if possible

  if (N == 0 || N1 == 0 || N2 == 0) return;
  if (J1 + N1 > N) return;

  J2 = J1 + 1;
  J3 = J1 + 2;
  J4 = J1 + 3;

  if (N1 == 1 && N2 == 1) {
    // Swap two 1-by-1 blocks.

    T11 = T[J1][J1];
    T22 = T[J2][J2];

    // Determine the transformation to perform the interchange.

    dlartg(T[J1][J2], T22 - T11, CS, SN, TEMP);

    // Apply transformation to the matrix T.

    if (J3 <= N) {
      drot(N - J1 - 1, T(J1, J3).asArray(), LDT, T(J2, J3).asArray(), LDT,
          CS.value, SN.value);
    }
    drot(J1 - 1, T(1, J1).asArray(), 1, T(1, J2).asArray(), 1, CS.value,
        SN.value);

    T[J1][J1] = T22;
    T[J2][J2] = T11;

    if (WANTQ) {
      // Accumulate transformation in the matrix Q.

      drot(N, Q(1, J1).asArray(), 1, Q(1, J2).asArray(), 1, CS.value, SN.value);
    }
  } else {
    // Swapping involves at least one 2-by-2 block.

    // Copy the diagonal block of order N1+N2 to the local array D
    // and compute its norm.

    ND = N1 + N2;
    dlacpy('Full', ND, ND, T(J1, J1), LDT, D, LDD);
    DNORM = dlange('Max', ND, ND, D, LDD, WORK);

    // Compute machine-dependent threshold for test for accepting
    // swap.

    EPS = dlamch('P');
    SMLNUM = dlamch('S') / EPS;
    THRESH = max(TEN * EPS * DNORM, SMLNUM);

    // Solve T11*X - X*T22 = scale*T12 for X.

    dlasy2(false, false, -1, N1, N2, D, LDD, D(N1 + 1, N1 + 1), LDD,
        D(1, N1 + 1), LDD, SCALE, X, LDX, XNORM, IERR);

    // Swap the adjacent diagonal blocks.

    K = N1 + N1 + N2 - 3;
    //  GO TO ( 10, 20, 30 )K;
    switch (K) {
      case 1:

        // N1 = 1, N2 = 2: generate elementary reflector H so that:

        // ( scale, X11, X12 ) H = ( 0, 0, * )

        U[1] = SCALE.value;
        U[2] = X[1][1];
        U[3] = X[1][2];
        dlarfg(3, U.box(3), U, 1, TAU);
        U[3] = ONE;
        T11 = T[J1][J1];

        // Perform swap provisionally on diagonal block in D.

        dlarfx('L', 3, 3, U, TAU.value, D, LDD, WORK);
        dlarfx('R', 3, 3, U, TAU.value, D, LDD, WORK);

        // Test whether to reject swap.

        if (max(D[3][1].abs(), max(D[3][2].abs(), (D[3][3] - T11).abs())) >
            THRESH) {
          // Exit with INFO = 1 if swap was rejected.
          INFO.value = 1;
          return;
        }

        // Accept swap: apply transformation to the entire matrix T.

        dlarfx('L', 3, N - J1 + 1, U, TAU.value, T(J1, J1), LDT, WORK);
        dlarfx('R', J2, 3, U, TAU.value, T(1, J1), LDT, WORK);

        T[J3][J1] = ZERO;
        T[J3][J2] = ZERO;
        T[J3][J3] = T11;

        if (WANTQ) {
          // Accumulate transformation in the matrix Q.

          dlarfx('R', N, 3, U, TAU.value, Q(1, J1), LDQ, WORK);
        }
        break;

      case 2:

        // N1 = 2, N2 = 1: generate elementary reflector H so that:

        // H (  -X11 ) = ( * )
        // (  -X21 ) = ( 0 )
        // ( scale ) = ( 0 )

        U[1] = -X[1][1];
        U[2] = -X[2][1];
        U[3] = SCALE.value;
        dlarfg(3, U.box(1), U(2), 1, TAU);
        U[1] = ONE;
        T33 = T[J3][J3];

        // Perform swap provisionally on diagonal block in D.

        dlarfx('L', 3, 3, U, TAU.value, D, LDD, WORK);
        dlarfx('R', 3, 3, U, TAU.value, D, LDD, WORK);

        // Test whether to reject swap.

        if (max(D[2][1].abs(), max(D[3][1].abs(), (D[1][1] - T33).abs())) >
            THRESH) {
          // Exit with INFO = 1 if swap was rejected.
          INFO.value = 1;
          return;
        }

        // Accept swap: apply transformation to the entire matrix T.

        dlarfx('R', J3, 3, U, TAU.value, T(1, J1), LDT, WORK);
        dlarfx('L', 3, N - J1, U, TAU.value, T(J1, J2), LDT, WORK);

        T[J1][J1] = T33;
        T[J2][J1] = ZERO;
        T[J3][J1] = ZERO;

        if (WANTQ) {
          // Accumulate transformation in the matrix Q.

          dlarfx('R', N, 3, U, TAU.value, Q(1, J1), LDQ, WORK);
        }
        break;

      case 3:

        // N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
        // that:
        //
        // H(2) H(1) (  -X11  -X12 ) = (  *  * )
        // (  -X21  -X22 )   (  0  * )
        // ( scale    0  )   (  0  0 )
        // (    0  scale )   (  0  0 )

        U1[1] = -X[1][1];
        U1[2] = -X[2][1];
        U1[3] = SCALE.value;
        dlarfg(3, U1.box(1), U1(2), 1, TAU1);
        U1[1] = ONE;

        TEMP.value = -TAU1.value * (X[1][2] + U1[2] * X[2][2]);
        U2[1] = -TEMP.value * U1[2] - X[2][2];
        U2[2] = -TEMP.value * U1[3];
        U2[3] = SCALE.value;
        dlarfg(3, U2.box(1), U2(2), 1, TAU2);
        U2[1] = ONE;

        // Perform swap provisionally on diagonal block in D.

        dlarfx('L', 3, 4, U1, TAU1.value, D, LDD, WORK);
        dlarfx('R', 4, 3, U1, TAU1.value, D, LDD, WORK);
        dlarfx('L', 3, 4, U2, TAU2.value, D(2, 1), LDD, WORK);
        dlarfx('R', 4, 3, U2, TAU2.value, D(1, 2), LDD, WORK);

        // Test whether to reject swap.

        if (max(
              max(D[3][1].abs(), D[3][2].abs()),
              max(D[4][1].abs(), D[4][2].abs()),
            ) >
            THRESH) {
          // Exit with INFO = 1 if swap was rejected.
          INFO.value = 1;
          return;
        }

        // Accept swap: apply transformation to the entire matrix T.

        dlarfx('L', 3, N - J1 + 1, U1, TAU1.value, T(J1, J1), LDT, WORK);
        dlarfx('R', J4, 3, U1, TAU1.value, T(1, J1), LDT, WORK);
        dlarfx('L', 3, N - J1 + 1, U2, TAU2.value, T(J2, J1), LDT, WORK);
        dlarfx('R', J4, 3, U2, TAU2.value, T(1, J2), LDT, WORK);

        T[J3][J1] = ZERO;
        T[J3][J2] = ZERO;
        T[J4][J1] = ZERO;
        T[J4][J2] = ZERO;

        if (WANTQ) {
          // Accumulate transformation in the matrix Q.

          dlarfx('R', N, 3, U1, TAU1.value, Q(1, J1), LDQ, WORK);
          dlarfx('R', N, 3, U2, TAU2.value, Q(1, J2), LDQ, WORK);
        }
        break;
    }

    if (N2 == 2) {
      // Standardize new 2-by-2 block T11

      dlanv2(T.box(J1, J1), T.box(J1, J2), T.box(J2, J1), T.box(J2, J2), WR1,
          WI1, WR2, WI2, CS, SN);
      drot(N - J1 - 1, T(J1, J1 + 2).asArray(), LDT, T(J2, J1 + 2).asArray(),
          LDT, CS.value, SN.value);
      drot(J1 - 1, T(1, J1).asArray(), 1, T(1, J2).asArray(), 1, CS.value,
          SN.value);
      if (WANTQ) {
        drot(N, Q(1, J1).asArray(), 1, Q(1, J2).asArray(), 1, CS.value,
            SN.value);
      }
    }

    if (N1 == 2) {
      // Standardize new 2-by-2 block T22

      J3 = J1 + N2;
      J4 = J3 + 1;
      dlanv2(T.box(J3, J3), T.box(J3, J4), T.box(J4, J3), T.box(J4, J4), WR1,
          WI1, WR2, WI2, CS, SN);
      if (J3 + 2 <= N) {
        drot(N - J3 - 1, T(J3, J3 + 2).asArray(), LDT, T(J4, J3 + 2).asArray(),
            LDT, CS.value, SN.value);
      }
      drot(J3 - 1, T(1, J3).asArray(), 1, T(1, J4).asArray(), 1, CS.value,
          SN.value);
      if (WANTQ) {
        drot(N, Q(1, J3).asArray(), 1, Q(1, J4).asArray(), 1, CS.value,
            SN.value);
      }
    }
  }
}
