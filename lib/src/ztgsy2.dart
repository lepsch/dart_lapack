// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgesc2.dart';
import 'package:dart_lapack/src/zgetc2.dart';
import 'package:dart_lapack/src/zlatdf.dart';

void ztgsy2(
  final String TRANS,
  final int IJOB,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> C_,
  final int LDC,
  final Matrix<Complex> D_,
  final int LDD,
  final Matrix<Complex> E_,
  final int LDE,
  final Matrix<Complex> F_,
  final int LDF,
  final Box<double> SCALE,
  final Box<double> RDSUM,
  final Box<double> RDSCAL,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final D = D_.having(ld: LDD);
  final E = E_.having(ld: LDE);
  final F = F_.having(ld: LDF);
  const ZERO = 0.0, ONE = 1.0, LDZ = 2;
  bool NOTRAN;
  int I, J, K;
  Complex ALPHA;
  final IPIV = Array<int>(LDZ), JPIV = Array<int>(LDZ);
  final RHS = Array<Complex>(LDZ), Z = Matrix<Complex>(LDZ, LDZ);
  final SCALOC = Box(0.0);
  final IERR = Box(0);

  // Decode and test input parameters

  INFO.value = 0;
  IERR.value = 0;
  NOTRAN = lsame(TRANS, 'N');
  if (!NOTRAN && !lsame(TRANS, 'C')) {
    INFO.value = -1;
  } else if (NOTRAN) {
    if ((IJOB < 0) || (IJOB > 2)) {
      INFO.value = -2;
    }
  }
  if (INFO.value == 0) {
    if (M <= 0) {
      INFO.value = -3;
    } else if (N <= 0) {
      INFO.value = -4;
    } else if (LDA < max(1, M)) {
      INFO.value = -6;
    } else if (LDB < max(1, N)) {
      INFO.value = -8;
    } else if (LDC < max(1, M)) {
      INFO.value = -10;
    } else if (LDD < max(1, M)) {
      INFO.value = -12;
    } else if (LDE < max(1, N)) {
      INFO.value = -14;
    } else if (LDF < max(1, M)) {
      INFO.value = -16;
    }
  }
  if (INFO.value != 0) {
    xerbla('ZTGSY2', -INFO.value);
    return;
  }

  if (NOTRAN) {
    // Solve (I, J) - system
    //    A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
    //    D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
    // for I = M, M - 1, ..., 1; J = 1, 2, ..., N

    SCALE.value = ONE;
    SCALOC.value = ONE;
    for (J = 1; J <= N; J++) {
      for (I = M; I >= 1; I--) {
        // Build 2 by 2 system

        Z[1][1] = A[I][I];
        Z[2][1] = D[I][I];
        Z[1][2] = -B[J][J];
        Z[2][2] = -E[J][J];

        // Set up right hand side(s)

        RHS[1] = C[I][J];
        RHS[2] = F[I][J];

        // Solve Z * x = RHS

        zgetc2(LDZ, Z, LDZ, IPIV, JPIV, IERR);
        if (IERR.value > 0) INFO.value = IERR.value;
        if (IJOB == 0) {
          zgesc2(LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
          if (SCALOC.value != ONE) {
            for (K = 1; K <= N; K++) {
              zscal(M, Complex(SCALOC.value, ZERO), C(1, K).asArray(), 1);
              zscal(M, Complex(SCALOC.value, ZERO), F(1, K).asArray(), 1);
            }
            SCALE.value *= SCALOC.value;
          }
        } else {
          zlatdf(IJOB, LDZ, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV);
        }

        // Unpack solution vector(s)

        C[I][J] = RHS[1];
        F[I][J] = RHS[2];

        // Substitute R(I, J) and L(I, J) into remaining equation.

        if (I > 1) {
          ALPHA = -RHS[1];
          zaxpy(I - 1, ALPHA, A(1, I).asArray(), 1, C(1, J).asArray(), 1);
          zaxpy(I - 1, ALPHA, D(1, I).asArray(), 1, F(1, J).asArray(), 1);
        }
        if (J < N) {
          zaxpy(N - J, RHS[2], B(J, J + 1).asArray(), LDB,
              C(I, J + 1).asArray(), LDC);
          zaxpy(N - J, RHS[2], E(J, J + 1).asArray(), LDE,
              F(I, J + 1).asArray(), LDF);
        }
      }
    }
  } else {
    // Solve transposed (I, J) - system:
    //    A(I, I)**H * R(I, J) + D(I, I)**H * L(J, J) = C(I, J)
    //    R(I, I) * B(J, J) + L(I, J) * E(J, J)   = -F(I, J)
    // for I = 1, 2, ..., M, J = N, N - 1, ..., 1

    SCALE.value = ONE;
    SCALOC.value = ONE;
    for (I = 1; I <= M; I++) {
      for (J = N; J >= 1; J--) {
        // Build 2 by 2 system Z**H

        Z[1][1] = A[I][I].conjugate();
        Z[2][1] = -B[J][J].conjugate();
        Z[1][2] = D[I][I].conjugate();
        Z[2][2] = -E[J][J].conjugate();

        // Set up right hand side(s)

        RHS[1] = C[I][J];
        RHS[2] = F[I][J];

        // Solve Z**H * x = RHS

        zgetc2(LDZ, Z, LDZ, IPIV, JPIV, IERR);
        if (IERR.value > 0) INFO.value = IERR.value;
        zgesc2(LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
        if (SCALOC.value != ONE) {
          for (K = 1; K <= N; K++) {
            zscal(M, Complex(SCALOC.value, ZERO), C(1, K).asArray(), 1);
            zscal(M, Complex(SCALOC.value, ZERO), F(1, K).asArray(), 1);
          }
          SCALE.value *= SCALOC.value;
        }

        // Unpack solution vector(s)

        C[I][J] = RHS[1];
        F[I][J] = RHS[2];

        // Substitute R(I, J) and L(I, J) into remaining equation.

        for (K = 1; K <= J - 1; K++) {
          F[I][K] +=
              RHS[1] * B[K][J].conjugate() + RHS[2] * E[K][J].conjugate();
        }
        for (K = I + 1; K <= M; K++) {
          C[K][J] -=
              A[I][K].conjugate() * RHS[1] + D[I][K].conjugate() * RHS[2];
        }
      }
    }
  }
}
