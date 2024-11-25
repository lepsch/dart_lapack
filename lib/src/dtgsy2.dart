// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/blas/dger.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgesc2.dart';
import 'package:dart_lapack/src/dgetc2.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dlatdf.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtgsy2(
  final String TRANS,
  final int IJOB,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> C_,
  final int LDC,
  final Matrix<double> D_,
  final int LDD,
  final Matrix<double> E_,
  final int LDE,
  final Matrix<double> F_,
  final int LDF,
  final Box<double> SCALE,
  final Box<double> RDSUM,
  final Box<double> RDSCAL,
  final Array<int> IWORK_,
  final Box<int> PQ,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final D = D_.having(ld: LDD);
  final E = E_.having(ld: LDE);
  final F = F_.having(ld: LDF);
  final IWORK = IWORK_.having();
  const LDZ = 8;
  const ZERO = 0.0, ONE = 1.0;
  bool NOTRAN;
  int I, IE, II, IS, ISP1, J, JE, JJ, JS, JSP1, K, MB, NB, P, Q, ZDIM;
  double ALPHA;
  final IPIV = Array<int>(LDZ), JPIV = Array<int>(LDZ);
  final RHS = Array<double>(LDZ), Z = Matrix<double>(LDZ, LDZ);
  final IERR = Box(0);
  final SCALOC = Box(0.0);

  // Decode and test input parameters

  INFO.value = 0;
  IERR.value = 0;
  NOTRAN = lsame(TRANS, 'N');
  if (!NOTRAN && !lsame(TRANS, 'T')) {
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
    xerbla('DTGSY2', -INFO.value);
    return;
  }

  // Determine block structure of A

  PQ.value = 0;
  P = 0;
  I = 1;
  while (true) {
    if (I > M) break;
    P++;
    IWORK[P] = I;
    if (I == M) break;
    if (A[I + 1][I] != ZERO) {
      I += 2;
    } else {
      I++;
    }
  }
  IWORK[P + 1] = M + 1;

  // Determine block structure of B

  Q = P + 1;
  J = 1;
  while (true) {
    if (J > N) break;
    Q++;
    IWORK[Q] = J;
    if (J == N) break;
    if (B[J + 1][J] != ZERO) {
      J += 2;
    } else {
      J++;
    }
  }
  IWORK[Q + 1] = N + 1;
  PQ.value = P * (Q - P - 1);

  if (NOTRAN) {
    // Solve (I, J) - subsystem
    // A[I][ I] * R(I, J) - L(I, J) * B[J][ J] = C[I][ J]
    // D[I][ I] * R(I, J) - L(I, J) * E[J][ J] = F[I][ J]
    // for I = P, P - 1, ..., 1; J = 1, 2, ..., Q

    SCALE.value = ONE;
    SCALOC.value = ONE;
    for (J = P + 2; J <= Q; J++) {
      JS = IWORK[J];
      JSP1 = JS + 1;
      JE = IWORK[J + 1] - 1;
      NB = JE - JS + 1;
      for (I = P; I >= 1; I--) {
        IS = IWORK[I];
        ISP1 = IS + 1;
        IE = IWORK[I + 1] - 1;
        MB = IE - IS + 1;
        ZDIM = MB * NB * 2;

        if ((MB == 1) && (NB == 1)) {
          // Build a 2-by-2 system Z * x = RHS

          Z[1][1] = A[IS][IS];
          Z[2][1] = D[IS][IS];
          Z[1][2] = -B[JS][JS];
          Z[2][2] = -E[JS][JS];

          // Set up right hand side(s)

          RHS[1] = C[IS][JS];
          RHS[2] = F[IS][JS];

          // Solve Z * x = RHS

          dgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR);
          if (IERR.value > 0) INFO.value = IERR.value;

          if (IJOB == 0) {
            dgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
            if (SCALOC.value != ONE) {
              for (K = 1; K <= N; K++) {
                dscal(M, SCALOC.value, C(1, K).asArray(), 1);
                dscal(M, SCALOC.value, F(1, K).asArray(), 1);
              }
              SCALE.value *= SCALOC.value;
            }
          } else {
            dlatdf(IJOB, ZDIM, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV);
          }

          // Unpack solution vector(s)

          C[IS][JS] = RHS[1];
          F[IS][JS] = RHS[2];

          // Substitute R(I, J) and L(I, J) into remaining
          // equation.

          if (I > 1) {
            ALPHA = -RHS[1];
            daxpy(IS - 1, ALPHA, A(1, IS).asArray(), 1, C(1, JS).asArray(), 1);
            daxpy(IS - 1, ALPHA, D(1, IS).asArray(), 1, F(1, JS).asArray(), 1);
          }
          if (J < Q) {
            daxpy(N - JE, RHS[2], B(JS, JE + 1).asArray(), LDB,
                C(IS, JE + 1).asArray(), LDC);
            daxpy(N - JE, RHS[2], E(JS, JE + 1).asArray(), LDE,
                F(IS, JE + 1).asArray(), LDF);
          }
        } else if ((MB == 1) && (NB == 2)) {
          // Build a 4-by-4 system Z * x = RHS

          Z[1][1] = A[IS][IS];
          Z[2][1] = ZERO;
          Z[3][1] = D[IS][IS];
          Z[4][1] = ZERO;

          Z[1][2] = ZERO;
          Z[2][2] = A[IS][IS];
          Z[3][2] = ZERO;
          Z[4][2] = D[IS][IS];

          Z[1][3] = -B[JS][JS];
          Z[2][3] = -B[JS][JSP1];
          Z[3][3] = -E[JS][JS];
          Z[4][3] = -E[JS][JSP1];

          Z[1][4] = -B[JSP1][JS];
          Z[2][4] = -B[JSP1][JSP1];
          Z[3][4] = ZERO;
          Z[4][4] = -E[JSP1][JSP1];

          // Set up right hand side(s)

          RHS[1] = C[IS][JS];
          RHS[2] = C[IS][JSP1];
          RHS[3] = F[IS][JS];
          RHS[4] = F[IS][JSP1];

          // Solve Z * x = RHS

          dgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR);
          if (IERR.value > 0) INFO.value = IERR.value;

          if (IJOB == 0) {
            dgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
            if (SCALOC.value != ONE) {
              for (K = 1; K <= N; K++) {
                dscal(M, SCALOC.value, C(1, K).asArray(), 1);
                dscal(M, SCALOC.value, F(1, K).asArray(), 1);
              }
              SCALE.value *= SCALOC.value;
            }
          } else {
            dlatdf(IJOB, ZDIM, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV);
          }

          // Unpack solution vector(s)

          C[IS][JS] = RHS[1];
          C[IS][JSP1] = RHS[2];
          F[IS][JS] = RHS[3];
          F[IS][JSP1] = RHS[4];

          // Substitute R(I, J) and L(I, J) into remaining
          // equation.

          if (I > 1) {
            dger(IS - 1, NB, -ONE, A(1, IS).asArray(), 1, RHS(1), 1, C(1, JS),
                LDC);
            dger(IS - 1, NB, -ONE, D(1, IS).asArray(), 1, RHS(1), 1, F(1, JS),
                LDF);
          }
          if (J < Q) {
            daxpy(N - JE, RHS[3], B(JS, JE + 1).asArray(), LDB,
                C(IS, JE + 1).asArray(), LDC);
            daxpy(N - JE, RHS[3], E(JS, JE + 1).asArray(), LDE,
                F(IS, JE + 1).asArray(), LDF);
            daxpy(N - JE, RHS[4], B(JSP1, JE + 1).asArray(), LDB,
                C(IS, JE + 1).asArray(), LDC);
            daxpy(N - JE, RHS[4], E(JSP1, JE + 1).asArray(), LDE,
                F(IS, JE + 1).asArray(), LDF);
          }
        } else if ((MB == 2) && (NB == 1)) {
          // Build a 4-by-4 system Z * x = RHS

          Z[1][1] = A[IS][IS];
          Z[2][1] = A[ISP1][IS];
          Z[3][1] = D[IS][IS];
          Z[4][1] = ZERO;

          Z[1][2] = A[IS][ISP1];
          Z[2][2] = A[ISP1][ISP1];
          Z[3][2] = D[IS][ISP1];
          Z[4][2] = D[ISP1][ISP1];

          Z[1][3] = -B[JS][JS];
          Z[2][3] = ZERO;
          Z[3][3] = -E[JS][JS];
          Z[4][3] = ZERO;

          Z[1][4] = ZERO;
          Z[2][4] = -B[JS][JS];
          Z[3][4] = ZERO;
          Z[4][4] = -E[JS][JS];

          // Set up right hand side(s)

          RHS[1] = C[IS][JS];
          RHS[2] = C[ISP1][JS];
          RHS[3] = F[IS][JS];
          RHS[4] = F[ISP1][JS];

          // Solve Z * x = RHS

          dgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR);
          if (IERR.value > 0) INFO.value = IERR.value;
          if (IJOB == 0) {
            dgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
            if (SCALOC.value != ONE) {
              for (K = 1; K <= N; K++) {
                dscal(M, SCALOC.value, C(1, K).asArray(), 1);
                dscal(M, SCALOC.value, F(1, K).asArray(), 1);
              }
              SCALE.value *= SCALOC.value;
            }
          } else {
            dlatdf(IJOB, ZDIM, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV);
          }

          // Unpack solution vector(s)

          C[IS][JS] = RHS[1];
          C[ISP1][JS] = RHS[2];
          F[IS][JS] = RHS[3];
          F[ISP1][JS] = RHS[4];

          // Substitute R(I, J) and L(I, J) into remaining
          // equation.

          if (I > 1) {
            dgemv('N', IS - 1, MB, -ONE, A(1, IS), LDA, RHS(1), 1, ONE,
                C(1, JS).asArray(), 1);
            dgemv('N', IS - 1, MB, -ONE, D(1, IS), LDD, RHS(1), 1, ONE,
                F(1, JS).asArray(), 1);
          }
          if (J < Q) {
            dger(MB, N - JE, ONE, RHS(3), 1, B(JS, JE + 1).asArray(), LDB,
                C(IS, JE + 1), LDC);
            dger(MB, N - JE, ONE, RHS(3), 1, E(JS, JE + 1).asArray(), LDE,
                F(IS, JE + 1), LDF);
          }
        } else if ((MB == 2) && (NB == 2)) {
          // Build an 8-by-8 system Z * x = RHS

          dlaset('F', LDZ, LDZ, ZERO, ZERO, Z, LDZ);

          Z[1][1] = A[IS][IS];
          Z[2][1] = A[ISP1][IS];
          Z[5][1] = D[IS][IS];

          Z[1][2] = A[IS][ISP1];
          Z[2][2] = A[ISP1][ISP1];
          Z[5][2] = D[IS][ISP1];
          Z[6][2] = D[ISP1][ISP1];

          Z[3][3] = A[IS][IS];
          Z[4][3] = A[ISP1][IS];
          Z[7][3] = D[IS][IS];

          Z[3][4] = A[IS][ISP1];
          Z[4][4] = A[ISP1][ISP1];
          Z[7][4] = D[IS][ISP1];
          Z[8][4] = D[ISP1][ISP1];

          Z[1][5] = -B[JS][JS];
          Z[3][5] = -B[JS][JSP1];
          Z[5][5] = -E[JS][JS];
          Z[7][5] = -E[JS][JSP1];

          Z[2][6] = -B[JS][JS];
          Z[4][6] = -B[JS][JSP1];
          Z[6][6] = -E[JS][JS];
          Z[8][6] = -E[JS][JSP1];

          Z[1][7] = -B[JSP1][JS];
          Z[3][7] = -B[JSP1][JSP1];
          Z[7][7] = -E[JSP1][JSP1];

          Z[2][8] = -B[JSP1][JS];
          Z[4][8] = -B[JSP1][JSP1];
          Z[8][8] = -E[JSP1][JSP1];

          // Set up right hand side(s)

          K = 1;
          II = MB * NB + 1;
          for (JJ = 0; JJ <= NB - 1; JJ++) {
            dcopy(MB, C(IS, JS + JJ).asArray(), 1, RHS(K), 1);
            dcopy(MB, F(IS, JS + JJ).asArray(), 1, RHS(II), 1);
            K += MB;
            II += MB;
          }

          // Solve Z * x = RHS

          dgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR);
          if (IERR.value > 0) INFO.value = IERR.value;
          if (IJOB == 0) {
            dgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
            if (SCALOC.value != ONE) {
              for (K = 1; K <= N; K++) {
                dscal(M, SCALOC.value, C(1, K).asArray(), 1);
                dscal(M, SCALOC.value, F(1, K).asArray(), 1);
              }
              SCALE.value *= SCALOC.value;
            }
          } else {
            dlatdf(IJOB, ZDIM, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV);
          }

          // Unpack solution vector(s)

          K = 1;
          II = MB * NB + 1;
          for (JJ = 0; JJ <= NB - 1; JJ++) {
            dcopy(MB, RHS(K), 1, C(IS, JS + JJ).asArray(), 1);
            dcopy(MB, RHS(II), 1, F(IS, JS + JJ).asArray(), 1);
            K += MB;
            II += MB;
          }

          // Substitute R(I, J) and L(I, J) into remaining
          // equation.

          if (I > 1) {
            dgemm('N', 'N', IS - 1, NB, MB, -ONE, A(1, IS), LDA,
                RHS(1).asMatrix(MB), MB, ONE, C(1, JS), LDC);
            dgemm('N', 'N', IS - 1, NB, MB, -ONE, D(1, IS), LDD,
                RHS(1).asMatrix(MB), MB, ONE, F(1, JS), LDF);
          }
          if (J < Q) {
            K = MB * NB + 1;
            dgemm('N', 'N', MB, N - JE, NB, ONE, RHS(K).asMatrix(MB), MB,
                B(JS, JE + 1), LDB, ONE, C(IS, JE + 1), LDC);
            dgemm('N', 'N', MB, N - JE, NB, ONE, RHS(K).asMatrix(MB), MB,
                E(JS, JE + 1), LDE, ONE, F(IS, JE + 1), LDF);
          }
        }
      }
    }
  } else {
    // Solve (I, J) - subsystem
    // A[I][ I]**T * R(I, J) + D[I][ I]**T * L(J, J)  =  C[I][ J]
    // R(I, I)  * B[J][ J] + L(I, J)  * E[J][ J]  = -F[I][ J]
    // for I = 1, 2, ..., P, J = Q, Q - 1, ..., 1

    SCALE.value = ONE;
    SCALOC.value = ONE;
    for (I = 1; I <= P; I++) {
      IS = IWORK[I];
      ISP1 = IS + 1;
      IE = IWORK[I + 1] - 1;
      MB = IE - IS + 1;
      for (J = Q; J >= P + 2; J--) {
        JS = IWORK[J];
        JSP1 = JS + 1;
        JE = IWORK[J + 1] - 1;
        NB = JE - JS + 1;
        ZDIM = MB * NB * 2;
        if ((MB == 1) && (NB == 1)) {
          // Build a 2-by-2 system Z**T * x = RHS

          Z[1][1] = A[IS][IS];
          Z[2][1] = -B[JS][JS];
          Z[1][2] = D[IS][IS];
          Z[2][2] = -E[JS][JS];

          // Set up right hand side(s)

          RHS[1] = C[IS][JS];
          RHS[2] = F[IS][JS];

          // Solve Z**T * x = RHS

          dgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR);
          if (IERR.value > 0) INFO.value = IERR.value;

          dgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
          if (SCALOC.value != ONE) {
            for (K = 1; K <= N; K++) {
              dscal(M, SCALOC.value, C(1, K).asArray(), 1);
              dscal(M, SCALOC.value, F(1, K).asArray(), 1);
            }
            SCALE.value *= SCALOC.value;
          }

          // Unpack solution vector(s)

          C[IS][JS] = RHS[1];
          F[IS][JS] = RHS[2];

          // Substitute R(I, J) and L(I, J) into remaining
          // equation.

          if (J > P + 2) {
            ALPHA = RHS[1];
            daxpy(
                JS - 1, ALPHA, B(1, JS).asArray(), 1, F(IS, 1).asArray(), LDF);
            ALPHA = RHS[2];
            daxpy(
                JS - 1, ALPHA, E(1, JS).asArray(), 1, F(IS, 1).asArray(), LDF);
          }
          if (I < P) {
            ALPHA = -RHS[1];
            daxpy(M - IE, ALPHA, A(IS, IE + 1).asArray(), LDA,
                C(IE + 1, JS).asArray(), 1);
            ALPHA = -RHS[2];
            daxpy(M - IE, ALPHA, D(IS, IE + 1).asArray(), LDD,
                C(IE + 1, JS).asArray(), 1);
          }
        } else if ((MB == 1) && (NB == 2)) {
          // Build a 4-by-4 system Z**T * x = RHS

          Z[1][1] = A[IS][IS];
          Z[2][1] = ZERO;
          Z[3][1] = -B[JS][JS];
          Z[4][1] = -B[JSP1][JS];

          Z[1][2] = ZERO;
          Z[2][2] = A[IS][IS];
          Z[3][2] = -B[JS][JSP1];
          Z[4][2] = -B[JSP1][JSP1];

          Z[1][3] = D[IS][IS];
          Z[2][3] = ZERO;
          Z[3][3] = -E[JS][JS];
          Z[4][3] = ZERO;

          Z[1][4] = ZERO;
          Z[2][4] = D[IS][IS];
          Z[3][4] = -E[JS][JSP1];
          Z[4][4] = -E[JSP1][JSP1];

          // Set up right hand side(s)

          RHS[1] = C[IS][JS];
          RHS[2] = C[IS][JSP1];
          RHS[3] = F[IS][JS];
          RHS[4] = F[IS][JSP1];

          // Solve Z**T * x = RHS

          dgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR);
          if (IERR.value > 0) INFO.value = IERR.value;
          dgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
          if (SCALOC.value != ONE) {
            for (K = 1; K <= N; K++) {
              dscal(M, SCALOC.value, C(1, K).asArray(), 1);
              dscal(M, SCALOC.value, F(1, K).asArray(), 1);
            }
            SCALE.value *= SCALOC.value;
          }

          // Unpack solution vector(s)

          C[IS][JS] = RHS[1];
          C[IS][JSP1] = RHS[2];
          F[IS][JS] = RHS[3];
          F[IS][JSP1] = RHS[4];

          // Substitute R(I, J) and L(I, J) into remaining
          // equation.

          if (J > P + 2) {
            daxpy(
                JS - 1, RHS[1], B(1, JS).asArray(), 1, F(IS, 1).asArray(), LDF);
            daxpy(JS - 1, RHS[2], B(1, JSP1).asArray(), 1, F(IS, 1).asArray(),
                LDF);
            daxpy(
                JS - 1, RHS[3], E(1, JS).asArray(), 1, F(IS, 1).asArray(), LDF);
            daxpy(JS - 1, RHS[4], E(1, JSP1).asArray(), 1, F(IS, 1).asArray(),
                LDF);
          }
          if (I < P) {
            dger(M - IE, NB, -ONE, A(IS, IE + 1).asArray(), LDA, RHS(1), 1,
                C(IE + 1, JS), LDC);
            dger(M - IE, NB, -ONE, D(IS, IE + 1).asArray(), LDD, RHS(3), 1,
                C(IE + 1, JS), LDC);
          }
        } else if ((MB == 2) && (NB == 1)) {
          // Build a 4-by-4 system Z**T * x = RHS

          Z[1][1] = A[IS][IS];
          Z[2][1] = A[IS][ISP1];
          Z[3][1] = -B[JS][JS];
          Z[4][1] = ZERO;

          Z[1][2] = A[ISP1][IS];
          Z[2][2] = A[ISP1][ISP1];
          Z[3][2] = ZERO;
          Z[4][2] = -B[JS][JS];

          Z[1][3] = D[IS][IS];
          Z[2][3] = D[IS][ISP1];
          Z[3][3] = -E[JS][JS];
          Z[4][3] = ZERO;

          Z[1][4] = ZERO;
          Z[2][4] = D[ISP1][ISP1];
          Z[3][4] = ZERO;
          Z[4][4] = -E[JS][JS];

          // Set up right hand side(s)

          RHS[1] = C[IS][JS];
          RHS[2] = C[ISP1][JS];
          RHS[3] = F[IS][JS];
          RHS[4] = F[ISP1][JS];

          // Solve Z**T * x = RHS

          dgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR);
          if (IERR.value > 0) INFO.value = IERR.value;

          dgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
          if (SCALOC.value != ONE) {
            for (K = 1; K <= N; K++) {
              dscal(M, SCALOC.value, C(1, K).asArray(), 1);
              dscal(M, SCALOC.value, F(1, K).asArray(), 1);
            }
            SCALE.value *= SCALOC.value;
          }

          // Unpack solution vector(s)

          C[IS][JS] = RHS[1];
          C[ISP1][JS] = RHS[2];
          F[IS][JS] = RHS[3];
          F[ISP1][JS] = RHS[4];

          // Substitute R(I, J) and L(I, J) into remaining
          // equation.

          if (J > P + 2) {
            dger(MB, JS - 1, ONE, RHS(1), 1, B(1, JS).asArray(), 1, F(IS, 1),
                LDF);
            dger(MB, JS - 1, ONE, RHS(3), 1, E(1, JS).asArray(), 1, F(IS, 1),
                LDF);
          }
          if (I < P) {
            dgemv('T', MB, M - IE, -ONE, A(IS, IE + 1), LDA, RHS(1), 1, ONE,
                C(IE + 1, JS).asArray(), 1);
            dgemv('T', MB, M - IE, -ONE, D(IS, IE + 1), LDD, RHS(3), 1, ONE,
                C(IE + 1, JS).asArray(), 1);
          }
        } else if ((MB == 2) && (NB == 2)) {
          // Build an 8-by-8 system Z**T * x = RHS

          dlaset('F', LDZ, LDZ, ZERO, ZERO, Z, LDZ);

          Z[1][1] = A[IS][IS];
          Z[2][1] = A[IS][ISP1];
          Z[5][1] = -B[JS][JS];
          Z[7][1] = -B[JSP1][JS];

          Z[1][2] = A[ISP1][IS];
          Z[2][2] = A[ISP1][ISP1];
          Z[6][2] = -B[JS][JS];
          Z[8][2] = -B[JSP1][JS];

          Z[3][3] = A[IS][IS];
          Z[4][3] = A[IS][ISP1];
          Z[5][3] = -B[JS][JSP1];
          Z[7][3] = -B[JSP1][JSP1];

          Z[3][4] = A[ISP1][IS];
          Z[4][4] = A[ISP1][ISP1];
          Z[6][4] = -B[JS][JSP1];
          Z[8][4] = -B[JSP1][JSP1];

          Z[1][5] = D[IS][IS];
          Z[2][5] = D[IS][ISP1];
          Z[5][5] = -E[JS][JS];

          Z[2][6] = D[ISP1][ISP1];
          Z[6][6] = -E[JS][JS];

          Z[3][7] = D[IS][IS];
          Z[4][7] = D[IS][ISP1];
          Z[5][7] = -E[JS][JSP1];
          Z[7][7] = -E[JSP1][JSP1];

          Z[4][8] = D[ISP1][ISP1];
          Z[6][8] = -E[JS][JSP1];
          Z[8][8] = -E[JSP1][JSP1];

          // Set up right hand side(s)

          K = 1;
          II = MB * NB + 1;
          for (JJ = 0; JJ <= NB - 1; JJ++) {
            dcopy(MB, C(IS, JS + JJ).asArray(), 1, RHS(K), 1);
            dcopy(MB, F(IS, JS + JJ).asArray(), 1, RHS(II), 1);
            K += MB;
            II += MB;
          }

          // Solve Z**T * x = RHS

          dgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR);
          if (IERR.value > 0) INFO.value = IERR.value;

          dgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC);
          if (SCALOC.value != ONE) {
            for (K = 1; K <= N; K++) {
              dscal(M, SCALOC.value, C(1, K).asArray(), 1);
              dscal(M, SCALOC.value, F(1, K).asArray(), 1);
            }
            SCALE.value *= SCALOC.value;
          }

          // Unpack solution vector(s)

          K = 1;
          II = MB * NB + 1;
          for (JJ = 0; JJ <= NB - 1; JJ++) {
            dcopy(MB, RHS(K), 1, C(IS, JS + JJ).asArray(), 1);
            dcopy(MB, RHS(II), 1, F(IS, JS + JJ).asArray(), 1);
            K += MB;
            II += MB;
          }

          // Substitute R(I, J) and L(I, J) into remaining
          // equation.

          if (J > P + 2) {
            dgemm('N', 'T', MB, JS - 1, NB, ONE, C(IS, JS), LDC, B(1, JS), LDB,
                ONE, F(IS, 1), LDF);
            dgemm('N', 'T', MB, JS - 1, NB, ONE, F(IS, JS), LDF, E(1, JS), LDE,
                ONE, F(IS, 1), LDF);
          }
          if (I < P) {
            dgemm('T', 'N', M - IE, NB, MB, -ONE, A(IS, IE + 1), LDA, C(IS, JS),
                LDC, ONE, C(IE + 1, JS), LDC);
            dgemm('T', 'N', M - IE, NB, MB, -ONE, D(IS, IE + 1), LDD, F(IS, JS),
                LDF, ONE, C(IE + 1, JS), LDC);
          }
        }
      }
    }
  }
}
