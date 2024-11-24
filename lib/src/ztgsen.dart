// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacn2.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlassq.dart';
import 'package:dart_lapack/src/ztgexc.dart';
import 'package:dart_lapack/src/ztgsyl.dart';

void ztgsen(
  final int IJOB,
  final bool WANTQ,
  final bool WANTZ,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Box<int> M,
  final Box<double> PL,
  final Box<double> PR,
  final Array<double> DIF_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final SELECT = SELECT_.having();
  final DIF = DIF_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  const IDIFJB = 3;
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, SWAP, WANTD, WANTD1, WANTD2, WANTP;
  int I, IJB, K, LIWMIN, LWMIN, MN2, N1, N2;
  double SAFMIN;
  Complex TEMP1, TEMP2;
  final ISAVE = Array<int>(3);
  final DSCALE = Box(0.0), DSUM = Box(0.0), RDSCAL = Box(0.0);
  final IERR = Box(0), KS = Box(0), KASE = Box(0);

  // Decode and test the input parameters

  INFO.value = 0;
  LQUERY = (LWORK == -1 || LIWORK == -1);

  if (IJOB < 0 || IJOB > 5) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDQ < 1 || (WANTQ && LDQ < N)) {
    INFO.value = -13;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -15;
  }

  if (INFO.value != 0) {
    xerbla('ZTGSEN', -INFO.value);
    return;
  }

  IERR.value = 0;

  WANTP = IJOB == 1 || IJOB >= 4;
  WANTD1 = IJOB == 2 || IJOB == 4;
  WANTD2 = IJOB == 3 || IJOB == 5;
  WANTD = WANTD1 || WANTD2;

  // Set M to the dimension of the specified pair of deflating
  // subspaces.

  M.value = 0;
  if (!LQUERY || IJOB != 0) {
    for (K = 1; K <= N; K++) {
      ALPHA[K] = A[K][K];
      BETA[K] = B[K][K];
      if (K < N) {
        if (SELECT[K]) M.value++;
      } else {
        if (SELECT[N]) M.value++;
      }
    }
  }

  if (IJOB == 1 || IJOB == 2 || IJOB == 4) {
    LWMIN = max(1, 2 * M.value * (N - M.value));
    LIWMIN = max(1, N + 2);
  } else if (IJOB == 3 || IJOB == 5) {
    LWMIN = max(1, 4 * M.value * (N - M.value));
    LIWMIN = max(1, max(2 * M.value * (N - M.value), N + 2));
  } else {
    LWMIN = 1;
    LIWMIN = 1;
  }

  WORK[1] = LWMIN.toComplex();
  IWORK[1] = LIWMIN;

  if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -21;
  } else if (LIWORK < LIWMIN && !LQUERY) {
    INFO.value = -23;
  }

  if (INFO.value != 0) {
    xerbla('ZTGSEN', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible.

  if (M.value == N || M.value == 0) {
    if (WANTP) {
      PL.value = ONE;
      PR.value = ONE;
    }
    if (WANTD) {
      DSCALE.value = ZERO;
      DSUM.value = ONE;
      for (I = 1; I <= N; I++) {
        zlassq(N, A(1, I).asArray(), 1, DSCALE, DSUM);
        zlassq(N, B(1, I).asArray(), 1, DSCALE, DSUM);
      }
      DIF[1] = DSCALE.value * sqrt(DSUM.value);
      DIF[2] = DIF[1];
    }
    WORK[1] = LWMIN.toComplex();
    IWORK[1] = LIWMIN;
    return;
  }

  // Get machine constant

  SAFMIN = dlamch('S');

  // Collect the selected blocks at the top-left corner of (A, B).

  KS.value = 0;
  for (K = 1; K <= N; K++) {
    SWAP = SELECT[K];
    if (SWAP) {
      KS.value++;

      // Swap the K-th block to position KS. Compute unitary Q
      // and Z that will swap adjacent diagonal blocks in (A, B).

      if (K != KS.value) {
        ztgexc(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, K, KS, IERR);
      }

      if (IERR.value > 0) {
        // Swap is rejected: exit.

        INFO.value = 1;
        if (WANTP) {
          PL.value = ZERO;
          PR.value = ZERO;
        }
        if (WANTD) {
          DIF[1] = ZERO;
          DIF[2] = ZERO;
        }
        WORK[1] = LWMIN.toComplex();
        IWORK[1] = LIWMIN;
        return;
      }
    }
  }
  if (WANTP) {
    // Solve generalized Sylvester equation for R and L:
    //            A11 * R - L * A22 = A12
    //            B11 * R - L * B22 = B12

    N1 = M.value;
    N2 = N - M.value;
    I = N1 + 1;
    zlacpy('Full', N1, N2, A(1, I), LDA, WORK.asMatrix(N1), N1);
    zlacpy('Full', N1, N2, B(1, I), LDB, WORK(N1 * N2 + 1).asMatrix(N1), N1);
    IJB = 0;
    ztgsyl(
        'N',
        IJB,
        N1,
        N2,
        A,
        LDA,
        A(I, I),
        LDA,
        WORK.asMatrix(N1),
        N1,
        B,
        LDB,
        B(I, I),
        LDB,
        WORK(N1 * N2 + 1).asMatrix(N1),
        N1,
        DSCALE,
        DIF(1),
        WORK(N1 * N2 * 2 + 1),
        LWORK - 2 * N1 * N2,
        IWORK,
        IERR);

    // Estimate the reciprocal of norms of "projections" onto
    // left and right eigenspaces

    RDSCAL.value = ZERO;
    DSUM.value = ONE;
    zlassq(N1 * N2, WORK, 1, RDSCAL, DSUM);
    PL.value = RDSCAL.value * sqrt(DSUM.value);
    if (PL.value == ZERO) {
      PL.value = ONE;
    } else {
      PL.value = DSCALE.value /
          (sqrt(DSCALE.value * DSCALE.value / PL.value + PL.value) *
              sqrt(PL.value));
    }
    RDSCAL.value = ZERO;
    DSUM.value = ONE;
    zlassq(N1 * N2, WORK(N1 * N2 + 1), 1, RDSCAL, DSUM);
    PR.value = RDSCAL.value * sqrt(DSUM.value);
    if (PR.value == ZERO) {
      PR.value = ONE;
    } else {
      PR.value = DSCALE.value /
          (sqrt(DSCALE.value * DSCALE.value / PR.value + PR.value) *
              sqrt(PR.value));
    }
  }
  if (WANTD) {
    // Compute estimates Difu and Difl.

    if (WANTD1) {
      N1 = M.value;
      N2 = N - M.value;
      I = N1 + 1;
      IJB = IDIFJB;

      // Frobenius norm-based Difu estimate.

      ztgsyl(
          'N',
          IJB,
          N1,
          N2,
          A,
          LDA,
          A(I, I),
          LDA,
          WORK.asMatrix(N1),
          N1,
          B,
          LDB,
          B(I, I),
          LDB,
          WORK(N1 * N2 + 1).asMatrix(N1),
          N1,
          DSCALE,
          DIF(1),
          WORK(N1 * N2 * 2 + 1),
          LWORK - 2 * N1 * N2,
          IWORK,
          IERR);

      // Frobenius norm-based Difl estimate.

      ztgsyl(
          'N',
          IJB,
          N2,
          N1,
          A(I, I),
          LDA,
          A,
          LDA,
          WORK.asMatrix(N2),
          N2,
          B(I, I),
          LDB,
          B,
          LDB,
          WORK(N1 * N2 + 1).asMatrix(N2),
          N2,
          DSCALE,
          DIF(2),
          WORK(N1 * N2 * 2 + 1),
          LWORK - 2 * N1 * N2,
          IWORK,
          IERR);
    } else {
      // Compute 1-norm-based estimates of Difu and Difl using
      // reversed communication with ZLACN2. In each step a
      // generalized Sylvester equation or a transposed variant
      // is solved.

      KASE.value = 0;
      N1 = M.value;
      N2 = N - M.value;
      I = N1 + 1;
      IJB = 0;
      MN2 = 2 * N1 * N2;

      // 1-norm-based estimate of Difu.

      while (true) {
        zlacn2(MN2, WORK(MN2 + 1), WORK, DIF(1), KASE, ISAVE);
        if (KASE.value == 0) break;
        if (KASE.value == 1) {
          // Solve generalized Sylvester equation

          ztgsyl(
              'N',
              IJB,
              N1,
              N2,
              A,
              LDA,
              A(I, I),
              LDA,
              WORK.asMatrix(N1),
              N1,
              B,
              LDB,
              B(I, I),
              LDB,
              WORK(N1 * N2 + 1).asMatrix(N1),
              N1,
              DSCALE,
              DIF(1),
              WORK(N1 * N2 * 2 + 1),
              LWORK - 2 * N1 * N2,
              IWORK,
              IERR);
        } else {
          // Solve the transposed variant.

          ztgsyl(
              'C',
              IJB,
              N1,
              N2,
              A,
              LDA,
              A(I, I),
              LDA,
              WORK.asMatrix(N1),
              N1,
              B,
              LDB,
              B(I, I),
              LDB,
              WORK(N1 * N2 + 1).asMatrix(N1),
              N1,
              DSCALE,
              DIF(1),
              WORK(N1 * N2 * 2 + 1),
              LWORK - 2 * N1 * N2,
              IWORK,
              IERR);
        }
      }
      DIF[1] = DSCALE.value / DIF[1];

      // 1-norm-based estimate of Difl.

      while (true) {
        zlacn2(MN2, WORK(MN2 + 1), WORK, DIF(2), KASE, ISAVE);
        if (KASE.value == 0) break;
        if (KASE.value == 1) {
          // Solve generalized Sylvester equation

          ztgsyl(
              'N',
              IJB,
              N2,
              N1,
              A(I, I),
              LDA,
              A,
              LDA,
              WORK.asMatrix(N2),
              N2,
              B(I, I),
              LDB,
              B,
              LDB,
              WORK(N1 * N2 + 1).asMatrix(N2),
              N2,
              DSCALE,
              DIF(2),
              WORK(N1 * N2 * 2 + 1),
              LWORK - 2 * N1 * N2,
              IWORK,
              IERR);
        } else {
          // Solve the transposed variant.

          ztgsyl(
              'C',
              IJB,
              N2,
              N1,
              A(I, I),
              LDA,
              A,
              LDA,
              WORK.asMatrix(N2),
              N2,
              B,
              LDB,
              B(I, I),
              LDB,
              WORK(N1 * N2 + 1).asMatrix(N2),
              N2,
              DSCALE,
              DIF(2),
              WORK(N1 * N2 * 2 + 1),
              LWORK - 2 * N1 * N2,
              IWORK,
              IERR);
        }
      }
      DIF[2] = DSCALE.value / DIF[2];
    }
  }

  // If B(K,K) is complex, make it real and positive (normalization
  // of the generalized Schur form) and Store the generalized
  // eigenvalues of reordered pair (A, B)

  for (K = 1; K <= N; K++) {
    DSCALE.value = B[K][K].abs();
    if (DSCALE.value > SAFMIN) {
      TEMP1 = (B[K][K] / DSCALE.value.toComplex()).conjugate();
      TEMP2 = B[K][K] / DSCALE.value.toComplex();
      B[K][K] = DSCALE.value.toComplex();
      zscal(N - K, TEMP1, B(K, K + 1).asArray(), LDB);
      zscal(N - K + 1, TEMP1, A(K, K).asArray(), LDA);
      if (WANTQ) zscal(N, TEMP2, Q(1, K).asArray(), 1);
    } else {
      B[K][K] = Complex.zero;
    }

    ALPHA[K] = A[K][K];
    BETA[K] = B[K][K];
  }

  WORK[1] = LWMIN.toComplex();
  IWORK[1] = LIWMIN;
}
