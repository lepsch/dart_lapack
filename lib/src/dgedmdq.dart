// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgedmd.dart';
import 'package:dart_lapack/src/dgeqrf.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorgqr.dart';
import 'package:dart_lapack/src/dormqr.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgedmdq(
  final String JOBS,
  final String JOBZ,
  final String JOBR,
  final String JOBQ,
  final String JOBT,
  final String JOBF,
  final int WHTSVD,
  final int M,
  final int N,
  final Matrix<double> F_,
  final int LDF,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> Y_,
  final int LDY,
  final int NRNK,
  final double TOL,
  final Box<int> K,
  final Array<double> REIG_,
  final Array<double> IMEIG_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> RES_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> V_,
  final int LDV,
  final Matrix<double> S_,
  final int LDS,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
  final F = F_.having(ld: LDF),
      X = X_.having(ld: LDX),
      Y = Y_.having(ld: LDY),
      Z = Z_.having(ld: LDZ),
      B = B_.having(ld: LDB),
      V = V_.having(ld: LDV),
      S = S_.having(ld: LDS),
      REIG = REIG_.having(),
      IMEIG = IMEIG_.having(),
      RES = RES_.having(),
      WORK = WORK_.having(),
      IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int IMINWR = 0,
      MLWDMD,
      MLWGQR,
      MLWMQR,
      MLWORK = 0,
      MLWQR,
      MINMN,
      OLWDMD,
      OLWGQR,
      OLWMQR,
      OLWORK = 0,
      OLWQR;
  bool LQUERY,
      SCCOLX,
      SCCOLY,
      WANTQ,
      WNTTRF,
      WNTRES,
      WNTVEC,
      WNTVCF,
      WNTVCQ,
      WNTREF,
      WNTEX;
  String JOBVL;
  final RDUMMY = Array<double>(2);
  final INFO1 = Box(0);

  // Test the input arguments
  WNTRES = lsame(JOBR, 'R');
  SCCOLX = lsame(JOBS, 'S') || lsame(JOBS, 'C');
  SCCOLY = lsame(JOBS, 'Y');
  WNTVEC = lsame(JOBZ, 'V');
  WNTVCF = lsame(JOBZ, 'F');
  WNTVCQ = lsame(JOBZ, 'Q');
  WNTREF = lsame(JOBF, 'R');
  WNTEX = lsame(JOBF, 'E');
  WANTQ = lsame(JOBQ, 'Q');
  WNTTRF = lsame(JOBT, 'R');
  MINMN = min(M, N);
  INFO.value = 0;
  LQUERY = ((LWORK == -1) || (LIWORK == -1));

  if (!(SCCOLX || SCCOLY || lsame(JOBS, 'N'))) {
    INFO.value = -1;
  } else if (!(WNTVEC || WNTVCF || WNTVCQ || lsame(JOBZ, 'N'))) {
    INFO.value = -2;
  } else if (!(WNTRES || lsame(JOBR, 'N')) || (WNTRES && lsame(JOBZ, 'N'))) {
    INFO.value = -3;
  } else if (!(WANTQ || lsame(JOBQ, 'N'))) {
    INFO.value = -4;
  } else if (!(WNTTRF || lsame(JOBT, 'N'))) {
    INFO.value = -5;
  } else if (!(WNTREF || WNTEX || lsame(JOBF, 'N'))) {
    INFO.value = -6;
  } else if (!((WHTSVD == 1) ||
      (WHTSVD == 2) ||
      (WHTSVD == 3) ||
      (WHTSVD == 4))) {
    INFO.value = -7;
  } else if (M < 0) {
    INFO.value = -8;
  } else if ((N < 0) || (N > M + 1)) {
    INFO.value = -9;
  } else if (LDF < M) {
    INFO.value = -11;
  } else if (LDX < MINMN) {
    INFO.value = -13;
  } else if (LDY < MINMN) {
    INFO.value = -15;
  } else if (!((NRNK == -2) || (NRNK == -1) || ((NRNK >= 1) && (NRNK <= N)))) {
    INFO.value = -16;
  } else if ((TOL < ZERO) || (TOL >= ONE)) {
    INFO.value = -17;
  } else if (LDZ < M) {
    INFO.value = -22;
  } else if ((WNTREF || WNTEX) && (LDB < MINMN)) {
    INFO.value = -25;
  } else if (LDV < N - 1) {
    INFO.value = -27;
  } else if (LDS < N - 1) {
    INFO.value = -29;
  }

  if (WNTVEC || WNTVCF || WNTVCQ) {
    JOBVL = 'V';
  } else {
    JOBVL = 'N';
  }
  if (INFO.value == 0) {
    // Compute the minimal and the optimal workspace
    // requirements. Simulate running the code and
    // determine minimal and optimal sizes of the
    // workspace at any moment of the run.
    if ((N == 0) || (N == 1)) {
      // All output except K is void. INFO=1 signals
      // the void input. In case of a workspace query,
      // the minimal workspace lengths are returned.
      if (LQUERY) {
        IWORK[1] = 1;
        WORK[1] = 2;
        WORK[2] = 2;
      } else {
        K.value = 0;
      }
      INFO.value = 1;
      return;
    }
    MLWQR = max(1, N); // Minimal workspace length for DGEQRF.
    MLWORK = MINMN + MLWQR;
    if (LQUERY) {
      dgeqrf(M, N, F, LDF, WORK, RDUMMY, -1, INFO1);
      OLWQR = RDUMMY[1].toInt();
      OLWORK = min(M, N) + OLWQR;
    }
    dgedmd(
        JOBS,
        JOBVL,
        JOBR,
        JOBF,
        WHTSVD,
        MINMN,
        N - 1,
        X,
        LDX,
        Y,
        LDY,
        NRNK,
        TOL,
        K,
        REIG,
        IMEIG,
        Z,
        LDZ,
        RES,
        B,
        LDB,
        V,
        LDV,
        S,
        LDS,
        WORK,
        -1,
        IWORK,
        LIWORK,
        INFO1);
    MLWDMD = WORK[1].toInt();
    MLWORK = max(MLWORK, MINMN + MLWDMD);
    IMINWR = IWORK[1];
    if (LQUERY) {
      OLWDMD = WORK[2].toInt();
      OLWORK = max(OLWORK, MINMN + OLWDMD);
    }
    if (WNTVEC || WNTVCF) {
      MLWMQR = max(1, N);
      MLWORK = max(MLWORK, MINMN + N - 1 + MLWMQR);
      if (LQUERY) {
        dormqr('L', 'N', M, N, MINMN, F, LDF, WORK, Z, LDZ, WORK, -1, INFO1);
        OLWMQR = WORK[1].toInt();
        OLWORK = max(OLWORK, MINMN + N - 1 + OLWMQR);
      }
    }
    if (WANTQ) {
      MLWGQR = N;
      MLWORK = max(MLWORK, MINMN + N - 1 + MLWGQR);
      if (LQUERY) {
        dorgqr(M, MINMN, MINMN, F, LDF, WORK, WORK, -1, INFO1);
        OLWGQR = WORK[1].toInt();
        OLWORK = max(OLWORK, MINMN + N - 1 + OLWGQR);
      }
    }
    IMINWR = max(1, IMINWR);
    MLWORK = max(2, MLWORK);
    if (LWORK < MLWORK && !LQUERY) INFO.value = -31;
    if (LIWORK < IMINWR && !LQUERY) INFO.value = -33;
  }
  if (INFO.value != 0) {
    xerbla('DGEDMDQ', -INFO.value);
    return;
  } else if (LQUERY) {
    // Return minimal and optimal workspace sizes
    IWORK[1] = IMINWR;
    WORK[1] = MLWORK.toDouble();
    WORK[2] = OLWORK.toDouble();
    return;
  }

  // Initial QR factorization that is used to represent the
  // snapshots as elements of lower dimensional subspace.
  // For large scale computation with M >>N , at this place
  // one can use an out of core QRF.

  dgeqrf(M, N, F, LDF, WORK, WORK(MINMN + 1), LWORK - MINMN, INFO1);

  // Define X and Y as the snapshots representations in the
  // orthogonal basis computed in the QR factorization.
  // X corresponds to the leading N-1 and Y to the trailing
  // N-1 snapshots.
  dlaset('L', MINMN, N - 1, ZERO, ZERO, X, LDX);
  dlacpy('U', MINMN, N - 1, F, LDF, X, LDX);
  dlacpy('A', MINMN, N - 1, F(1, 2), LDF, Y, LDY);
  if (M >= 3) {
    dlaset('L', MINMN - 2, N - 2, ZERO, ZERO, Y(3, 1), LDY);
  }

  // Compute the DMD of the projected snapshot pairs (X,Y)
  dgedmd(
      JOBS,
      JOBVL,
      JOBR,
      JOBF,
      WHTSVD,
      MINMN,
      N - 1,
      X,
      LDX,
      Y,
      LDY,
      NRNK,
      TOL,
      K,
      REIG,
      IMEIG,
      Z,
      LDZ,
      RES,
      B,
      LDB,
      V,
      LDV,
      S,
      LDS,
      WORK(MINMN + 1),
      LWORK - MINMN,
      IWORK,
      LIWORK,
      INFO1);
  if (INFO1.value == 2 || INFO1.value == 3) {
    // Return with error code. See DGEDMD for details.
    INFO.value = INFO1.value;
    return;
  } else {
    INFO.value = INFO1.value;
  }

  // The Ritz vectors (Koopman modes) can be explicitly
  // formed or returned in factored form.
  if (WNTVEC) {
    // Compute the eigenvectors explicitly.
    if (M > MINMN) {
      dlaset('A', M - MINMN, K.value, ZERO, ZERO, Z(MINMN + 1, 1), LDZ);
    }
    dormqr('L', 'N', M, K.value, MINMN, F, LDF, WORK, Z, LDZ, WORK(MINMN + N),
        LWORK - (MINMN + N - 1), INFO1);
  } else if (WNTVCF) {
    // Return the Ritz vectors (eigenvectors) in factored
    // form Z*V, where Z contains orthonormal matrix (the
    // product of Q from the initial QR factorization and
    // the SVD/POD_basis returned by DGEDMD in X) and the
    // second factor (the eigenvectors of the Rayleigh
    // quotient) is in the array V, as returned by DGEDMD.
    dlacpy('A', N, K.value, X, LDX, Z, LDZ);
    if (M > N) dlaset('A', M - N, K.value, ZERO, ZERO, Z(N + 1, 1), LDZ);
    dormqr('L', 'N', M, K.value, MINMN, F, LDF, WORK, Z, LDZ, WORK(MINMN + N),
        LWORK - (MINMN + N - 1), INFO1);
  }

  // Some optional output variables:
  //
  // The upper triangular factor R in the initial QR
  // factorization is optionally returned in the array Y.
  // This is useful if this to DGEDMDQ is to be
  // followed by a streaming DMD that is implemented in a
  // QR compressed form.
  if (WNTTRF) {
    // Return the upper triangular R in Y
    dlaset('A', MINMN, N, ZERO, ZERO, Y, LDY);
    dlacpy('U', MINMN, N, F, LDF, Y, LDY);
  }

  //  The orthonormal/orthogonal factor Q in the initial QR
  //  factorization is optionally returned in the array F.
  //  Same as with the triangular factor above, this is
  //  useful in a streaming DMD.
  if (WANTQ) {
    // Q overwrites F
    dorgqr(M, MINMN, MINMN, F, LDF, WORK, WORK(MINMN + N),
        LWORK - (MINMN + N - 1), INFO1);
  }
}
