// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgedmd.dart';
import 'package:lapack/src/zgeqrf.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zungqr.dart';
import 'package:lapack/src/zunmqr.dart';

void zgedmdq(
  final String JOBS,
  final String JOBZ,
  final String JOBR,
  final String JOBQ,
  final String JOBT,
  final String JOBF,
  final int WHTSVD,
  final int M,
  final int N,
  final Matrix<Complex> F_,
  final int LDF,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> Y_,
  final int LDY,
  final int NRNK,
  final double TOL,
  final Box<int> K,
  final Array<Complex> EIGS_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<double> RES_,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> S_,
  final int LDS,
  final Array<Complex> ZWORK_,
  final int LZWORK,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
//  -- LAPACK driver routine                                           --

//  -- LAPACK is a software package provided by University of          --
//  -- Tennessee, University of California Berkeley, University of     --
//  -- Colorado Denver and NAG Ltd..                                   --
  final F = F_.having(ld: LDF);
  final X = X_.having(ld: LDX);
  final Y = Y_.having(ld: LDY);
  final EIGS = EIGS_.having();
  final Z = Z_.having(ld: LDZ);
  final RES = RES_.having();
  final B = B_.having(ld: LDB);
  final V = V_.having(ld: LDV);
  final S = S_.having(ld: LDS);
  final ZWORK = ZWORK_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();

  const ONE = 1.0, ZERO = 0.0;
  int IMINWR = 0,
      MINMN,
      MLRWRK = 0,
      MLWDMD,
      MLWGQR,
      MLWMQR,
      MLWORK = 0,
      MLWQR,
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
  LQUERY = ((LZWORK == -1) || (LWORK == -1) || (LIWORK == -1));

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
    INFO.value = -21;
  } else if ((WNTREF || WNTEX) && (LDB < MINMN)) {
    INFO.value = -24;
  } else if (LDV < N - 1) {
    INFO.value = -26;
  } else if (LDS < N - 1) {
    INFO.value = -28;
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
        ZWORK[1] = 2.toComplex();
        ZWORK[2] = 2.toComplex();
        WORK[1] = 2;
        WORK[2] = 2;
      } else {
        K.value = 0;
      }
      INFO.value = 1;
      return;
    }

    MLRWRK = 2;
    MLWORK = 2;
    OLWORK = 2;
    IMINWR = 1;
    MLWQR = max(1, N); // Minimal workspace length for zgeqrf.
    MLWORK = max(MLWORK, MINMN + MLWQR);

    if (LQUERY) {
      zgeqrf(M, N, F, LDF, ZWORK, ZWORK, -1, INFO1);
      OLWQR = ZWORK[1].toInt();
      OLWORK = max(OLWORK, MINMN + OLWQR);
    }
    zgedmd(
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
        EIGS,
        Z,
        LDZ,
        RES,
        B,
        LDB,
        V,
        LDV,
        S,
        LDS,
        ZWORK,
        -1,
        WORK,
        -1,
        IWORK,
        -1,
        INFO1);
    MLWDMD = ZWORK[1].toInt();
    MLWORK = max(MLWORK, MINMN + MLWDMD);
    MLRWRK = max(MLRWRK, WORK[1].toInt());
    IMINWR = max(IMINWR, IWORK[1]);
    if (LQUERY) {
      OLWDMD = ZWORK[2].toInt();
      OLWORK = max(OLWORK, MINMN + OLWDMD);
    }
    if (WNTVEC || WNTVCF) {
      MLWMQR = max(1, N);
      MLWORK = max(MLWORK, MINMN + MLWMQR);
      if (LQUERY) {
        zunmqr('L', 'N', M, N, MINMN, F, LDF, ZWORK, Z, LDZ, ZWORK, -1, INFO1);
        OLWMQR = ZWORK[1].toInt();
        OLWORK = max(OLWORK, MINMN + OLWMQR);
      }
    }
    if (WANTQ) {
      MLWGQR = max(1, N);
      MLWORK = max(MLWORK, MINMN + MLWGQR);
      if (LQUERY) {
        zungqr(M, MINMN, MINMN, F, LDF, ZWORK, ZWORK, -1, INFO1);
        OLWGQR = ZWORK[1].toInt();
        OLWORK = max(OLWORK, MINMN + OLWGQR);
      }
    }
    if (LIWORK < IMINWR && (!LQUERY)) INFO.value = -34;
    if (LWORK < MLRWRK && (!LQUERY)) INFO.value = -32;
    if (LZWORK < MLWORK && (!LQUERY)) INFO.value = -30;
  }
  if (INFO.value != 0) {
    xerbla('ZGEDMDQ', -INFO.value);
    return;
  } else if (LQUERY) {
    // Return minimal and optimal workspace sizes
    IWORK[1] = IMINWR;
    ZWORK[1] = MLWORK.toComplex();
    ZWORK[2] = OLWORK.toComplex();
    WORK[1] = MLRWRK.toDouble();
    WORK[2] = MLRWRK.toDouble();
    return;
  }

  // Initial QR factorization that is used to represent the
  // snapshots as elements of lower dimensional subspace.
  // For large scale computation with M >> N, at this place
  // one can use an out of core QRF.

  zgeqrf(M, N, F, LDF, ZWORK, ZWORK(MINMN + 1), LZWORK - MINMN, INFO1);

  // Define X and Y as the snapshots representations in the
  // orthogonal basis computed in the QR factorization.
  // X corresponds to the leading N-1 and Y to the trailing
  // N-1 snapshots.
  zlaset('L', MINMN, N - 1, Complex.zero, Complex.zero, X, LDX);
  zlacpy('U', MINMN, N - 1, F, LDF, X, LDX);
  zlacpy('A', MINMN, N - 1, F(1, 2), LDF, Y, LDY);
  if (M >= 3) {
    zlaset('L', MINMN - 2, N - 2, Complex.zero, Complex.zero, Y(3, 1), LDY);
  }

  // Compute the DMD of the projected snapshot pairs (X,Y)
  zgedmd(
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
      EIGS,
      Z,
      LDZ,
      RES,
      B,
      LDB,
      V,
      LDV,
      S,
      LDS,
      ZWORK(MINMN + 1),
      LZWORK - MINMN,
      WORK,
      LWORK,
      IWORK,
      LIWORK,
      INFO1);
  if (INFO1.value == 2 || INFO1.value == 3) {
    // Return with error code. See zgedmd for details.
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
      zlaset('A', M - MINMN, K.value, Complex.zero, Complex.zero,
          Z(MINMN + 1, 1), LDZ);
    }
    zunmqr('L', 'N', M, K.value, MINMN, F, LDF, ZWORK, Z, LDZ, ZWORK(MINMN + 1),
        LZWORK - MINMN, INFO1);
  } else if (WNTVCF) {
    // Return the Ritz vectors (eigenvectors) in factored
    // form Z*V, where Z contains orthonormal matrix (the
    // product of Q from the initial QR factorization and
    // the SVD/POD_basis returned by zgedmd in X) and the
    // second factor (the eigenvectors of the Rayleigh
    // quotient) is in the array V, as returned by zgedmd.
    zlacpy('A', N, K.value, X, LDX, Z, LDZ);
    if (M > N) {
      zlaset('A', M - N, K.value, Complex.zero, Complex.zero, Z(N + 1, 1), LDZ);
    }
    zunmqr('L', 'N', M, K.value, MINMN, F, LDF, ZWORK, Z, LDZ, ZWORK(MINMN + 1),
        LZWORK - MINMN, INFO1);
  }

  // Some optional output variables:

  // The upper triangular factor R in the initial QR
  // factorization is optionally returned in the array Y.
  // This is useful if this to ZGEDMDQ is to be
  // followed by a streaming DMD that is implemented in a
  // QR compressed form.
  if (WNTTRF) {
    // Return the upper triangular R in Y
    zlaset('A', MINMN, N, Complex.zero, Complex.zero, Y, LDY);
    zlacpy('U', MINMN, N, F, LDF, Y, LDY);
  }

  // The orthonormal/unitary factor Q in the initial QR
  // factorization is optionally returned in the array F.
  // Same as with the triangular factor above, this is
  // useful in a streaming DMD.
  if (WANTQ) {
    // Q overwrites F
    zungqr(M, MINMN, MINMN, F, LDF, ZWORK, ZWORK(MINMN + 1), LZWORK - MINMN,
        INFO1);
  }
}
