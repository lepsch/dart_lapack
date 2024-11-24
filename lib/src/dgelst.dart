// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgelqt.dart';
import 'package:lapack/src/dgemlqt.dart';
import 'package:lapack/src/dgemqrt.dart';
import 'package:lapack/src/dgeqrt.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dtrtrs.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgelst(
  final String TRANS,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, TPSD = false;
  int BROW,
      I,
      IASCL,
      IBSCL,
      J,
      LWOPT = 0,
      MN,
      MNNRHS = 0,
      NB = 0,
      NBMIN,
      SCLLEN;
  double ANRM, BIGNUM, BNRM, SMLNUM;
  final RWORK = Array<double>(1);
  // ..
  // .. External Functions ..
  //- bool               lsame;
  //- int                ILAENV;
  //- double             DLAMCH, DLANGE;
  // EXTERNAL lsame, ILAENV, DLAMCH, DLANGE
  // ..
  // .. External Subroutines ..
  // EXTERNAL DGELQT, DGEQRT, DGEMLQT, DGEMQRT, DLASCL, DLASET, DTRTRS, XERBLA
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC DBLE, MAX, MIN

  // Test the input arguments.

  INFO.value = 0;
  MN = min(M, N);
  LQUERY = (LWORK == -1);
  if (!(lsame(TRANS, 'N') || lsame(TRANS, 'T'))) {
    INFO.value = -1;
  } else if (M < 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, M)) {
    INFO.value = -6;
  } else if (LDB < max(1, max(M, N))) {
    INFO.value = -8;
  } else if (LWORK < max(1, MN + max(MN, NRHS)) && !LQUERY) {
    INFO.value = -10;
  }

  // Figure out optimal block size and optimal workspace size

  if (INFO.value == 0 || INFO.value == -10) {
    TPSD = true;
    if (lsame(TRANS, 'N')) TPSD = false;

    NB = ilaenv(1, 'DGELST', ' ', M, N, -1, -1);

    MNNRHS = max(MN, NRHS);
    LWOPT = max(1, (MN + MNNRHS) * NB);
    WORK[1] = LWOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DGELST', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (min(M, min(N, NRHS)) == 0) {
    dlaset('Full', max(M, N), NRHS, ZERO, ZERO, B, LDB);
    WORK[1] = LWOPT.toDouble();
    return;
  }

  // *GEQRT and *GELQT routines cannot accept NB larger than min(M,N)

  if (NB > MN) NB = MN;

  // Determine the block size from the supplied LWORK
  // ( at this stage we know that LWORK >= (minimum required workspace,
  // but it may be less than optimal)

  NB = min(NB, LWORK ~/ (MN + MNNRHS));

  // The minimum value of NB, when blocked code is used

  NBMIN = max(2, ilaenv(2, 'DGELST', ' ', M, N, -1, -1));

  if (NB < NBMIN) {
    NB = 1;
  }

  // Get machine parameters

  SMLNUM = dlamch('S') / dlamch('P');
  BIGNUM = ONE / SMLNUM;

  // Scale A, B if max element outside range [SMLNUM,BIGNUM]

  ANRM = dlange('M', M, N, A, LDA, RWORK);
  IASCL = 0;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    dlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO);
    IASCL = 1;
  } else if (ANRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    dlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO);
    IASCL = 2;
  } else if (ANRM == ZERO) {
    // Matrix all zero. Return zero solution.

    dlaset('Full', max(M, N), NRHS, ZERO, ZERO, B, LDB);
    WORK[1] = LWOPT.toDouble();
    return;
  }

  BROW = M;
  if (TPSD) BROW = N;
  BNRM = dlange('M', BROW, NRHS, B, LDB, RWORK);
  IBSCL = 0;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    dlascl('G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO);
    IBSCL = 1;
  } else if (BNRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    dlascl('G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO);
    IBSCL = 2;
  }

  if (M >= N) {
    // M > N:
    // Compute the blocked QR factorization of A,
    // using the compact WY representation of Q,
    // workspace at least N, optimally N*NB.

    dgeqrt(M, N, NB, A, LDA, WORK(1).asMatrix(NB), NB, WORK(MN * NB + 1), INFO);

    if (!TPSD) {
      // M > N, A is not transposed:
      // Overdetermined system of equations,
      // least-squares problem, min || A * X - B ||.

      // Compute B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS),
      // using the compact WY representation of Q,
      // workspace at least NRHS, optimally NRHS*NB.

      dgemqrt('Left', 'Transpose', M, NRHS, N, NB, A, LDA, WORK(1).asMatrix(NB),
          NB, B, LDB, WORK(MN * NB + 1), INFO);

      // Compute B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

      dtrtrs(
          'Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      SCLLEN = N;
    } else {
      // M > N, A is transposed:
      // Underdetermined system of equations,
      // minimum norm solution of A**T * X = B.

      // Compute B := inv(R**T) * B in two row blocks of B.

      // Block 1: B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

      dtrtrs('Upper', 'Transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      // Block 2: Zero out all rows below the N-th row in B:
      // B(N+1:M,1:NRHS) = ZERO

      for (J = 1; J <= NRHS; J++) {
        for (I = N + 1; I <= M; I++) {
          B[I][J] = ZERO;
        }
      }

      // Compute B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS),
      // using the compact WY representation of Q,
      // workspace at least NRHS, optimally NRHS*NB.

      dgemqrt('Left', 'No transpose', M, NRHS, N, NB, A, LDA,
          WORK(1).asMatrix(NB), NB, B, LDB, WORK(MN * NB + 1), INFO);

      SCLLEN = M;
    }
  } else {
    // M < N:
    // Compute the blocked LQ factorization of A,
    // using the compact WY representation of Q,
    // workspace at least M, optimally M*NB.

    dgelqt(M, N, NB, A, LDA, WORK(1).asMatrix(NB), NB, WORK(MN * NB + 1), INFO);

    if (!TPSD) {
      // M < N, A is not transposed:
      // Underdetermined system of equations,
      // minimum norm solution of A * X = B.

      // Compute B := inv(L) * B in two row blocks of B.

      // Block 1: B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

      dtrtrs(
          'Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      // Block 2: Zero out all rows below the M-th row in B:
      // B(M+1:N,1:NRHS) = ZERO

      for (J = 1; J <= NRHS; J++) {
        for (I = M + 1; I <= N; I++) {
          B[I][J] = ZERO;
        }
      }

      // Compute B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS),
      // using the compact WY representation of Q,
      // workspace at least NRHS, optimally NRHS*NB.

      dgemlqt('Left', 'Transpose', N, NRHS, M, NB, A, LDA, WORK(1).asMatrix(NB),
          NB, B, LDB, WORK(MN * NB + 1), INFO);

      SCLLEN = N;
    } else {
      // M < N, A is transposed:
      // Overdetermined system of equations,
      // least-squares problem, min || A**T * X - B ||.

      // Compute B(1:N,1:NRHS) := Q * B(1:N,1:NRHS),
      // using the compact WY representation of Q,
      // workspace at least NRHS, optimally NRHS*NB.

      dgemlqt('Left', 'No transpose', N, NRHS, M, NB, A, LDA,
          WORK(1).asMatrix(NB), NB, B, LDB, WORK(MN * NB + 1), INFO);

      // Compute B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

      dtrtrs('Lower', 'Transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      SCLLEN = M;
    }
  }

  // Undo scaling

  if (IASCL == 1) {
    dlascl('G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO);
  } else if (IASCL == 2) {
    dlascl('G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO);
  }
  if (IBSCL == 1) {
    dlascl('G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO);
  } else if (IBSCL == 2) {
    dlascl('G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO);
  }

  WORK[1] = LWOPT.toDouble();
}
