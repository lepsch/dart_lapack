// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgeqrf.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgelqf.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/ztrtrs.dart';
import 'package:dart_lapack/src/zunmlq.dart';
import 'package:dart_lapack/src/zunmqr.dart';

void zgels(
  final String TRANS,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
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
  int BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE = 0;
  double ANRM, BIGNUM, BNRM, SMLNUM;
  final RWORK = Array<double>(1);

  // Test the input arguments.

  INFO.value = 0;
  MN = min(M, N);
  LQUERY = (LWORK == -1);
  if (!(lsame(TRANS, 'N') || lsame(TRANS, 'C'))) {
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

  // Figure out optimal block size

  if (INFO.value == 0 || INFO.value == -10) {
    TPSD = true;
    if (lsame(TRANS, 'N')) TPSD = false;

    if (M >= N) {
      NB = ilaenv(1, 'ZGEQRF', ' ', M, N, -1, -1);
      if (TPSD) {
        NB = max(NB, ilaenv(1, 'ZUNMQR', 'LN', M, NRHS, N, -1));
      } else {
        NB = max(NB, ilaenv(1, 'ZUNMQR', 'LC', M, NRHS, N, -1));
      }
    } else {
      NB = ilaenv(1, 'ZGELQF', ' ', M, N, -1, -1);
      if (TPSD) {
        NB = max(NB, ilaenv(1, 'ZUNMLQ', 'LC', N, NRHS, M, -1));
      } else {
        NB = max(NB, ilaenv(1, 'ZUNMLQ', 'LN', N, NRHS, M, -1));
      }
    }

    WSIZE = max(1, MN + max(MN, NRHS) * NB);
    WORK[1] = WSIZE.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZGELS', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (min(M, min(N, NRHS)) == 0) {
    zlaset('Full', max(M, N), NRHS, Complex.zero, Complex.zero, B, LDB);
    return;
  }

  // Get machine parameters

  SMLNUM = dlamch('S') / dlamch('P');
  BIGNUM = ONE / SMLNUM;

  // Scale A, B if max element outside range [SMLNUM,BIGNUM]

  ANRM = zlange('M', M, N, A, LDA, RWORK);
  IASCL = 0;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO);
    IASCL = 1;
  } else if (ANRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO);
    IASCL = 2;
  } else if (ANRM == ZERO) {
    // Matrix all zero. Return zero solution.

    zlaset('F', max(M, N), NRHS, Complex.zero, Complex.zero, B, LDB);
    WORK[1] = WSIZE.toComplex();
    return;
  }

  BROW = M;
  if (TPSD) BROW = N;
  BNRM = zlange('M', BROW, NRHS, B, LDB, RWORK);
  IBSCL = 0;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    zlascl('G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO);
    IBSCL = 1;
  } else if (BNRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    zlascl('G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO);
    IBSCL = 2;
  }

  if (M >= N) {
    // compute QR factorization of A

    zgeqrf(M, N, A, LDA, WORK(1), WORK(MN + 1), LWORK - MN, INFO);

    // workspace at least N, optimally N*NB

    if (!TPSD) {
      // Least-Squares Problem min || A * X - B ||

      // B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)

      zunmqr('Left', 'Conjugate transpose', M, NRHS, N, A, LDA, WORK(1), B, LDB,
          WORK(MN + 1), LWORK - MN, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

      ztrtrs(
          'Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      SCLLEN = N;
    } else {
      // Underdetermined system of equations A**T * X = B

      // B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS)

      ztrtrs('Upper', 'Conjugate transpose', 'Non-unit', N, NRHS, A, LDA, B,
          LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      // B(N+1:M,1:NRHS) = ZERO

      for (J = 1; J <= NRHS; J++) {
        for (I = N + 1; I <= M; I++) {
          B[I][J] = Complex.zero;
        }
      }

      // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

      zunmqr('Left', 'No transpose', M, NRHS, N, A, LDA, WORK(1), B, LDB,
          WORK(MN + 1), LWORK - MN, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      SCLLEN = M;
    }
  } else {
    // Compute LQ factorization of A

    zgelqf(M, N, A, LDA, WORK(1), WORK(MN + 1), LWORK - MN, INFO);

    // workspace at least M, optimally M*NB.

    if (!TPSD) {
      // underdetermined system of equations A * X = B

      // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

      ztrtrs(
          'Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      // B(M+1:N,1:NRHS) = 0

      for (J = 1; J <= NRHS; J++) {
        for (I = M + 1; I <= N; I++) {
          B[I][J] = Complex.zero;
        }
      }

      // B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS)

      zunmlq('Left', 'Conjugate transpose', N, NRHS, M, A, LDA, WORK(1), B, LDB,
          WORK(MN + 1), LWORK - MN, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      SCLLEN = N;
    } else {
      // overdetermined system min || A**H * X - B ||

      // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

      zunmlq('Left', 'No transpose', N, NRHS, M, A, LDA, WORK(1), B, LDB,
          WORK(MN + 1), LWORK - MN, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      // B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS)

      ztrtrs('Lower', 'Conjugate transpose', 'Non-unit', M, NRHS, A, LDA, B,
          LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      SCLLEN = M;
    }
  }

  // Undo scaling

  if (IASCL == 1) {
    zlascl('G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO);
  } else if (IASCL == 2) {
    zlascl('G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO);
  }
  if (IBSCL == 1) {
    zlascl('G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO);
  } else if (IBSCL == 2) {
    zlascl('G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO);
  }

  WORK[1] = WSIZE.toComplex();
}
