import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgelqf.dart';
import 'package:lapack/src/dgeqrf.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dormlq.dart';
import 'package:lapack/src/dormqr.dart';
import 'package:lapack/src/dtrtrs.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgels(
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
  int BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE = 0;
  double ANRM, BIGNUM, BNRM, SMLNUM;
  final RWORK = Array<double>(1);

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

  // Figure out optimal block size

  if (INFO.value == 0 || INFO.value == -10) {
    TPSD = true;
    if (lsame(TRANS, 'N')) TPSD = false;

    if (M >= N) {
      NB = ilaenv(1, 'DGEQRF', ' ', M, N, -1, -1);
      if (TPSD) {
        NB = max(NB, ilaenv(1, 'DORMQR', 'LN', M, NRHS, N, -1));
      } else {
        NB = max(NB, ilaenv(1, 'DORMQR', 'LT', M, NRHS, N, -1));
      }
    } else {
      NB = ilaenv(1, 'DGELQF', ' ', M, N, -1, -1);
      if (TPSD) {
        NB = max(NB, ilaenv(1, 'DORMLQ', 'LT', N, NRHS, M, -1));
      } else {
        NB = max(NB, ilaenv(1, 'DORMLQ', 'LN', N, NRHS, M, -1));
      }
    }

    WSIZE = max(1, MN + max(MN, NRHS) * NB);
    WORK[1] = WSIZE.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DGELS ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (min(M, min(N, NRHS)) == 0) {
    dlaset('Full', max(M, N), NRHS, ZERO, ZERO, B, LDB);
    return;
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

    dlaset('F', max(M, N), NRHS, ZERO, ZERO, B, LDB);
    WORK[1] = WSIZE.toDouble();
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
    // compute QR factorization of A

    dgeqrf(M, N, A, LDA, WORK(1), WORK(MN + 1), LWORK - MN, INFO);

    // workspace at least N, optimally N*NB

    if (!TPSD) {
      // Least-Squares Problem min || A * X - B ||

      // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

      dormqr('Left', 'Transpose', M, NRHS, N, A, LDA, WORK(1), B, LDB,
          WORK(MN + 1), LWORK - MN, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

      dtrtrs(
          'Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      SCLLEN = N;
    } else {
      // Underdetermined system of equations A**T * X = B

      // B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

      dtrtrs('Upper', 'Transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      // B(N+1:M,1:NRHS) = ZERO

      for (J = 1; J <= NRHS; J++) {
        for (I = N + 1; I <= M; I++) {
          B[I][J] = ZERO;
        }
      }

      // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

      dormqr('Left', 'No transpose', M, NRHS, N, A, LDA, WORK(1), B, LDB,
          WORK(MN + 1), LWORK - MN, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      SCLLEN = M;
    }
  } else {
    // Compute LQ factorization of A

    dgelqf(M, N, A, LDA, WORK(1), WORK(MN + 1), LWORK - MN, INFO);

    // workspace at least M, optimally M*NB.

    if (!TPSD) {
      // underdetermined system of equations A * X = B

      // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

      dtrtrs(
          'Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      // B(M+1:N,1:NRHS) = 0

      for (J = 1; J <= NRHS; J++) {
        for (I = M + 1; I <= N; I++) {
          B[I][J] = ZERO;
        }
      }

      // B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)

      dormlq('Left', 'Transpose', N, NRHS, M, A, LDA, WORK(1), B, LDB,
          WORK(MN + 1), LWORK - MN, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      SCLLEN = N;
    } else {
      // overdetermined system min || A**T * X - B ||

      // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

      dormlq('Left', 'No transpose', N, NRHS, M, A, LDA, WORK(1), B, LDB,
          WORK(MN + 1), LWORK - MN, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      // B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

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

  WORK[1] = WSIZE.toDouble();
}
