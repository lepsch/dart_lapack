import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgelq.dart';
import 'package:lapack/src/dgemlq.dart';
import 'package:lapack/src/dgemqr.dart';
import 'package:lapack/src/dgeqr.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dtrtrs.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgetsls(
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
  bool LQUERY, TRAN;
  int I,
      IASCL,
      IBSCL,
      J,
      MAXMN,
      BROW,
      SCLLEN,
      TSZO = 0,
      TSZM = 0,
      LWO = 0,
      LWM = 0,
      LW1,
      LW2,
      WSIZEO = 0,
      WSIZEM = 0;
  double ANRM, BIGNUM, BNRM, SMLNUM;
  final TQ = Array<double>(5), WORKQ = Array<double>(1);
  final INFO2 = Box(0);

  // Test the input arguments.

  INFO.value = 0;
  MAXMN = max(M, N);
  TRAN = lsame(TRANS, 'T');

  LQUERY = (LWORK == -1 || LWORK == -2);
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
  }

  if (INFO.value == 0) {
    // Determine the optimum and minimum LWORK

    if (min(M, min(N, NRHS)) == 0) {
      WSIZEM = 1;
      WSIZEO = 1;
    } else if (M >= N) {
      dgeqr(M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2);
      TSZO = TQ[1].toInt();
      LWO = WORKQ[1].toInt();
      dgemqr(
          'L', TRANS, M, NRHS, N, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2);
      LWO = max(LWO, WORKQ[1]).toInt();
      dgeqr(M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2);
      TSZM = TQ[1].toInt();
      LWM = WORKQ[1].toInt();
      dgemqr(
          'L', TRANS, M, NRHS, N, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2);
      LWM = max(LWM, WORKQ[1]).toInt();
      WSIZEO = TSZO + LWO;
      WSIZEM = TSZM + LWM;
    } else {
      dgelq(M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2);
      TSZO = TQ[1].toInt();
      LWO = WORKQ[1].toInt();
      dgemlq(
          'L', TRANS, N, NRHS, M, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2);
      LWO = max(LWO, WORKQ[1]).toInt();
      dgelq(M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2);
      TSZM = TQ[1].toInt();
      LWM = WORKQ[1].toInt();
      dgemlq(
          'L', TRANS, N, NRHS, M, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2);
      LWM = max(LWM, WORKQ[1]).toInt();
      WSIZEO = TSZO + LWO;
      WSIZEM = TSZM + LWM;
    }

    if ((LWORK < WSIZEM) && !LQUERY) {
      INFO.value = -10;
    }

    WORK[1] = WSIZEO.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DGETSLS', -INFO.value);
    return;
  }
  if (LQUERY) {
    if (LWORK == -2) WORK[1] = WSIZEM.toDouble();
    return;
  }
  if (LWORK < WSIZEO) {
    LW1 = TSZM;
    LW2 = LWM;
  } else {
    LW1 = TSZO;
    LW2 = LWO;
  }

  // Quick return if possible

  if (min(M, min(N, NRHS)) == 0) {
    dlaset('FULL', max(M, N), NRHS, ZERO, ZERO, B, LDB);
    return;
  }

  // Get machine parameters

  SMLNUM = dlamch('S') / dlamch('P');
  BIGNUM = ONE / SMLNUM;

  // Scale A, B if max element outside range [SMLNUM,BIGNUM]

  ANRM = dlange('M', M, N, A, LDA, WORK);
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

    dlaset('F', MAXMN, NRHS, ZERO, ZERO, B, LDB);
    WORK[1] = (TSZO + LWO).toDouble();
    return;
  }

  BROW = M;
  if (TRAN) {
    BROW = N;
  }
  BNRM = dlange('M', BROW, NRHS, B, LDB, WORK);
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

    dgeqr(M, N, A, LDA, WORK(LW2 + 1), LW1, WORK(1), LW2, INFO);
    if (!TRAN) {
      // Least-Squares Problem min || A * X - B ||

      // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

      dgemqr('L', 'T', M, NRHS, N, A, LDA, WORK(LW2 + 1), LW1, B, LDB, WORK(1),
          LW2, INFO);

      // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

      dtrtrs('U', 'N', 'N', N, NRHS, A, LDA, B, LDB, INFO);
      if (INFO.value > 0) {
        return;
      }
      SCLLEN = N;
    } else {
      // Overdetermined system of equations A**T * X = B

      // B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

      dtrtrs('U', 'T', 'N', N, NRHS, A, LDA, B, LDB, INFO);

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

      dgemqr('L', 'N', M, NRHS, N, A, LDA, WORK(LW2 + 1), LW1, B, LDB, WORK(1),
          LW2, INFO);

      SCLLEN = M;
    }
  } else {
    // Compute LQ factorization of A

    dgelq(M, N, A, LDA, WORK(LW2 + 1), LW1, WORK(1), LW2, INFO);

    // workspace at least M, optimally M*NB.

    if (!TRAN) {
      // underdetermined system of equations A * X = B

      // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

      dtrtrs('L', 'N', 'N', M, NRHS, A, LDA, B, LDB, INFO);

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

      dgemlq('L', 'T', N, NRHS, M, A, LDA, WORK(LW2 + 1), LW1, B, LDB, WORK(1),
          LW2, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      SCLLEN = N;
    } else {
      // overdetermined system min || A**T * X - B ||

      // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

      dgemlq('L', 'N', N, NRHS, M, A, LDA, WORK(LW2 + 1), LW1, B, LDB, WORK(1),
          LW2, INFO);

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
  WORK[1] = (TSZO + LWO).toDouble();
}
