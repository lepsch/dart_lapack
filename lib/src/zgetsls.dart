import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgelq.dart';
import 'package:lapack/src/zgemlq.dart';
import 'package:lapack/src/zgemqr.dart';
import 'package:lapack/src/zgeqr.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztrtrs.dart';

void zgetsls(
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
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final WORK = WORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
  final INFO2 = Box(0);
  final DUM = Array<double>(1);
  final TQ = Array<Complex>(5), WORKQ = Array<Complex>(1);

  // Test the input arguments.

  INFO.value = 0;
  MAXMN = max(M, N);
  TRAN = lsame(TRANS, 'C');

  LQUERY = (LWORK == -1 || LWORK == -2);
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
  }

  if (INFO.value == 0) {
    // Determine the optimum and minimum LWORK

    if (min(M, min(N, NRHS)) == 0) {
      WSIZEO = 1;
      WSIZEM = 1;
    } else if (M >= N) {
      zgeqr(M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2);
      TSZO = TQ[1].toInt();
      LWO = WORKQ[1].toInt();
      zgemqr(
          'L', TRANS, M, NRHS, N, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2);
      LWO = max(LWO, WORKQ[1].toInt());
      zgeqr(M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2);
      TSZM = TQ[1].toInt();
      LWM = WORKQ[1].toInt();
      zgemqr(
          'L', TRANS, M, NRHS, N, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2);
      LWM = max(LWM, WORKQ[1].toInt());
      WSIZEO = TSZO + LWO;
      WSIZEM = TSZM + LWM;
    } else {
      zgelq(M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2);
      TSZO = TQ[1].toInt();
      LWO = WORKQ[1].toInt();
      zgemlq(
          'L', TRANS, N, NRHS, M, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2);
      LWO = max(LWO, WORKQ[1].toInt());
      zgelq(M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2);
      TSZM = TQ[1].toInt();
      LWM = WORKQ[1].toInt();
      zgemlq(
          'L', TRANS, N, NRHS, M, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2);
      LWM = max(LWM, WORKQ[1].toInt());
      WSIZEO = TSZO + LWO;
      WSIZEM = TSZM + LWM;
    }

    if ((LWORK < WSIZEM) && (!LQUERY)) {
      INFO.value = -10;
    }

    WORK[1] = WSIZEO.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZGETSLS', -INFO.value);
    return;
  }
  if (LQUERY) {
    if (LWORK == -2) WORK[1] = WSIZEM.toComplex();
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
    zlaset('FULL', max(M, N), NRHS, Complex.zero, Complex.zero, B, LDB);
    return;
  }

  // Get machine parameters

  SMLNUM = dlamch('S') / dlamch('P');
  BIGNUM = ONE / SMLNUM;

  // Scale A, B if max element outside range [SMLNUM,BIGNUM]

  ANRM = zlange('M', M, N, A, LDA, DUM);
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

    zlaset('F', MAXMN, NRHS, Complex.zero, Complex.zero, B, LDB);
    WORK[1] = (TSZO + LWO).toComplex();
    return;
  }

  BROW = M;
  if (TRAN) {
    BROW = N;
  }
  BNRM = zlange('M', BROW, NRHS, B, LDB, DUM);
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

    zgeqr(M, N, A, LDA, WORK(LW2 + 1), LW1, WORK(1), LW2, INFO);
    if (!TRAN) {
      // Least-Squares Problem min || A * X - B ||

      // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

      zgemqr('L', 'C', M, NRHS, N, A, LDA, WORK(LW2 + 1), LW1, B, LDB, WORK(1),
          LW2, INFO);

      // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

      ztrtrs('U', 'N', 'N', N, NRHS, A, LDA, B, LDB, INFO);
      if (INFO.value > 0) {
        return;
      }
      SCLLEN = N;
    } else {
      // Overdetermined system of equations A**T * X = B

      // B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

      ztrtrs('U', 'C', 'N', N, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      // B(N+1:M,1:NRHS) = Complex.zero

      for (J = 1; J <= NRHS; J++) {
        // 20
        for (I = N + 1; I <= M; I++) {
          // 10
          B[I][J] = Complex.zero;
        } // 10
      } // 20

      // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

      zgemqr('L', 'N', M, NRHS, N, A, LDA, WORK(LW2 + 1), LW1, B, LDB, WORK(1),
          LW2, INFO);

      SCLLEN = M;
    }
  } else {
    // Compute LQ factorization of A

    zgelq(M, N, A, LDA, WORK(LW2 + 1), LW1, WORK(1), LW2, INFO);

    // workspace at least M, optimally M*NB.

    if (!TRAN) {
      // underdetermined system of equations A * X = B

      // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

      ztrtrs('L', 'N', 'N', M, NRHS, A, LDA, B, LDB, INFO);

      if (INFO.value > 0) {
        return;
      }

      // B(M+1:N,1:NRHS) = 0

      for (J = 1; J <= NRHS; J++) {
        // 40
        for (I = M + 1; I <= N; I++) {
          // 30
          B[I][J] = Complex.zero;
        } // 30
      } // 40

      // B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)

      zgemlq('L', 'C', N, NRHS, M, A, LDA, WORK(LW2 + 1), LW1, B, LDB, WORK(1),
          LW2, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      SCLLEN = N;
    } else {
      // overdetermined system min || A**T * X - B ||

      // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

      zgemlq('L', 'N', N, NRHS, M, A, LDA, WORK(LW2 + 1), LW1, B, LDB, WORK(1),
          LW2, INFO);

      // workspace at least NRHS, optimally NRHS*NB

      // B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

      ztrtrs('L', 'C', 'N', M, NRHS, A, LDA, B, LDB, INFO);

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

  WORK[1] = (TSZO + LWO).toComplex();
}
