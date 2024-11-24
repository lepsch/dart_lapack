// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dbdsvdx.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgeqrf.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgebrd.dart';
import 'package:dart_lapack/src/zgelqf.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zunmbr.dart';
import 'package:dart_lapack/src/zunmlq.dart';
import 'package:dart_lapack/src/zunmqr.dart';

void zgesvdx(
  final String JOBU,
  final String JOBVT,
  final String RANGE,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final Box<int> NS,
  final Array<double> S_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> VT_,
  final int LDVT,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final WORK = WORK_.having();
  final S = S_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  String JOBZ, RNGTGK;
  bool ALLS, INDS, LQUERY, VALS, WANTU, WANTVT;
  int I,
      ID,
      IE,
      ILQF,
      ILTGK,
      IQRF,
      ISCL,
      ITAU,
      ITAUP,
      ITAUQ,
      ITEMP,
      ITEMPR,
      ITGKZ,
      IUTGK,
      J,
      K,
      MAXWRK = 0,
      MINMN,
      MINWRK,
      MNTHR = 0;
  double ANRM, BIGNUM, EPS, SMLNUM;
  final DUM = Array<double>(1);
  final IERR = Box(0);

  // Test the input arguments.

  NS.value = 0;
  INFO.value = 0;
  LQUERY = (LWORK == -1);
  MINMN = min(M, N);

  WANTU = lsame(JOBU, 'V');
  WANTVT = lsame(JOBVT, 'V');
  if (WANTU || WANTVT) {
    JOBZ = 'V';
  } else {
    JOBZ = 'N';
  }
  ALLS = lsame(RANGE, 'A');
  VALS = lsame(RANGE, 'V');
  INDS = lsame(RANGE, 'I');

  INFO.value = 0;
  if (!lsame(JOBU, 'V') && !lsame(JOBU, 'N')) {
    INFO.value = -1;
  } else if (!lsame(JOBVT, 'V') && !lsame(JOBVT, 'N')) {
    INFO.value = -2;
  } else if (!(ALLS || VALS || INDS)) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (M > LDA) {
    INFO.value = -7;
  } else if (MINMN > 0) {
    if (VALS) {
      if (VL < ZERO) {
        INFO.value = -8;
      } else if (VU <= VL) {
        INFO.value = -9;
      }
    } else if (INDS) {
      if (IL < 1 || IL > max(1, MINMN)) {
        INFO.value = -10;
      } else if (IU < min(MINMN, IL) || IU > MINMN) {
        INFO.value = -11;
      }
    }
    if (INFO.value == 0) {
      if (WANTU && LDU < M) {
        INFO.value = -15;
      } else if (WANTVT) {
        if (INDS) {
          if (LDVT < IU - IL + 1) {
            INFO.value = -17;
          }
        } else if (LDVT < MINMN) {
          INFO.value = -17;
        }
      }
    }
  }

  // Compute workspace
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ILAENV.)

  if (INFO.value == 0) {
    MINWRK = 1;
    MAXWRK = 1;
    if (MINMN > 0) {
      if (M >= N) {
        MNTHR = ilaenv(6, 'ZGESVD', JOBU + JOBVT, M, N, 0, 0);
        if (M >= MNTHR) {
          // Path 1 (M much larger than N)

          MINWRK = N * (N + 5);
          MAXWRK = N + N * ilaenv(1, 'ZGEQRF', ' ', M, N, -1, -1);
          MAXWRK = max(MAXWRK,
              N * N + 2 * N + 2 * N * ilaenv(1, 'ZGEBRD', ' ', N, N, -1, -1));
          if (WANTU || WANTVT) {
            MAXWRK = max(MAXWRK,
                N * N + 2 * N + N * ilaenv(1, 'ZUNMQR', 'LN', N, N, N, -1));
          }
        } else {
          // Path 2 (M at least N, but not much larger)

          MINWRK = 3 * N + M;
          MAXWRK = 2 * N + (M + N) * ilaenv(1, 'ZGEBRD', ' ', M, N, -1, -1);
          if (WANTU || WANTVT) {
            MAXWRK =
                max(MAXWRK, 2 * N + N * ilaenv(1, 'ZUNMQR', 'LN', N, N, N, -1));
          }
        }
      } else {
        MNTHR = ilaenv(6, 'ZGESVD', JOBU + JOBVT, M, N, 0, 0);
        if (N >= MNTHR) {
          // Path 1t (N much larger than M)

          MINWRK = M * (M + 5);
          MAXWRK = M + M * ilaenv(1, 'ZGELQF', ' ', M, N, -1, -1);
          MAXWRK = max(MAXWRK,
              M * M + 2 * M + 2 * M * ilaenv(1, 'ZGEBRD', ' ', M, M, -1, -1));
          if (WANTU || WANTVT) {
            MAXWRK = max(MAXWRK,
                M * M + 2 * M + M * ilaenv(1, 'ZUNMQR', 'LN', M, M, M, -1));
          }
        } else {
          // Path 2t (N greater than M, but not much larger)

          MINWRK = 3 * M + N;
          MAXWRK = 2 * M + (M + N) * ilaenv(1, 'ZGEBRD', ' ', M, N, -1, -1);
          if (WANTU || WANTVT) {
            MAXWRK =
                max(MAXWRK, 2 * M + M * ilaenv(1, 'ZUNMQR', 'LN', M, M, M, -1));
          }
        }
      }
    }
    MAXWRK = max(MAXWRK, MINWRK);
    WORK[1] = MAXWRK.toComplex();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -19;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZGESVDX', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    return;
  }

  // Set singular values indices accord to RANGE='A'.

  if (ALLS) {
    RNGTGK = 'I';
    ILTGK = 1;
    IUTGK = min(M, N);
  } else if (INDS) {
    RNGTGK = 'I';
    ILTGK = IL;
    IUTGK = IU;
  } else {
    RNGTGK = 'V';
    ILTGK = 0;
    IUTGK = 0;
  }

  // Get machine constants

  EPS = dlamch('P');
  SMLNUM = sqrt(dlamch('S')) / EPS;
  BIGNUM = ONE / SMLNUM;

  // Scale A if max element outside range [SMLNUM,BIGNUM]

  ANRM = zlange('M', M, N, A, LDA, DUM);
  ISCL = 0;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    ISCL = 1;
    zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO);
  } else if (ANRM > BIGNUM) {
    ISCL = 1;
    zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO);
  }

  if (M >= N) {
    // A has at least as many rows as columns. If A has sufficiently
    // more rows than columns, first reduce A using the QR
    // decomposition.

    if (M >= MNTHR) {
      // Path 1 (M much larger than N):
      // A = Q * R = Q * ( QB * B * PB**T )
      //           = Q * ( QB * ( UB * S * VB**T ) * PB**T )
      // U = Q * QB * UB; V**T = VB**T * PB**T

      // Compute A=Q*R
      // (Workspace: need 2*N, prefer N+N*NB)

      ITAU = 1;
      ITEMP = ITAU + N;
      zgeqrf(M, N, A, LDA, WORK(ITAU), WORK(ITEMP), LWORK - ITEMP + 1, INFO);

      // Copy R into WORK and bidiagonalize it:
      // (Workspace: need N*N+3*N, prefer N*N+N+2*N*NB)

      IQRF = ITEMP;
      ITAUQ = ITEMP + N * N;
      ITAUP = ITAUQ + N;
      ITEMP = ITAUP + N;
      ID = 1;
      IE = ID + N;
      ITGKZ = IE + N;
      zlacpy('U', N, N, A, LDA, WORK(IQRF).asMatrix(N), N);
      zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
          WORK(IQRF + 1).asMatrix(N), N);
      zgebrd(N, N, WORK(IQRF).asMatrix(N), N, RWORK(ID), RWORK(IE), WORK(ITAUQ),
          WORK(ITAUP), WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      ITEMPR = ITGKZ + N * (N * 2 + 1);

      // Solve eigenvalue problem TGK*Z *=S.
      // (Workspace: need 2*N*N+14*N)

      dbdsvdx(
          'U',
          JOBZ,
          RNGTGK,
          N,
          RWORK(ID),
          RWORK(IE),
          VL,
          VU,
          ILTGK,
          IUTGK,
          NS,
          S,
          RWORK(ITGKZ).asMatrix(N * 2),
          N * 2,
          RWORK(ITEMPR),
          IWORK,
          INFO);

      // If needed, compute left singular vectors.

      if (WANTU) {
        K = ITGKZ;
        for (I = 1; I <= NS.value; I++) {
          for (J = 1; J <= N; J++) {
            U[J][I] = Complex(RWORK[K], ZERO);
            K++;
          }
          K += N;
        }
        zlaset(
            'A', M - N, NS.value, Complex.zero, Complex.zero, U(N + 1, 1), LDU);

        // Call ZUNMBR to compute QB*UB.
        // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

        zunmbr('Q', 'L', 'N', N, NS.value, N, WORK(IQRF).asMatrix(N), N,
            WORK(ITAUQ), U, LDU, WORK(ITEMP), LWORK - ITEMP + 1, INFO);

        // Call ZUNMQR to compute Q*(QB*UB).
        // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

        zunmqr('L', 'N', M, NS.value, N, A, LDA, WORK(ITAU), U, LDU,
            WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      }

      // If needed, compute right singular vectors.

      if (WANTVT) {
        K = ITGKZ + N;
        for (I = 1; I <= NS.value; I++) {
          for (J = 1; J <= N; J++) {
            VT[I][J] = Complex(RWORK[K], ZERO);
            K++;
          }
          K += N;
        }

        // Call ZUNMBR to compute VB**T * PB**T
        // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

        zunmbr('P', 'R', 'C', NS.value, N, N, WORK(IQRF).asMatrix(N), N,
            WORK(ITAUP), VT, LDVT, WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      }
    } else {
      // Path 2 (M at least N, but not much larger)
      // Reduce A to bidiagonal form without QR decomposition
      // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
      // U = QB * UB; V**T = VB**T * PB**T

      // Bidiagonalize A
      // (Workspace: need 2*N+M, prefer 2*N+(M+N)*NB)

      ITAUQ = 1;
      ITAUP = ITAUQ + N;
      ITEMP = ITAUP + N;
      ID = 1;
      IE = ID + N;
      ITGKZ = IE + N;
      zgebrd(M, N, A, LDA, RWORK(ID), RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
          WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      ITEMPR = ITGKZ + N * (N * 2 + 1);

      // Solve eigenvalue problem TGK*Z *=S.
      // (Workspace: need 2*N*N+14*N)

      dbdsvdx(
          'U',
          JOBZ,
          RNGTGK,
          N,
          RWORK(ID),
          RWORK(IE),
          VL,
          VU,
          ILTGK,
          IUTGK,
          NS,
          S,
          RWORK(ITGKZ).asMatrix(N * 2),
          N * 2,
          RWORK(ITEMPR),
          IWORK,
          INFO);

      // If needed, compute left singular vectors.

      if (WANTU) {
        K = ITGKZ;
        for (I = 1; I <= NS.value; I++) {
          for (J = 1; J <= N; J++) {
            U[J][I] = Complex(RWORK[K], ZERO);
            K++;
          }
          K += N;
        }
        zlaset(
            'A', M - N, NS.value, Complex.zero, Complex.zero, U(N + 1, 1), LDU);

        // Call ZUNMBR to compute QB*UB.
        // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

        zunmbr('Q', 'L', 'N', M, NS.value, N, A, LDA, WORK(ITAUQ), U, LDU,
            WORK(ITEMP), LWORK - ITEMP + 1, IERR);
      }

      // If needed, compute right singular vectors.

      if (WANTVT) {
        K = ITGKZ + N;
        for (I = 1; I <= NS.value; I++) {
          for (J = 1; J <= N; J++) {
            VT[I][J] = Complex(RWORK[K], ZERO);
            K++;
          }
          K += N;
        }

        // Call ZUNMBR to compute VB**T * PB**T
        // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

        zunmbr('P', 'R', 'C', NS.value, N, N, A, LDA, WORK(ITAUP), VT, LDVT,
            WORK(ITEMP), LWORK - ITEMP + 1, IERR);
      }
    }
  } else {
    // A has more columns than rows. If A has sufficiently more
    // columns than rows, first reduce A using the LQ decomposition.

    if (N >= MNTHR) {
      // Path 1t (N much larger than M):
      // A = L * Q = ( QB * B * PB**T ) * Q
      //           = ( QB * ( UB * S * VB**T ) * PB**T ) * Q
      // U = QB * UB ; V**T = VB**T * PB**T * Q

      // Compute A=L*Q
      // (Workspace: need 2*M, prefer M+M*NB)

      ITAU = 1;
      ITEMP = ITAU + M;
      zgelqf(M, N, A, LDA, WORK(ITAU), WORK(ITEMP), LWORK - ITEMP + 1, INFO);

      // Copy L into WORK and bidiagonalize it:
      // (Workspace in WORK( ITEMP ): need M*M+3*M, prefer M*M+M+2*M*NB)

      ILQF = ITEMP;
      ITAUQ = ILQF + M * M;
      ITAUP = ITAUQ + M;
      ITEMP = ITAUP + M;
      ID = 1;
      IE = ID + M;
      ITGKZ = IE + M;
      zlacpy('L', M, M, A, LDA, WORK(ILQF).asMatrix(M), M);
      zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
          WORK(ILQF + M).asMatrix(M), M);
      zgebrd(M, M, WORK(ILQF).asMatrix(M), M, RWORK(ID), RWORK(IE), WORK(ITAUQ),
          WORK(ITAUP), WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      ITEMPR = ITGKZ + M * (M * 2 + 1);

      // Solve eigenvalue problem TGK*Z *=S.
      // (Workspace: need 2*M*M+14*M)

      dbdsvdx(
          'U',
          JOBZ,
          RNGTGK,
          M,
          RWORK(ID),
          RWORK(IE),
          VL,
          VU,
          ILTGK,
          IUTGK,
          NS,
          S,
          RWORK(ITGKZ).asMatrix(M * 2),
          M * 2,
          RWORK(ITEMPR),
          IWORK,
          INFO);

      // If needed, compute left singular vectors.

      if (WANTU) {
        K = ITGKZ;
        for (I = 1; I <= NS.value; I++) {
          for (J = 1; J <= M; J++) {
            U[J][I] = Complex(RWORK[K], ZERO);
            K++;
          }
          K += M;
        }

        // Call ZUNMBR to compute QB*UB.
        // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

        zunmbr('Q', 'L', 'N', M, NS.value, M, WORK(ILQF).asMatrix(M), M,
            WORK(ITAUQ), U, LDU, WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      }

      // If needed, compute right singular vectors.

      if (WANTVT) {
        K = ITGKZ + M;
        for (I = 1; I <= NS.value; I++) {
          for (J = 1; J <= M; J++) {
            VT[I][J] = Complex(RWORK[K], ZERO);
            K++;
          }
          K += M;
        }
        zlaset('A', NS.value, N - M, Complex.zero, Complex.zero, VT(1, M + 1),
            LDVT);

        // Call ZUNMBR to compute (VB**T)*(PB**T)
        // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

        zunmbr('P', 'R', 'C', NS.value, M, M, WORK(ILQF).asMatrix(M), M,
            WORK(ITAUP), VT, LDVT, WORK(ITEMP), LWORK - ITEMP + 1, INFO);

        // Call ZUNMLQ to compute ((VB**T)*(PB**T))*Q.
        // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

        zunmlq('R', 'N', NS.value, N, M, A, LDA, WORK(ITAU), VT, LDVT,
            WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      }
    } else {
      // Path 2t (N greater than M, but not much larger)
      // Reduce to bidiagonal form without LQ decomposition
      // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
      // U = QB * UB; V**T = VB**T * PB**T

      // Bidiagonalize A
      // (Workspace: need 2*M+N, prefer 2*M+(M+N)*NB)

      ITAUQ = 1;
      ITAUP = ITAUQ + M;
      ITEMP = ITAUP + M;
      ID = 1;
      IE = ID + M;
      ITGKZ = IE + M;
      zgebrd(M, N, A, LDA, RWORK(ID), RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
          WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      ITEMPR = ITGKZ + M * (M * 2 + 1);

      // Solve eigenvalue problem TGK*Z *=S.
      // (Workspace: need 2*M*M+14*M)

      dbdsvdx(
          'L',
          JOBZ,
          RNGTGK,
          M,
          RWORK(ID),
          RWORK(IE),
          VL,
          VU,
          ILTGK,
          IUTGK,
          NS,
          S,
          RWORK(ITGKZ).asMatrix(M * 2),
          M * 2,
          RWORK(ITEMPR),
          IWORK,
          INFO);

      // If needed, compute left singular vectors.

      if (WANTU) {
        K = ITGKZ;
        for (I = 1; I <= NS.value; I++) {
          for (J = 1; J <= M; J++) {
            U[J][I] = Complex(RWORK[K], ZERO);
            K++;
          }
          K += M;
        }

        // Call ZUNMBR to compute QB*UB.
        // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

        zunmbr('Q', 'L', 'N', M, NS.value, N, A, LDA, WORK(ITAUQ), U, LDU,
            WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      }

      // If needed, compute right singular vectors.

      if (WANTVT) {
        K = ITGKZ + M;
        for (I = 1; I <= NS.value; I++) {
          for (J = 1; J <= M; J++) {
            VT[I][J] = Complex(RWORK[K], ZERO);
            K++;
          }
          K += M;
        }
        zlaset('A', NS.value, N - M, Complex.zero, Complex.zero, VT(1, M + 1),
            LDVT);

        // Call ZUNMBR to compute VB**T * PB**T
        // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

        zunmbr('P', 'R', 'C', NS.value, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
            WORK(ITEMP), LWORK - ITEMP + 1, INFO);
      }
    }
  }

  // Undo scaling if necessary

  if (ISCL == 1) {
    if (ANRM > BIGNUM) {
      dlascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, INFO);
    }
    if (ANRM < SMLNUM) {
      dlascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, INFO);
    }
  }

  // Return optimal workspace in WORK(1)

  WORK[1] = MAXWRK.toComplex();
}
