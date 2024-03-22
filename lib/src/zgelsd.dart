import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/qr/ll/zgeqrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgebrd.dart';
import 'package:lapack/src/zgelqf.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlalsd.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zunmbr.dart';
import 'package:lapack/src/zunmlq.dart';
import 'package:lapack/src/zunmqr.dart';

void zgelsd(
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> S_,
  final double RCOND,
  final Box<int> RANK,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final S = S_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  bool LQUERY;
  int IASCL,
      IBSCL,
      IE,
      IL,
      ITAU,
      ITAUP,
      ITAUQ,
      LDWORK,
      LIWORK = 0,
      LRWORK = 0,
      MAXMN,
      MAXWRK = 0,
      MINMN,
      MINWRK,
      MM,
      MNTHR = 0,
      NLVL,
      NRWORK,
      NWORK,
      SMLSIZ = 0;
  double ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM;

  // Test the input arguments.

  INFO.value = 0;
  MINMN = min(M, N);
  MAXMN = max(M, N);
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDB < max(1, MAXMN)) {
    INFO.value = -7;
  }

  // Compute workspace.
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ILAENV.)

  if (INFO.value == 0) {
    MINWRK = 1;
    MAXWRK = 1;
    LIWORK = 1;
    LRWORK = 1;
    if (MINMN > 0) {
      SMLSIZ = ilaenv(9, 'ZGELSD', ' ', 0, 0, 0, 0);
      MNTHR = ilaenv(6, 'ZGELSD', ' ', M, N, NRHS, -1);
      NLVL = max((log(MINMN / (SMLSIZ + 1)) ~/ log(TWO)) + 1, 0);
      LIWORK = 3 * MINMN * NLVL + 11 * MINMN;
      MM = M;
      if (M >= N && M >= MNTHR) {
        // Path 1a - overdetermined, with many more rows than
        //           columns.

        MM = N;
        MAXWRK = max(MAXWRK, N * ilaenv(1, 'ZGEQRF', ' ', M, N, -1, -1));
        MAXWRK = max(MAXWRK, NRHS * ilaenv(1, 'ZUNMQR', 'LC', M, NRHS, N, -1));
      }
      if (M >= N) {
        // Path 1 - overdetermined or exactly determined.

        LRWORK = 10 * N +
            2 * N * SMLSIZ +
            8 * N * NLVL +
            3 * SMLSIZ * NRHS +
            max(pow(SMLSIZ + 1, 2).toInt(), N * (1 + NRHS) + 2 * NRHS);
        MAXWRK = max(
            MAXWRK, 2 * N + (MM + N) * ilaenv(1, 'ZGEBRD', ' ', MM, N, -1, -1));
        MAXWRK = max(
            MAXWRK, 2 * N + NRHS * ilaenv(1, 'ZUNMBR', 'QLC', MM, NRHS, N, -1));
        MAXWRK = max(MAXWRK,
            2 * N + (N - 1) * ilaenv(1, 'ZUNMBR', 'PLN', N, NRHS, N, -1));
        MAXWRK = max(MAXWRK, 2 * N + N * NRHS);
        MINWRK = max(2 * N + MM, 2 * N + N * NRHS);
      }
      if (N > M) {
        LRWORK = 10 * M +
            2 * M * SMLSIZ +
            8 * M * NLVL +
            3 * SMLSIZ * NRHS +
            max(pow(SMLSIZ + 1, 2).toInt(), N * (1 + NRHS) + 2 * NRHS);
        if (N >= MNTHR) {
          // Path 2a - underdetermined, with many more columns
          //           than rows.

          MAXWRK = M + M * ilaenv(1, 'ZGELQF', ' ', M, N, -1, -1);
          MAXWRK = max(MAXWRK,
              M * M + 4 * M + 2 * M * ilaenv(1, 'ZGEBRD', ' ', M, M, -1, -1));
          MAXWRK = max(
              MAXWRK,
              M * M +
                  4 * M +
                  NRHS * ilaenv(1, 'ZUNMBR', 'QLC', M, NRHS, M, -1));
          MAXWRK = max(
              MAXWRK,
              M * M +
                  4 * M +
                  (M - 1) * ilaenv(1, 'ZUNMLQ', 'LC', N, NRHS, M, -1));
          if (NRHS > 1) {
            MAXWRK = max(MAXWRK, M * M + M + M * NRHS);
          } else {
            MAXWRK = max(MAXWRK, M * M + 2 * M);
          }
          MAXWRK = max(MAXWRK, M * M + 4 * M + M * NRHS);
          // XXX: Ensure the Path 2a case below is triggered.  The workspace
          // calculation should use queries for all routines eventually.
          MAXWRK = max(MAXWRK,
              4 * M + M * M + max(max(M, 2 * M - 4), max(NRHS, N - 3 * M)));
        } else {
          // Path 2 - underdetermined.

          MAXWRK = 2 * M + (N + M) * ilaenv(1, 'ZGEBRD', ' ', M, N, -1, -1);
          MAXWRK = max(MAXWRK,
              2 * M + NRHS * ilaenv(1, 'ZUNMBR', 'QLC', M, NRHS, M, -1));
          MAXWRK = max(
              MAXWRK, 2 * M + M * ilaenv(1, 'ZUNMBR', 'PLN', N, NRHS, M, -1));
          MAXWRK = max(MAXWRK, 2 * M + M * NRHS);
        }
        MINWRK = max(2 * M + N, 2 * M + M * NRHS);
      }
    }
    MINWRK = min(MINWRK, MAXWRK);
    WORK[1] = MAXWRK.toComplex();
    IWORK[1] = LIWORK;
    RWORK[1] = LRWORK.toDouble();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -12;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZGELSD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible.

  if (M == 0 || N == 0) {
    RANK.value = 0;
    return;
  }

  // Get machine parameters.

  EPS = dlamch('P');
  SFMIN = dlamch('S');
  SMLNUM = SFMIN / EPS;
  BIGNUM = ONE / SMLNUM;

  // Scale A if max entry outside range [SMLNUM,BIGNUM].

  ANRM = zlange('M', M, N, A, LDA, RWORK);
  IASCL = 0;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO);
    IASCL = 1;
  } else if (ANRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM.

    zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO);
    IASCL = 2;
  } else if (ANRM == ZERO) {
    // Matrix all zero. Return zero solution.

    zlaset('F', max(M, N), NRHS, Complex.zero, Complex.zero, B, LDB);
    dlaset('F', MINMN, 1, ZERO, ZERO, S.asMatrix(1), 1);
    RANK.value = 0;

    WORK[1] = MAXWRK.toComplex();
    IWORK[1] = LIWORK;
    RWORK[1] = LRWORK.toDouble();
    return;
  }

  // Scale B if max entry outside range [SMLNUM,BIGNUM].

  BNRM = zlange('M', M, NRHS, B, LDB, RWORK);
  IBSCL = 0;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM.

    zlascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO);
    IBSCL = 1;
  } else if (BNRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM.

    zlascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO);
    IBSCL = 2;
  }

  // If M < N make sure B(M+1:N,:) = 0

  if (M < N) {
    zlaset('F', N - M, NRHS, Complex.zero, Complex.zero, B(M + 1, 1), LDB);
  }

  // Overdetermined case.

  if (M >= N) {
    // Path 1 - overdetermined or exactly determined.

    MM = M;
    if (M >= MNTHR) {
      // Path 1a - overdetermined, with many more rows than columns

      MM = N;
      ITAU = 1;
      NWORK = ITAU + N;

      // Compute A=Q*R.
      // (RWorkspace: need N)
      // (CWorkspace: need N, prefer N*NB)

      zgeqrf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, INFO);

      // Multiply B by transpose(Q).
      // (RWorkspace: need N)
      // (CWorkspace: need NRHS, prefer NRHS*NB)

      zunmqr('L', 'C', M, NRHS, N, A, LDA, WORK(ITAU), B, LDB, WORK(NWORK),
          LWORK - NWORK + 1, INFO);

      // Zero out below R.

      if (N > 1) {
        zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
      }
    }

    ITAUQ = 1;
    ITAUP = ITAUQ + N;
    NWORK = ITAUP + N;
    IE = 1;
    NRWORK = IE + N;

    // Bidiagonalize R in A.
    // (RWorkspace: need N)
    // (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)

    zgebrd(MM, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(NWORK),
        LWORK - NWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors of R.
    // (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)

    zunmbr('Q', 'L', 'C', MM, NRHS, N, A, LDA, WORK(ITAUQ), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);

    // Solve the bidiagonal least squares problem.

    zlalsd('U', SMLSIZ, N, NRHS, S, RWORK(IE), B, LDB, RCOND, RANK, WORK(NWORK),
        RWORK(NRWORK), IWORK, INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toComplex();
      IWORK[1] = LIWORK;
      RWORK[1] = LRWORK.toDouble();
      return;
    }

    // Multiply B by right bidiagonalizing vectors of R.

    zunmbr('P', 'L', 'N', N, NRHS, N, A, LDA, WORK(ITAUP), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);
  } else if (N >= MNTHR &&
      LWORK >= 4 * M + M * M + max(max(M, 2 * M - 4), max(NRHS, N - 3 * M))) {
    // Path 2a - underdetermined, with many more columns than rows
    // and sufficient workspace for an efficient algorithm.

    LDWORK = M;
    if (LWORK >=
        max(4 * M + M * LDA + max(max(M, 2 * M - 4), max(NRHS, N - 3 * M)),
            M * LDA + M + M * NRHS)) LDWORK = LDA;
    ITAU = 1;
    NWORK = M + 1;

    // Compute A=L*Q.
    // (CWorkspace: need 2*M, prefer M+M*NB)

    zgelqf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, INFO);
    IL = NWORK;

    // Copy L to WORK(IL), zeroing out above its diagonal.

    zlacpy('L', M, M, A, LDA, WORK(IL).asMatrix(LDWORK), LDWORK);
    zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
        WORK(IL + LDWORK).asMatrix(LDWORK), LDWORK);
    ITAUQ = IL + LDWORK * M;
    ITAUP = ITAUQ + M;
    NWORK = ITAUP + M;
    IE = 1;
    NRWORK = IE + M;

    // Bidiagonalize L in WORK(IL).
    // (RWorkspace: need M)
    // (CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB)

    zgebrd(M, M, WORK(IL).asMatrix(LDWORK), LDWORK, S, RWORK(IE), WORK(ITAUQ),
        WORK(ITAUP), WORK(NWORK), LWORK - NWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors of L.
    // (CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)

    zunmbr('Q', 'L', 'C', M, NRHS, M, WORK(IL).asMatrix(LDWORK), LDWORK,
        WORK(ITAUQ), B, LDB, WORK(NWORK), LWORK - NWORK + 1, INFO);

    // Solve the bidiagonal least squares problem.

    zlalsd('U', SMLSIZ, M, NRHS, S, RWORK(IE), B, LDB, RCOND, RANK, WORK(NWORK),
        RWORK(NRWORK), IWORK, INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toComplex();
      IWORK[1] = LIWORK;
      RWORK[1] = LRWORK.toDouble();
      return;
    }

    // Multiply B by right bidiagonalizing vectors of L.

    zunmbr('P', 'L', 'N', M, NRHS, M, WORK(IL).asMatrix(LDWORK), LDWORK,
        WORK(ITAUP), B, LDB, WORK(NWORK), LWORK - NWORK + 1, INFO);

    // Zero out below first M rows of B.

    zlaset('F', N - M, NRHS, Complex.zero, Complex.zero, B(M + 1, 1), LDB);
    NWORK = ITAU + M;

    // Multiply transpose(Q) by B.
    // (CWorkspace: need NRHS, prefer NRHS*NB)

    zunmlq('L', 'C', N, NRHS, M, A, LDA, WORK(ITAU), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);
  } else {
    // Path 2 - remaining underdetermined cases.

    ITAUQ = 1;
    ITAUP = ITAUQ + M;
    NWORK = ITAUP + M;
    IE = 1;
    NRWORK = IE + M;

    // Bidiagonalize A.
    // (RWorkspace: need M)
    // (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)

    zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(NWORK),
        LWORK - NWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors.
    // (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)

    zunmbr('Q', 'L', 'C', M, NRHS, N, A, LDA, WORK(ITAUQ), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);

    // Solve the bidiagonal least squares problem.

    zlalsd('L', SMLSIZ, M, NRHS, S, RWORK(IE), B, LDB, RCOND, RANK, WORK(NWORK),
        RWORK(NRWORK), IWORK, INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toComplex();
      IWORK[1] = LIWORK;
      RWORK[1] = LRWORK.toDouble();
      return;
    }

    // Multiply B by right bidiagonalizing vectors of A.

    zunmbr('P', 'L', 'N', N, NRHS, M, A, LDA, WORK(ITAUP), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);
  }

  // Undo scaling.

  if (IASCL == 1) {
    zlascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO);
    dlascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, INFO);
  } else if (IASCL == 2) {
    zlascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO);
    dlascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, INFO);
  }
  if (IBSCL == 1) {
    zlascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO);
  } else if (IBSCL == 2) {
    zlascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO);
  }
  WORK[1] = MAXWRK.toComplex();
  IWORK[1] = LIWORK;
  RWORK[1] = LRWORK.toDouble();
}
