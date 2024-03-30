import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgelqf.dart';
import 'package:lapack/src/dgeqp3.dart';
import 'package:lapack/src/dgeqrf.dart';
import 'package:lapack/src/dgesvj.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/dlaswp.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/dormlq.dart';
import 'package:lapack/src/dormqr.dart';
import 'package:lapack/src/dpocon.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgejsv(
  final String JOBA,
  final String JOBU,
  final String JOBV,
  final String JOBR,
  final String JOBT,
  final String JOBP,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> SVA_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> V_,
  final int LDV,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final SVA = SVA_.having();
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  double AATMAX,
      AATMIN,
      BIG,
      BIG1,
      COND_OK,
      CONDR1,
      CONDR2,
      ENTRA,
      ENTRAT,
      EPSLN,
      MAXPRJ,
      SCALEM,
      SCONDA,
      SFMIN,
      SMALL,
      USCAL1,
      USCAL2;
  int N1 = 0, NR, NUMRANK, p, q, WARNING;
  bool ALMORT,
      DEFR,
      ERREST,
      GOSCAL,
      JRACC,
      KILL,
      LSVEC,
      L2ABER,
      L2KILL,
      L2PERT,
      L2RANK,
      L2TRAN,
      NOSCAL,
      ROWPIV,
      RSVEC,
      TRANSP;
  final AAPP = Box(0.0), AAQQ = Box(0.0), XSC = Box(0.0), TEMP1 = Box(0.0);
  final IERR = Box(0);

  // Test the input arguments

  LSVEC = lsame(JOBU, 'U') || lsame(JOBU, 'F');
  JRACC = lsame(JOBV, 'J');
  RSVEC = lsame(JOBV, 'V') || JRACC;
  ROWPIV = lsame(JOBA, 'F') || lsame(JOBA, 'G');
  L2RANK = lsame(JOBA, 'R');
  L2ABER = lsame(JOBA, 'A');
  ERREST = lsame(JOBA, 'E') || lsame(JOBA, 'G');
  L2TRAN = lsame(JOBT, 'T');
  L2KILL = lsame(JOBR, 'R');
  DEFR = lsame(JOBR, 'N');
  L2PERT = lsame(JOBP, 'P');

  if (!(ROWPIV || L2RANK || L2ABER || ERREST || lsame(JOBA, 'C'))) {
    INFO.value = -1;
  } else if (!(LSVEC || lsame(JOBU, 'N') || lsame(JOBU, 'W'))) {
    INFO.value = -2;
  } else if (!(RSVEC || lsame(JOBV, 'N') || lsame(JOBV, 'W')) ||
      (JRACC && !LSVEC)) {
    INFO.value = -3;
  } else if (!(L2KILL || DEFR)) {
    INFO.value = -4;
  } else if (!(L2TRAN || lsame(JOBT, 'N'))) {
    INFO.value = -5;
  } else if (!(L2PERT || lsame(JOBP, 'N'))) {
    INFO.value = -6;
  } else if (M < 0) {
    INFO.value = -7;
  } else if ((N < 0) || (N > M)) {
    INFO.value = -8;
  } else if (LDA < M) {
    INFO.value = -10;
  } else if (LSVEC && (LDU < M)) {
    INFO.value = -13;
  } else if (RSVEC && (LDV < N)) {
    INFO.value = -15;
  } else if ((!(LSVEC || RSVEC || ERREST) &&
          (LWORK < max(7, max(4 * N + 1, 2 * M + N)))) ||
      (!(LSVEC || RSVEC) &&
          ERREST &&
          (LWORK < max(7, max(4 * N + N * N, 2 * M + N)))) ||
      (LSVEC && !RSVEC && (LWORK < max(7, max(2 * M + N, 4 * N + 1)))) ||
      (RSVEC && !LSVEC && (LWORK < max(7, max(2 * M + N, 4 * N + 1)))) ||
      (LSVEC &&
          RSVEC &&
          !JRACC &&
          (LWORK < max(2 * M + N, 6 * N + 2 * N * N))) ||
      (LSVEC &&
          RSVEC &&
          JRACC &&
          LWORK < max(2 * M + N, max(4 * N + N * N, 2 * N + N * N + 6)))) {
    INFO.value = -17;
  } else {
    INFO.value = 0;
  }

  if (INFO.value != 0) {
    xerbla('DGEJSV', -INFO.value);
    return;
  }

  // Quick return for void matrix (Y3K safe)
  if ((M == 0) || (N == 0)) {
    for (var i = 1; i <= 3; i++) {
      IWORK[i] = 0;
    }
    for (var i = 1; i <= 7; i++) {
      WORK[i] = 0;
    }
    return;
  }

  // Determine whether the matrix U should be M x N or M x M

  if (LSVEC) {
    N1 = N;
    if (lsame(JOBU, 'F')) N1 = M;
  }

  // Set numerical parameters

  // NOTE: Make sure dlamch() does not fail on the target architecture.

  EPSLN = dlamch('Epsilon');
  SFMIN = dlamch('SafeMinimum');
  SMALL = SFMIN / EPSLN;
  BIG = dlamch('O');
  // BIG   = ONE / SFMIN

  // Initialize SVA[1:N] = diag( ||A e_i||_2 )_1^N

  // If necessary, scale SVA() to protect the largest norm from
  // overflow. It is possible that this scaling pushes the smallest
  // column norm left from the underflow threshold (extreme case).

  SCALEM = ONE / sqrt(M * N);
  NOSCAL = true;
  GOSCAL = true;
  for (p = 1; p <= N; p++) {
    AAPP.value = ZERO;
    AAQQ.value = ONE;
    dlassq(M, A(1, p).asArray(), 1, AAPP, AAQQ);
    if (AAPP.value > BIG) {
      INFO.value = -9;
      xerbla('DGEJSV', -INFO.value);
      return;
    }
    AAQQ.value = sqrt(AAQQ.value);
    if ((AAPP.value < (BIG / AAQQ.value)) && NOSCAL) {
      SVA[p] = AAPP.value * AAQQ.value;
    } else {
      NOSCAL = false;
      SVA[p] = AAPP.value * (AAQQ.value * SCALEM);
      if (GOSCAL) {
        GOSCAL = false;
        dscal(p - 1, SCALEM, SVA, 1);
      }
    }
  }

  if (NOSCAL) SCALEM = ONE;

  AAPP.value = ZERO;
  AAQQ.value = BIG;
  for (p = 1; p <= N; p++) {
    AAPP.value = max(AAPP.value, SVA[p]);
    if (SVA[p] != ZERO) AAQQ.value = min(AAQQ.value, SVA[p]);
  }

  // Quick return for zero M x N matrix
  if (AAPP.value == ZERO) {
    if (LSVEC) dlaset('G', M, N1, ZERO, ONE, U, LDU);
    if (RSVEC) dlaset('G', N, N, ZERO, ONE, V, LDV);
    WORK[1] = ONE;
    WORK[2] = ONE;
    if (ERREST) WORK[3] = ONE;
    if (LSVEC && RSVEC) {
      WORK[4] = ONE;
      WORK[5] = ONE;
    }
    if (L2TRAN) {
      WORK[6] = ZERO;
      WORK[7] = ZERO;
    }
    IWORK[1] = 0;
    IWORK[2] = 0;
    IWORK[3] = 0;
    return;
  }

  // Issue warning if denormalized column norms detected. Override the
  // high relative accuracy request. Issue licence to kill columns
  // (set them to zero) whose norm is less than sigma_max / BIG (roughly).
  WARNING = 0;
  if (AAQQ.value <= SFMIN) {
    L2RANK = true;
    L2KILL = true;
    WARNING = 1;
  }

  // Quick return for one-column matrix
  if (N == 1) {
    if (LSVEC) {
      dlascl('G', 0, 0, SVA[1], SCALEM, M, 1, A(1, 1), LDA, IERR);
      dlacpy('A', M, 1, A, LDA, U, LDU);
      // computing all M left singular vectors of the M x 1 matrix
      if (N1 != N) {
        dgeqrf(M, N, U, LDU, WORK, WORK(N + 1), LWORK - N, IERR);
        dorgqr(M, N1, 1, U, LDU, WORK, WORK(N + 1), LWORK - N, IERR);
        dcopy(M, A(1, 1).asArray(), 1, U(1, 1).asArray(), 1);
      }
    }
    if (RSVEC) {
      V[1][1] = ONE;
    }
    if (SVA[1] < (BIG * SCALEM)) {
      SVA[1] /= SCALEM;
      SCALEM = ONE;
    }
    WORK[1] = ONE / SCALEM;
    WORK[2] = ONE;
    if (SVA[1] != ZERO) {
      IWORK[1] = 1;
      if ((SVA[1] / SCALEM) >= SFMIN) {
        IWORK[2] = 1;
      } else {
        IWORK[2] = 0;
      }
    } else {
      IWORK[1] = 0;
      IWORK[2] = 0;
    }
    IWORK[3] = 0;
    if (ERREST) WORK[3] = ONE;
    if (LSVEC && RSVEC) {
      WORK[4] = ONE;
      WORK[5] = ONE;
    }
    if (L2TRAN) {
      WORK[6] = ZERO;
      WORK[7] = ZERO;
    }
    return;
  }

  TRANSP = false;
  L2TRAN = L2TRAN && (M == N);

  AATMAX = -ONE;
  AATMIN = BIG;
  if (ROWPIV || L2TRAN) {
    // Compute the row norms, needed to determine row pivoting sequence
    // (in the case of heavily row weighted A, row pivoting is strongly
    // advised) and to collect information needed to compare the
    // structures of A * A^t and A^t * A (in the case L2TRAN == true ).

    if (L2TRAN) {
      for (p = 1; p <= M; p++) {
        XSC.value = ZERO;
        TEMP1.value = ONE;
        dlassq(N, A(p, 1).asArray(), LDA, XSC, TEMP1);
        // DLASSQ gets both the ell_2 and the ell_infinity norm
        // in one pass through the vector
        WORK[M + N + p] = XSC.value * SCALEM;
        WORK[N + p] = XSC.value * (SCALEM * sqrt(TEMP1.value));
        AATMAX = max(AATMAX, WORK[N + p]);
        if (WORK[N + p] != ZERO) AATMIN = min(AATMIN, WORK[N + p]);
      }
    } else {
      for (p = 1; p <= M; p++) {
        WORK[M + N + p] =
            SCALEM * A[p][idamax(N, A(p, 1).asArray(), LDA)].abs();
        AATMAX = max(AATMAX, WORK[M + N + p]);
        AATMIN = min(AATMIN, WORK[M + N + p]);
      }
    }
  }

  // For square matrix A try to determine whether A^t  would be  better
  // input for the preconditioned Jacobi SVD, with faster convergence.
  // The decision is based on an O(N) function of the vector of column
  // and row norms of A, based on the Shannon entropy. This should give
  // the right choice in most cases when the difference actually matters.
  // It may fail and pick the slower converging side.

  ENTRA = ZERO;
  ENTRAT = ZERO;
  if (L2TRAN) {
    XSC.value = ZERO;
    TEMP1.value = ONE;
    dlassq(N, SVA, 1, XSC, TEMP1);
    TEMP1.value = ONE / TEMP1.value;

    ENTRA = ZERO;
    for (p = 1; p <= N; p++) {
      BIG1 = (pow((SVA[p] / XSC.value), 2)) * TEMP1.value;
      if (BIG1 != ZERO) ENTRA += BIG1 * log(BIG1);
    }
    ENTRA = -ENTRA / log(N);

    // Now, SVA().^2/Trace(A^t * A) is a point in the probability simplex.
    // It is derived from the diagonal of  A^t * A.  Do the same with the
    // diagonal of A * A^t, compute the entropy of the corresponding
    // probability distribution. Note that A * A^t and A^t * A have the
    // same trace.

    ENTRAT = ZERO;
    for (p = N + 1; p <= N + M; p++) {
      BIG1 = (pow((WORK[p] / XSC.value), 2)) * TEMP1.value;
      if (BIG1 != ZERO) ENTRAT += BIG1 * log(BIG1);
    }
    ENTRAT = -ENTRAT / log(M);

    // Analyze the entropies and decide A or A^t. Smaller entropy
    // usually means better input for the algorithm.

    TRANSP = (ENTRAT < ENTRA);

    // If A^t is better than A, transpose A.

    if (TRANSP) {
      // In an optimal implementation, this trivial transpose
      // should be replaced with faster transpose.
      for (p = 1; p <= N - 1; p++) {
        for (q = p + 1; q <= N; q++) {
          TEMP1.value = A[q][p];
          A[q][p] = A[p][q];
          A[p][q] = TEMP1.value;
        }
      }
      for (p = 1; p <= N; p++) {
        WORK[M + N + p] = SVA[p];
        SVA[p] = WORK[N + p];
      }
      TEMP1.value = AAPP.value;
      AAPP.value = AATMAX;
      AATMAX = TEMP1.value;
      TEMP1.value = AAQQ.value;
      AAQQ.value = AATMIN;
      AATMIN = TEMP1.value;
      KILL = LSVEC;
      LSVEC = RSVEC;
      RSVEC = KILL;
      if (LSVEC) N1 = N;

      ROWPIV = true;
    }
  }
  // END IF L2TRAN

  // Scale the matrix so that its maximal singular value remains less
  // than sqrt(BIG) -- the matrix is scaled so that its maximal column
  // has Euclidean norm equal to sqrt(BIG/N). The only reason to keep
  // sqrt(BIG) instead of BIG is the fact that DGEJSV uses LAPACK and
  // BLAS routines that, in some implementations, are not capable of
  // working in the full interval [SFMIN,BIG] and that they may provoke
  // overflows in the intermediate results. If the singular values spread
  // from SFMIN to BIG, then DGESVJ will compute them. So, in that case,
  // one should use DGESVJ instead of DGEJSV.

  BIG1 = sqrt(BIG);
  TEMP1.value = sqrt(BIG / N);

  dlascl('G', 0, 0, AAPP.value, TEMP1.value, N, 1, SVA.asMatrix(N), N, IERR);
  if (AAQQ.value > (AAPP.value * SFMIN)) {
    AAQQ.value = (AAQQ.value / AAPP.value) * TEMP1.value;
  } else {
    AAQQ.value = (AAQQ.value * TEMP1.value) / AAPP.value;
  }
  TEMP1.value *= SCALEM;
  dlascl('G', 0, 0, AAPP.value, TEMP1.value, M, N, A, LDA, IERR);

  // To undo scaling at the end of this procedure, multiply the
  // computed singular values with USCAL2 / USCAL1.

  USCAL1 = TEMP1.value;
  USCAL2 = AAPP.value;

  if (L2KILL) {
    // L2KILL enforces computation of nonzero singular values in
    // the restricted range of condition number of the initial A,
    // sigma_max(A) / sigma_min(A) approx. sqrt(BIG)/sqrt(SFMIN).
    XSC.value = sqrt(SFMIN);
  } else {
    XSC.value = SMALL;

    // Now, if the condition number of A is too big,
    // sigma_max(A) / sigma_min(A) > sqrt(BIG/N) * EPSLN / SFMIN,
    // as a precaution measure, the full SVD is computed using DGESVJ
    // with accumulated Jacobi rotations. This provides numerically
    // more robust computation, at the cost of slightly increased run
    // time. Depending on the concrete implementation of BLAS and LAPACK
    // (i.e. how they behave in presence of extreme ill-conditioning) the
    // implementor may decide to remove this switch.
    if ((AAQQ.value < sqrt(SFMIN)) && LSVEC && RSVEC) {
      JRACC = true;
    }
  }
  if (AAQQ.value < XSC.value) {
    for (p = 1; p <= N; p++) {
      if (SVA[p] < XSC.value) {
        dlaset('A', M, 1, ZERO, ZERO, A(1, p), LDA);
        SVA[p] = ZERO;
      }
    }
  }

  // Preconditioning using QR factorization with pivoting

  if (ROWPIV) {
    // Optional row permutation (Bjoerck row pivoting):
    // A result by Cox and Higham shows that the Bjoerck's
    // row pivoting combined with standard column pivoting
    // has similar effect as Powell-Reid complete pivoting.
    // The ell-infinity norms of A are made nonincreasing.
    for (p = 1; p <= M - 1; p++) {
      q = idamax(M - p + 1, WORK(M + N + p), 1) + p - 1;
      IWORK[2 * N + p] = q;
      if (p != q) {
        TEMP1.value = WORK[M + N + p];
        WORK[M + N + p] = WORK[M + N + q];
        WORK[M + N + q] = TEMP1.value;
      }
    }
    dlaswp(N, A, LDA, 1, M - 1, IWORK(2 * N + 1), 1);
  }

  // End of the preparation phase (scaling, optional sorting and
  // transposing, optional flushing of small columns).

  // Preconditioning

  // If the full SVD is needed, the right singular vectors are computed
  // from a matrix equation, and for that we need theoretical analysis
  // of the Businger-Golub pivoting. So we use DGEQP3 as the first RR QRF.
  // In all other cases the first RR QRF can be chosen by other criteria
  // (eg speed by replacing global with restricted window pivoting, such
  // as in SGEQPX from TOMS # 782). Good results will be obtained using
  // SGEQPX with properly (!) chosen numerical parameters.
  // Any improvement of DGEQP3 improves overall performance of DGEJSV.

  // A * P1 = Q1 * [ R1^t 0]^t:
  for (p = 1; p <= N; p++) {
    // .. all columns are free columns
    IWORK[p] = 0;
  }
  dgeqp3(M, N, A, LDA, IWORK, WORK, WORK(N + 1), LWORK - N, IERR);

  // The upper triangular matrix R1 from the first QRF is inspected for
  // rank deficiency and possibilities for deflation, or possible
  // ill-conditioning. Depending on the user specified flag L2RANK,
  // the procedure explores possibilities to reduce the numerical
  // rank by inspecting the computed upper triangular factor. If
  // L2RANK or L2ABER are up, then DGEJSV will compute the SVD of
  // A + dA, where ||dA|| <= f(M,N)*EPSLN.

  NR = 1;
  if (L2ABER) {
    // Standard absolute error bound suffices. All sigma_i with
    // sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
    // aggressive enforcement of lower numerical rank by introducing a
    // backward error of the order of N*EPSLN*||A||.
    TEMP1.value = sqrt(N) * EPSLN;
    for (p = 2; p <= N; p++) {
      if (A[p][p].abs() >= (TEMP1.value * A[1][1].abs())) {
        NR++;
      } else {
        break;
      }
    }
  } else if (L2RANK) {
    // .. similarly as above, only slightly more gentle (less aggressive).
    // Sudden drop on the diagonal of R1 is used as the criterion for
    // close-to-rank-deficient.
    TEMP1.value = sqrt(SFMIN);
    for (p = 2; p <= N; p++) {
      if ((A[p][p].abs() < (EPSLN * A[p - 1][p - 1].abs())) ||
          (A[p][p].abs() < SMALL) ||
          (L2KILL && (A[p][p].abs() < TEMP1.value))) break;
      NR++;
    }
  } else {
    // The goal is high relative accuracy. However, if the matrix
    // has high scaled condition number the relative accuracy is in
    // general not feasible. Later on, a condition number estimator
    // will be deployed to estimate the scaled condition number.
    // Here we just remove the underflowed part of the triangular
    // factor. This prevents the situation in which the code is
    // working hard to get the accuracy not warranted by the data.
    TEMP1.value = sqrt(SFMIN);
    for (p = 2; p <= N; p++) {
      if ((A[p][p].abs() < SMALL) ||
          (L2KILL && (A[p][p].abs() < TEMP1.value))) {
        break;
      }
      NR++;
    }
  }

  ALMORT = false;
  if (NR == N) {
    MAXPRJ = ONE;
    for (p = 2; p <= N; p++) {
      TEMP1.value = A[p][p].abs() / SVA[IWORK[p]];
      MAXPRJ = min(MAXPRJ, TEMP1.value);
    }
    if (pow(MAXPRJ, 2) >= ONE - N * EPSLN) ALMORT = true;
  }

  SCONDA = -ONE;
  CONDR1 = -ONE;
  CONDR2 = -ONE;

  if (ERREST) {
    if (N == NR) {
      if (RSVEC) {
        // .. V is available as workspace
        dlacpy('U', N, N, A, LDA, V, LDV);
        for (p = 1; p <= N; p++) {
          TEMP1.value = SVA[IWORK[p]];
          dscal(p, ONE / TEMP1.value, V(1, p).asArray(), 1);
        }
        dpocon('U', N, V, LDV, ONE, TEMP1, WORK(N + 1), IWORK(2 * N + M + 1),
            IERR);
      } else if (LSVEC) {
        // .. U is available as workspace
        dlacpy('U', N, N, A, LDA, U, LDU);
        for (p = 1; p <= N; p++) {
          TEMP1.value = SVA[IWORK[p]];
          dscal(p, ONE / TEMP1.value, U(1, p).asArray(), 1);
        }
        dpocon('U', N, U, LDU, ONE, TEMP1, WORK(N + 1), IWORK(2 * N + M + 1),
            IERR);
      } else {
        dlacpy('U', N, N, A, LDA, WORK(N + 1).asMatrix(N), N);
        for (p = 1; p <= N; p++) {
          TEMP1.value = SVA[IWORK[p]];
          dscal(p, ONE / TEMP1.value, WORK(N + (p - 1) * N + 1), 1);
        }
        // .. the columns of R are scaled to have unit Euclidean lengths.
        dpocon('U', N, WORK(N + 1).asMatrix(N), N, ONE, TEMP1,
            WORK(N + N * N + 1), IWORK(2 * N + M + 1), IERR);
      }
      SCONDA = ONE / sqrt(TEMP1.value);
      // SCONDA is an estimate of sqrt(||(R^t * R)^(-1)||_1).
      // N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
    } else {
      SCONDA = -ONE;
    }
  }

  L2PERT = L2PERT && ((A[1][1] / A[NR][NR]).abs() > sqrt(BIG1));
  // If there is no violent scaling, artificial perturbation is not needed.

  // Phase 3:

  if (!(RSVEC || LSVEC)) {
    // Singular Values only

    // .. transpose A[1:NR][1:N]
    for (p = 1; p <= min(N - 1, NR); p++) {
      dcopy(N - p, A(p, p + 1).asArray(), LDA, A(p + 1, p).asArray(), 1);
    }

    // The following two DO-loops introduce small relative perturbation
    // into the strict upper triangle of the lower triangular matrix.
    // Small entries below the main diagonal are also changed.
    // This modification is useful if the computing environment does not
    // provide/allow FLUSH TO ZERO underflow, for it prevents many
    // annoying denormalized numbers in case of strongly scaled matrices.
    // The perturbation is structured so that it does not introduce any
    // new perturbation of the singular values, and it does not destroy
    // the job done by the preconditioner.
    // The licence for this perturbation is in the variable L2PERT, which
    // should be false if FLUSH TO ZERO underflow is active.

    if (!ALMORT) {
      if (L2PERT) {
        // XSC.value = sqrt(SMALL)
        XSC.value = EPSLN / N;
        for (q = 1; q <= NR; q++) {
          TEMP1.value = XSC.value * A[q][q].abs();
          for (p = 1; p <= N; p++) {
            if (((p > q) && (A[p][q].abs() <= TEMP1.value)) || (p < q)) {
              A[p][q] = sign(TEMP1.value, A[p][q]);
            }
          }
        }
      } else {
        dlaset('U', NR - 1, NR - 1, ZERO, ZERO, A(1, 2), LDA);
      }

      // .. second preconditioning using the QR factorization

      dgeqrf(N, NR, A, LDA, WORK, WORK(N + 1), LWORK - N, IERR);

      // .. and transpose upper to lower triangular
      for (p = 1; p <= NR - 1; p++) {
        dcopy(NR - p, A(p, p + 1).asArray(), LDA, A(p + 1, p).asArray(), 1);
      }
    }

    // Row-cyclic Jacobi SVD algorithm with column pivoting

    // .. again some perturbation (a "background noise") is added
    // to drown denormals
    if (L2PERT) {
      // XSC.value = sqrt(SMALL)
      XSC.value = EPSLN / N;
      for (q = 1; q <= NR; q++) {
        TEMP1.value = XSC.value * A[q][q].abs();
        for (p = 1; p <= NR; p++) {
          if (((p > q) && (A[p][q].abs() <= TEMP1.value)) || (p < q)) {
            A[p][q] = sign(TEMP1.value, A[p][q]);
          }
        }
      }
    } else {
      dlaset('U', NR - 1, NR - 1, ZERO, ZERO, A(1, 2), LDA);
    }

    // .. and one-sided Jacobi rotations are started on a lower
    // triangular matrix (plus perturbation which is ignored in
    // the part which destroys triangular form (confusing?!))

    dgesvj(
        'L', 'NoU', 'NoV', NR, NR, A, LDA, SVA, N, V, LDV, WORK, LWORK, INFO);

    SCALEM = WORK[1];
    NUMRANK = WORK[2].round();
  } else if (RSVEC && !LSVEC) {
    // -> Singular Values and Right Singular Vectors <-

    if (ALMORT) {
      // .. in this case NR equals N
      for (p = 1; p <= NR; p++) {
        dcopy(N - p + 1, A(p, p).asArray(), LDA, V(p, p).asArray(), 1);
      }
      dlaset('Upper', NR - 1, NR - 1, ZERO, ZERO, V(1, 2), LDV);

      dgesvj('L', 'U', 'N', N, NR, V, LDV, SVA, NR, A, LDA, WORK, LWORK, INFO);
      SCALEM = WORK[1];
      NUMRANK = WORK[2].round();
    } else {
      // .. two more QR factorizations ( one QRF is not enough, two require
      // accumulated product of Jacobi rotations, three are perfect )

      dlaset('Lower', NR - 1, NR - 1, ZERO, ZERO, A(2, 1), LDA);
      dgelqf(NR, N, A, LDA, WORK, WORK(N + 1), LWORK - N, IERR);
      dlacpy('Lower', NR, NR, A, LDA, V, LDV);
      dlaset('Upper', NR - 1, NR - 1, ZERO, ZERO, V(1, 2), LDV);
      dgeqrf(NR, NR, V, LDV, WORK(N + 1), WORK(2 * N + 1), LWORK - 2 * N, IERR);
      for (p = 1; p <= NR; p++) {
        dcopy(NR - p + 1, V(p, p).asArray(), LDV, V(p, p).asArray(), 1);
      }
      dlaset('Upper', NR - 1, NR - 1, ZERO, ZERO, V(1, 2), LDV);

      dgesvj('Lower', 'U', 'N', NR, NR, V, LDV, SVA, NR, U, LDU, WORK(N + 1),
          LWORK, INFO);
      SCALEM = WORK[N + 1];
      NUMRANK = WORK[N + 2].round();
      if (NR < N) {
        dlaset('A', N - NR, NR, ZERO, ZERO, V(NR + 1, 1), LDV);
        dlaset('A', NR, N - NR, ZERO, ZERO, V(1, NR + 1), LDV);
        dlaset('A', N - NR, N - NR, ZERO, ONE, V(NR + 1, NR + 1), LDV);
      }

      dormlq('Left', 'Transpose', N, N, NR, A, LDA, WORK, V, LDV, WORK(N + 1),
          LWORK - N, IERR);
    }

    for (p = 1; p <= N; p++) {
      dcopy(N, V(p, 1).asArray(), LDV, A(IWORK[p], 1).asArray(), LDA);
    }
    dlacpy('All', N, N, A, LDA, V, LDV);

    if (TRANSP) {
      dlacpy('All', N, N, V, LDV, U, LDU);
    }
  } else if (LSVEC && !RSVEC) {
    // .. Singular Values and Left Singular Vectors                 ..

    // .. second preconditioning step to avoid need to accumulate
    // Jacobi rotations in the Jacobi iterations.
    for (p = 1; p <= NR; p++) {
      dcopy(N - p + 1, A(p, p).asArray(), LDA, U(p, p).asArray(), 1);
    }
    dlaset('Upper', NR - 1, NR - 1, ZERO, ZERO, U(1, 2), LDU);

    dgeqrf(N, NR, U, LDU, WORK(N + 1), WORK(2 * N + 1), LWORK - 2 * N, IERR);

    for (p = 1; p <= NR - 1; p++) {
      dcopy(NR - p, U(p, p + 1).asArray(), LDU, U(p + 1, p).asArray(), 1);
    }
    dlaset('Upper', NR - 1, NR - 1, ZERO, ZERO, U(1, 2), LDU);

    dgesvj('Lower', 'U', 'N', NR, NR, U, LDU, SVA, NR, A, LDA, WORK(N + 1),
        LWORK - N, INFO);
    SCALEM = WORK[N + 1];
    NUMRANK = WORK[N + 2].round();

    if (NR < M) {
      dlaset('A', M - NR, NR, ZERO, ZERO, U(NR + 1, 1), LDU);
      if (NR < N1) {
        dlaset('A', NR, N1 - NR, ZERO, ZERO, U(1, NR + 1), LDU);
        dlaset('A', M - NR, N1 - NR, ZERO, ONE, U(NR + 1, NR + 1), LDU);
      }
    }

    dormqr('Left', 'No Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N + 1),
        LWORK - N, IERR);

    if (ROWPIV) dlaswp(N1, U, LDU, 1, M - 1, IWORK(2 * N + 1), -1);

    for (p = 1; p <= N1; p++) {
      XSC.value = ONE / dnrm2(M, U(1, p).asArray(), 1);
      dscal(M, XSC.value, U(1, p).asArray(), 1);
    }

    if (TRANSP) {
      dlacpy('All', N, N, U, LDU, V, LDV);
    }
  } else {
    // .. Full SVD ..

    if (!JRACC) {
      if (!ALMORT) {
        // Second Preconditioning Step (QRF [with pivoting])
        // Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
        // equivalent to an LQF CALL. Since in many libraries the QRF
        // seems to be better optimized than the LQF, we do explicit
        // transpose and use the QRF. This is subject to changes in an
        // optimized implementation of DGEJSV.

        for (p = 1; p <= NR; p++) {
          dcopy(N - p + 1, A(p, p).asArray(), LDA, V(p, p).asArray(), 1);
        }

        // .. the following two loops perturb small entries to avoid
        // denormals in the second QR factorization, where they are
        // as good as zeros. This is done to avoid painfully slow
        // computation with denormals. The relative size of the perturbation
        // is a parameter that can be changed by the implementer.
        // This perturbation device will be obsolete on machines with
        // properly implemented arithmetic.
        // To switch it off, set L2PERT= false To remove it from  the
        // code, remove the action under L2PERT= true , leave the ELSE part.
        // The following two loops should be blocked and fused with the
        // transposed copy above.

        if (L2PERT) {
          XSC.value = sqrt(SMALL);
          for (q = 1; q <= NR; q++) {
            TEMP1.value = XSC.value * V[q][q].abs();
            for (p = 1; p <= N; p++) {
              if ((p > q) && (V[p][q].abs() <= TEMP1.value) || (p < q)) {
                V[p][q] = sign(TEMP1.value, V[p][q]);
              }
              if (p < q) V[p][q] = -V[p][q];
            }
          }
        } else {
          dlaset('U', NR - 1, NR - 1, ZERO, ZERO, V(1, 2), LDV);
        }

        // Estimate the row scaled condition number of R1
        // (If R1 is rectangular, N > NR, then the condition number
        // of the leading NR x NR submatrix is estimated.)

        dlacpy('L', NR, NR, V, LDV, WORK(2 * N + 1).asMatrix(NR), NR);
        for (p = 1; p <= NR; p++) {
          TEMP1.value = dnrm2(NR - p + 1, WORK(2 * N + (p - 1) * NR + p), 1);
          dscal(
              NR - p + 1, ONE / TEMP1.value, WORK(2 * N + (p - 1) * NR + p), 1);
        }
        dpocon('Lower', NR, WORK(2 * N + 1).asMatrix(NR), NR, ONE, TEMP1,
            WORK(2 * N + NR * NR + 1), IWORK(M + 2 * N + 1), IERR);
        CONDR1 = ONE / sqrt(TEMP1.value);
        // .. here need a second opinion on the condition number
        // .. then assume worst case scenario
        // R1 is OK for inverse <=> CONDR1 < N
        // more conservative    <=> CONDR1 < sqrt(N)

        COND_OK = sqrt(NR);
        // [TP] COND_OK is a tuning parameter.

        if (CONDR1 < COND_OK) {
          // .. the second QRF without pivoting. Note: in an optimized
          // implementation, this QRF should be implemented as the QRF
          // of a lower triangular matrix.
          // R1^t = Q2 * R2
          dgeqrf(
              N, NR, V, LDV, WORK(N + 1), WORK(2 * N + 1), LWORK - 2 * N, IERR);

          if (L2PERT) {
            XSC.value = sqrt(SMALL) / EPSLN;
            for (p = 2; p <= NR; p++) {
              for (q = 1; q <= p - 1; q++) {
                TEMP1.value = XSC.value * min(V[p][p].abs(), V[q][q].abs());
                if (V[q][p].abs() <= TEMP1.value) {
                  V[q][p] = sign(TEMP1.value, V[q][p]);
                }
              }
            }
          }

          if (NR != N) {
            dlacpy('A', N, NR, V, LDV, WORK(2 * N + 1).asMatrix(N), N);
          }
          // .. save ...

          // .. this transposed copy should be better than naive
          for (p = 1; p <= NR - 1; p++) {
            dcopy(NR - p, V(p, p + 1).asArray(), LDV, V(p + 1, p).asArray(), 1);
          }

          CONDR2 = CONDR1;
        } else {
          // .. ill-conditioned case: second QRF with pivoting
          // Note that windowed pivoting would be equally good
          // numerically, and more run-time efficient. So, in
          // an optimal implementation, the next call to DGEQP3
          // should be replaced with eg. CALL SGEQPX (ACM TOMS #782)
          // with properly (carefully) chosen parameters.

          // R1^t * P2 = Q2 * R2
          for (p = 1; p <= NR; p++) {
            IWORK[N + p] = 0;
          }
          dgeqp3(N, NR, V, LDV, IWORK(N + 1), WORK(N + 1), WORK(2 * N + 1),
              LWORK - 2 * N, IERR);
          // CALL DGEQRF( N, NR, V, LDV, WORK[N+1], WORK[2*N+1], LWORK-2*N, IERR )
          if (L2PERT) {
            XSC.value = sqrt(SMALL);
            for (p = 2; p <= NR; p++) {
              for (q = 1; q <= p - 1; q++) {
                TEMP1.value = XSC.value * min(V[p][p].abs(), V[q][q].abs());
                if (V[q][p].abs() <= TEMP1.value) {
                  V[q][p] = sign(TEMP1.value, V[q][p]);
                }
              }
            }
          }

          dlacpy('A', N, NR, V, LDV, WORK(2 * N + 1).asMatrix(N), N);

          if (L2PERT) {
            XSC.value = sqrt(SMALL);
            for (p = 2; p <= NR; p++) {
              for (q = 1; q <= p - 1; q++) {
                TEMP1.value = XSC.value * min(V[p][p].abs(), V[q][q].abs());
                V[p][q] = -sign(TEMP1.value, V[q][p]);
              }
            }
          } else {
            dlaset('L', NR - 1, NR - 1, ZERO, ZERO, V(2, 1), LDV);
          }
          // Now, compute R2 = L3 * Q3, the LQ factorization.
          dgelqf(NR, NR, V, LDV, WORK(2 * N + N * NR + 1),
              WORK(2 * N + N * NR + NR + 1), LWORK - 2 * N - N * NR - NR, IERR);
          // .. and estimate the condition number
          dlacpy('L', NR, NR, V, LDV,
              WORK(2 * N + N * NR + NR + 1).asMatrix(NR), NR);
          for (p = 1; p <= NR; p++) {
            TEMP1.value = dnrm2(p, WORK(2 * N + N * NR + NR + p), NR);
            dscal(p, ONE / TEMP1.value, WORK(2 * N + N * NR + NR + p), NR);
          }
          dpocon(
              'L',
              NR,
              WORK(2 * N + N * NR + NR + 1).asMatrix(NR),
              NR,
              ONE,
              TEMP1,
              WORK(2 * N + N * NR + NR + NR * NR + 1),
              IWORK(M + 2 * N + 1),
              IERR);
          CONDR2 = ONE / sqrt(TEMP1.value);

          if (CONDR2 >= COND_OK) {
            // .. save the Householder vectors used for Q3
            // (this overwrites the copy of R2, as it will not be
            // needed in this branch, but it does not overwrite the
            // Huseholder vectors of Q2.).
            dlacpy('U', NR, NR, V, LDV, WORK(2 * N + 1).asMatrix(N), N);
            // .. and the rest of the information on Q3 is in
            // WORK[2*N+N*NR+1:2*N+N*NR+N]
          }
        }

        if (L2PERT) {
          XSC.value = sqrt(SMALL);
          for (q = 2; q <= NR; q++) {
            TEMP1.value = XSC.value * V[q][q];
            for (p = 1; p <= q - 1; p++) {
              // V[p][q] = - sign( TEMP1.value, V[q][p] )
              V[p][q] = -sign(TEMP1.value, V[p][q]);
            }
          }
        } else {
          dlaset('U', NR - 1, NR - 1, ZERO, ZERO, V(1, 2), LDV);
        }

        // Second preconditioning finished; continue with Jacobi SVD
        // The input matrix is lower triangular.

        // Recover the right singular vectors as solution of a well
        // conditioned triangular matrix equation.

        if (CONDR1 < COND_OK) {
          dgesvj('L', 'U', 'N', NR, NR, V, LDV, SVA, NR, U, LDU,
              WORK(2 * N + N * NR + NR + 1), LWORK - 2 * N - N * NR - NR, INFO);
          SCALEM = WORK[2 * N + N * NR + NR + 1];
          NUMRANK = WORK[2 * N + N * NR + NR + 2].round();
          for (p = 1; p <= NR; p++) {
            dcopy(NR, V(1, p).asArray(), 1, U(1, p).asArray(), 1);
            dscal(NR, SVA[p], V(1, p).asArray(), 1);
          }

          // .. pick the right matrix equation and solve it

          if (NR == N) {
            // .. best case, R1 is inverted. The solution of this matrix
            // equation is Q2*V2 = the product of the Jacobi rotations
            // used in DGESVJ, premultiplied with the orthogonal matrix
            // from the second QR factorization.
            dtrsm('L', 'U', 'N', 'N', NR, NR, ONE, A, LDA, V, LDV);
          } else {
            // .. R1 is well conditioned, but non-square. Transpose(R2)
            // is inverted to get the product of the Jacobi rotations
            // used in DGESVJ. The Q-factor from the second QR
            // factorization is then built in explicitly.
            dtrsm('L', 'U', 'T', 'N', NR, NR, ONE, WORK(2 * N + 1).asMatrix(N),
                N, V, LDV);
            if (NR < N) {
              dlaset('A', N - NR, NR, ZERO, ZERO, V(NR + 1, 1), LDV);
              dlaset('A', NR, N - NR, ZERO, ZERO, V(1, NR + 1), LDV);
              dlaset('A', N - NR, N - NR, ZERO, ONE, V(NR + 1, NR + 1), LDV);
            }
            dormqr(
                'L',
                'N',
                N,
                N,
                NR,
                WORK(2 * N + 1).asMatrix(N),
                N,
                WORK(N + 1),
                V,
                LDV,
                WORK(2 * N + N * NR + NR + 1),
                LWORK - 2 * N - N * NR - NR,
                IERR);
          }
        } else if (CONDR2 < COND_OK) {
          // .. the input matrix A is very likely a relative of
          // the Kahan matrix :)
          // The matrix R2 is inverted. The solution of the matrix equation
          // is Q3^T*V3 = the product of the Jacobi rotations (applied to
          // the lower triangular L3 from the LQ factorization of
          // R2=L3*Q3), pre-multiplied with the transposed Q3.
          dgesvj('L', 'U', 'N', NR, NR, V, LDV, SVA, NR, U, LDU,
              WORK(2 * N + N * NR + NR + 1), LWORK - 2 * N - N * NR - NR, INFO);
          SCALEM = WORK[2 * N + N * NR + NR + 1];
          NUMRANK = WORK[2 * N + N * NR + NR + 2].round();
          for (p = 1; p <= NR; p++) {
            dcopy(NR, V(1, p).asArray(), 1, U(1, p).asArray(), 1);
            dscal(NR, SVA[p], U(1, p).asArray(), 1);
          }
          dtrsm('L', 'U', 'N', 'N', NR, NR, ONE, WORK(2 * N + 1).asMatrix(N), N,
              U, LDU);
          // .. apply the permutation from the second QR factorization
          for (q = 1; q <= NR; q++) {
            for (p = 1; p <= NR; p++) {
              WORK[2 * N + N * NR + NR + IWORK[N + p]] = U[p][q];
            }
            for (p = 1; p <= NR; p++) {
              U[p][q] = WORK[2 * N + N * NR + NR + p];
            }
          }
          if (NR < N) {
            dlaset('A', N - NR, NR, ZERO, ZERO, V(NR + 1, 1), LDV);
            dlaset('A', NR, N - NR, ZERO, ZERO, V(1, NR + 1), LDV);
            dlaset('A', N - NR, N - NR, ZERO, ONE, V(NR + 1, NR + 1), LDV);
          }
          dormqr(
              'L',
              'N',
              N,
              N,
              NR,
              WORK(2 * N + 1).asMatrix(N),
              N,
              WORK(N + 1),
              V,
              LDV,
              WORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              IERR);
        } else {
          // Last line of defense.
          // .. This is a rather pathological case: no scaled condition
          // improvement after two pivoted QR factorizations. Other
          // possibility is that the rank revealing QR factorization
          // or the condition estimator has failed, or the COND_OK
          // is set very close to ONE (which is unnecessary). Normally,
          // this branch should never be executed, but in rare cases of
          // failure of the RRQR or condition estimator, the last line of
          // defense ensures that DGEJSV completes the task.
          // Compute the full SVD of L3 using DGESVJ with explicit
          // accumulation of Jacobi rotations.
          dgesvj('L', 'U', 'V', NR, NR, V, LDV, SVA, NR, U, LDU,
              WORK(2 * N + N * NR + NR + 1), LWORK - 2 * N - N * NR - NR, INFO);
          SCALEM = WORK[2 * N + N * NR + NR + 1];
          NUMRANK = WORK[2 * N + N * NR + NR + 2].round();
          if (NR < N) {
            dlaset('A', N - NR, NR, ZERO, ZERO, V(NR + 1, 1), LDV);
            dlaset('A', NR, N - NR, ZERO, ZERO, V(1, NR + 1), LDV);
            dlaset('A', N - NR, N - NR, ZERO, ONE, V(NR + 1, NR + 1), LDV);
          }
          dormqr(
              'L',
              'N',
              N,
              N,
              NR,
              WORK(2 * N + 1).asMatrix(N),
              N,
              WORK(N + 1),
              V,
              LDV,
              WORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              IERR);

          dormlq(
              'L',
              'T',
              NR,
              NR,
              NR,
              WORK(2 * N + 1).asMatrix(N),
              N,
              WORK(2 * N + N * NR + 1),
              U,
              LDU,
              WORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              IERR);
          for (q = 1; q <= NR; q++) {
            for (p = 1; p <= NR; p++) {
              WORK[2 * N + N * NR + NR + IWORK[N + p]] = U[p][q];
            }
            for (p = 1; p <= NR; p++) {
              U[p][q] = WORK[2 * N + N * NR + NR + p];
            }
          }
        }

        // Permute the rows of V using the (column) permutation from the
        // first QRF. Also, scale the columns to make them unit in
        // Euclidean norm. This applies to all cases.

        TEMP1.value = sqrt(N) * EPSLN;
        for (q = 1; q <= N; q++) {
          for (p = 1; p <= N; p++) {
            WORK[2 * N + N * NR + NR + IWORK[p]] = V[p][q];
          }
          for (p = 1; p <= N; p++) {
            V[p][q] = WORK[2 * N + N * NR + NR + p];
          }
          XSC.value = ONE / dnrm2(N, V(1, q).asArray(), 1);
          if ((XSC.value < (ONE - TEMP1.value)) ||
              (XSC.value > (ONE + TEMP1.value))) {
            dscal(N, XSC.value, V(1, q).asArray(), 1);
          }
        }
        // At this moment, V contains the right singular vectors of A.
        // Next, assemble the left singular vector matrix U (M x N).
        if (NR < M) {
          dlaset('A', M - NR, NR, ZERO, ZERO, U(NR + 1, 1), LDU);
          if (NR < N1) {
            dlaset('A', NR, N1 - NR, ZERO, ZERO, U(1, NR + 1), LDU);
            dlaset('A', M - NR, N1 - NR, ZERO, ONE, U(NR + 1, NR + 1), LDU);
          }
        }

        // The Q matrix from the first QRF is built into the left singular
        // matrix U. This applies to all cases.

        dormqr('Left', 'No_Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N + 1),
            LWORK - N, IERR);

        // The columns of U are normalized. The cost is O(M*N) flops.
        TEMP1.value = sqrt(M) * EPSLN;
        for (p = 1; p <= NR; p++) {
          XSC.value = ONE / dnrm2(M, U(1, p).asArray(), 1);
          if ((XSC.value < (ONE - TEMP1.value)) ||
              (XSC.value > (ONE + TEMP1.value))) {
            dscal(M, XSC.value, U(1, p).asArray(), 1);
          }
        }

        // If the initial QRF is computed with row pivoting, the left
        // singular vectors must be adjusted.

        if (ROWPIV) dlaswp(N1, U, LDU, 1, M - 1, IWORK(2 * N + 1), -1);
      } else {
        // .. the initial matrix A has almost orthogonal columns and
        // the second QRF is not needed

        dlacpy('Upper', N, N, A, LDA, WORK(N + 1).asMatrix(N), N);
        if (L2PERT) {
          XSC.value = sqrt(SMALL);
          for (p = 2; p <= N; p++) {
            TEMP1.value = XSC.value * WORK[N + (p - 1) * N + p];
            for (q = 1; q <= p - 1; q++) {
              WORK[N + (q - 1) * N + p] =
                  -sign(TEMP1.value, WORK[N + (p - 1) * N + q]);
            }
          }
        } else {
          dlaset('Lower', N - 1, N - 1, ZERO, ZERO, WORK(N + 2).asMatrix(N), N);
        }

        dgesvj('Upper', 'U', 'N', N, N, WORK(N + 1).asMatrix(N), N, SVA, N, U,
            LDU, WORK(N + N * N + 1), LWORK - N - N * N, INFO);

        SCALEM = WORK[N + N * N + 1];
        NUMRANK = WORK[N + N * N + 2].round();
        for (p = 1; p <= N; p++) {
          dcopy(N, WORK(N + (p - 1) * N + 1), 1, U(1, p).asArray(), 1);
          dscal(N, SVA[p], WORK(N + (p - 1) * N + 1), 1);
        }

        dtrsm('Left', 'Upper', 'NoTrans', 'No UD', N, N, ONE, A, LDA,
            WORK(N + 1).asMatrix(N), N);
        for (p = 1; p <= N; p++) {
          dcopy(N, WORK(N + p), N, V(IWORK[p], 1).asArray(), LDV);
        }
        TEMP1.value = sqrt(N) * EPSLN;
        for (p = 1; p <= N; p++) {
          XSC.value = ONE / dnrm2(N, V(1, p).asArray(), 1);
          if ((XSC.value < (ONE - TEMP1.value)) ||
              (XSC.value > (ONE + TEMP1.value))) {
            dscal(N, XSC.value, V(1, p).asArray(), 1);
          }
        }

        // Assemble the left singular vector matrix U (M x N).

        if (N < M) {
          dlaset('A', M - N, N, ZERO, ZERO, U(N + 1, 1), LDU);
          if (N < N1) {
            dlaset('A', N, N1 - N, ZERO, ZERO, U(1, N + 1), LDU);
            dlaset('A', M - N, N1 - N, ZERO, ONE, U(N + 1, N + 1), LDU);
          }
        }
        dormqr('Left', 'No Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N + 1),
            LWORK - N, IERR);
        TEMP1.value = sqrt(M) * EPSLN;
        for (p = 1; p <= N1; p++) {
          XSC.value = ONE / dnrm2(M, U(1, p).asArray(), 1);
          if ((XSC.value < (ONE - TEMP1.value)) ||
              (XSC.value > (ONE + TEMP1.value))) {
            dscal(M, XSC.value, U(1, p).asArray(), 1);
          }
        }

        if (ROWPIV) dlaswp(N1, U, LDU, 1, M - 1, IWORK(2 * N + 1), -1);
      }

      // end of the  >> almost orthogonal case <<  in the full SVD
    } else {
      // This branch deploys a preconditioned Jacobi SVD with explicitly
      // accumulated rotations. It is included as optional, mainly for
      // experimental purposes. It does perform well, and can also be used.
      // In this implementation, this branch will be automatically activated
      // if the  condition number sigma_max(A) / sigma_min(A) is predicted
      // to be greater than the overflow threshold. This is because the
      // a posteriori computation of the singular vectors assumes robust
      // implementation of BLAS and some LAPACK procedures, capable of working
      // in presence of extreme values. Since that is not always the case, ...

      for (p = 1; p <= NR; p++) {
        dcopy(N - p + 1, A(p, p).asArray(), LDA, V(p, p).asArray(), 1);
      }

      if (L2PERT) {
        XSC.value = sqrt(SMALL / EPSLN);
        for (q = 1; q <= NR; q++) {
          TEMP1.value = XSC.value * V[q][q].abs();
          for (p = 1; p <= N; p++) {
            if ((p > q) && (V[p][q].abs() <= TEMP1.value) || (p < q)) {
              V[p][q] = sign(TEMP1.value, V[p][q]);
            }
            if (p < q) V[p][q] = -V[p][q];
          }
        }
      } else {
        dlaset('U', NR - 1, NR - 1, ZERO, ZERO, V(1, 2), LDV);
      }
      dgeqrf(N, NR, V, LDV, WORK(N + 1), WORK(2 * N + 1), LWORK - 2 * N, IERR);
      dlacpy('L', N, NR, V, LDV, WORK(2 * N + 1).asMatrix(N), N);

      for (p = 1; p <= NR; p++) {
        dcopy(NR - p + 1, V(p, p).asArray(), LDV, U(p, p).asArray(), 1);
      }

      if (L2PERT) {
        XSC.value = sqrt(SMALL / EPSLN);
        for (q = 2; q <= NR; q++) {
          for (p = 1; p <= q - 1; p++) {
            TEMP1.value = XSC.value * min(U[p][p].abs(), U[q][q].abs());
            U[p][q] = -sign(TEMP1.value, U[q][p]);
          }
        }
      } else {
        dlaset('U', NR - 1, NR - 1, ZERO, ZERO, U(1, 2), LDU);
      }
      dgesvj('G', 'U', 'V', NR, NR, U, LDU, SVA, N, V, LDV,
          WORK(2 * N + N * NR + 1), LWORK - 2 * N - N * NR, INFO);
      SCALEM = WORK[2 * N + N * NR + 1];
      NUMRANK = WORK[2 * N + N * NR + 2].round();

      if (NR < N) {
        dlaset('A', N - NR, NR, ZERO, ZERO, V(NR + 1, 1), LDV);
        dlaset('A', NR, N - NR, ZERO, ZERO, V(1, NR + 1), LDV);
        dlaset('A', N - NR, N - NR, ZERO, ONE, V(NR + 1, NR + 1), LDV);
      }
      dormqr(
          'L',
          'N',
          N,
          N,
          NR,
          WORK(2 * N + 1).asMatrix(N),
          N,
          WORK(N + 1),
          V,
          LDV,
          WORK(2 * N + N * NR + NR + 1),
          LWORK - 2 * N - N * NR - NR,
          IERR);

      // Permute the rows of V using the (column) permutation from the
      // first QRF. Also, scale the columns to make them unit in
      // Euclidean norm. This applies to all cases.

      TEMP1.value = sqrt(N) * EPSLN;
      for (q = 1; q <= N; q++) {
        for (p = 1; p <= N; p++) {
          WORK[2 * N + N * NR + NR + IWORK[p]] = V[p][q];
        }
        for (p = 1; p <= N; p++) {
          V[p][q] = WORK[2 * N + N * NR + NR + p];
        }
        XSC.value = ONE / dnrm2(N, V(1, q).asArray(), 1);
        if ((XSC.value < (ONE - TEMP1.value)) ||
            (XSC.value > (ONE + TEMP1.value))) {
          dscal(N, XSC.value, V(1, q).asArray(), 1);
        }
      }

      // At this moment, V contains the right singular vectors of A.
      // Next, assemble the left singular vector matrix U (M x N).

      if (NR < M) {
        dlaset('A', M - NR, NR, ZERO, ZERO, U(NR + 1, 1), LDU);
        if (NR < N1) {
          dlaset('A', NR, N1 - NR, ZERO, ZERO, U(1, NR + 1), LDU);
          dlaset('A', M - NR, N1 - NR, ZERO, ONE, U(NR + 1, NR + 1), LDU);
        }
      }

      dormqr('Left', 'No Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N + 1),
          LWORK - N, IERR);

      if (ROWPIV) dlaswp(N1, U, LDU, 1, M - 1, IWORK(2 * N + 1), -1);
    }
    if (TRANSP) {
      // .. swap U and V because the procedure worked on A^t
      for (p = 1; p <= N; p++) {
        dswap(N, U(1, p).asArray(), 1, V(1, p).asArray(), 1);
      }
    }
  }
  // end of the full SVD

  // Undo scaling, if necessary (and possible)

  if (USCAL2 <= (BIG / SVA[1]) * USCAL1) {
    dlascl('G', 0, 0, USCAL1, USCAL2, NR, 1, SVA.asMatrix(N), N, IERR);
    USCAL1 = ONE;
    USCAL2 = ONE;
  }

  if (NR < N) {
    for (p = NR + 1; p <= N; p++) {
      SVA[p] = ZERO;
    }
  }

  WORK[1] = USCAL2 * SCALEM;
  WORK[2] = USCAL1;
  if (ERREST) WORK[3] = SCONDA;
  if (LSVEC && RSVEC) {
    WORK[4] = CONDR1;
    WORK[5] = CONDR2;
  }
  if (L2TRAN) {
    WORK[6] = ENTRA;
    WORK[7] = ENTRAT;
  }

  IWORK[1] = NR;
  IWORK[2] = NUMRANK;
  IWORK[3] = WARNING;
}
