import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/drotm.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgsvj0.dart';
import 'package:lapack/src/dgsvj1.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgesvj(
  final String JOBA,
  final String JOBU,
  final String JOBV,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> SVA_,
  final int MV,
  final Matrix<double> V_,
  final int LDV,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final SVA = SVA_.having();
  final V = V_.having(ld: LDV);
  final WORK = WORK_.having();
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  const NSWEEP = 30;
  double AAPP0,
      AAPQ,
      APOAQ,
      AQOAP,
      BIG,
      BIGTHETA,
      CS,
      CTOL,
      EPSLN,
      // LARGE,
      MXAAPQ = 0,
      MXSINJ = 0,
      ROOTBIG,
      ROOTEPS,
      ROOTSFMIN,
      ROOTTOL,
      SKL,
      SFMIN,
      SMALL,
      SN,
      THETA,
      THSIGN,
      TOL;
  int BLSKIP,
      EMPTSW,
      i,
      ibr,
      igl = 0,
      IJBLSK,
      ir1,
      ISWROT = 0,
      jbc,
      jgl,
      KBL,
      LKAHEAD,
      MVL = 0,
      N2,
      N34,
      N4,
      NBL,
      NOTROT,
      p,
      PSKIPPED,
      q,
      ROWSKIP,
      SWBAND,
      MINMN,
      LWMIN;
  bool APPLV,
      GOSCALE,
      LOWER,
      LQUERY,
      LSVEC,
      NOSCALE,
      ROTOK,
      RSVEC,
      UCTOL,
      UPPER;
  final FASTR = Array<double>(5);
  final IERR = Box(0);
  final AAPP = Box(0.0), AAQQ = Box(0.0), TEMP1 = Box(0.0), T = Box(0.0);

  // Test the input arguments

  LSVEC = lsame(JOBU, 'U');
  UCTOL = lsame(JOBU, 'C');
  RSVEC = lsame(JOBV, 'V');
  APPLV = lsame(JOBV, 'A');
  UPPER = lsame(JOBA, 'U');
  LOWER = lsame(JOBA, 'L');

  MINMN = min(M, N);
  if (MINMN == 0) {
    LWMIN = 1;
  } else {
    LWMIN = max(6, M + N);
  }

  LQUERY = (LWORK == -1);
  if (!(UPPER || LOWER || lsame(JOBA, 'G'))) {
    INFO.value = -1;
  } else if (!(LSVEC || UCTOL || lsame(JOBU, 'N'))) {
    INFO.value = -2;
  } else if (!(RSVEC || APPLV || lsame(JOBV, 'N'))) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if ((N < 0) || (N > M)) {
    INFO.value = -5;
  } else if (LDA < M) {
    INFO.value = -7;
  } else if (MV < 0) {
    INFO.value = -9;
  } else if ((RSVEC && (LDV < N)) || (APPLV && (LDV < MV))) {
    INFO.value = -11;
  } else if (UCTOL && (WORK[1] <= ONE)) {
    INFO.value = -12;
  } else if (LWORK < LWMIN && (!LQUERY)) {
    INFO.value = -13;
  } else {
    INFO.value = 0;
  }

  // #:(
  if (INFO.value != 0) {
    xerbla('DGESVJ', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = LWMIN.toDouble();
    return;
  }

  // #:) Quick return for void matrix

  if (MINMN == 0) return;

  // Set numerical parameters
  // The stopping criterion for Jacobi rotations is

  // max_{i<>j}|A[:][i]^T.value * A[:][j]|/(||A[:][i]||*||A[:][j]||) < CTOL*EPS

  // where EPS is the round-off and CTOL is defined as follows:

  if (UCTOL) {
    // ... user controlled
    CTOL = WORK[1];
  } else {
    // ... default
    if (LSVEC || RSVEC || APPLV) {
      CTOL = sqrt(M.toDouble());
    } else {
      CTOL = M.toDouble();
    }
  }
  // ... and the machine dependent parameters are
  // [!]  (Make sure that dlamch() works properly on the target machine.)

  EPSLN = dlamch('Epsilon');
  ROOTEPS = sqrt(EPSLN);
  SFMIN = dlamch('SafeMinimum');
  ROOTSFMIN = sqrt(SFMIN);
  SMALL = SFMIN / EPSLN;
  BIG = dlamch('Overflow');
  // BIG         = ONE    / SFMIN
  ROOTBIG = ONE / ROOTSFMIN;
  // LARGE = BIG / sqrt((M * N).toDouble());
  BIGTHETA = ONE / ROOTEPS;

  TOL = CTOL * EPSLN;
  ROOTTOL = sqrt(TOL);

  if (M.toDouble() * EPSLN >= ONE) {
    INFO.value = -4;
    xerbla('DGESVJ', -INFO.value);
    return;
  }

  // Initialize the right singular vector matrix.

  if (RSVEC) {
    MVL = N;
    dlaset('A', MVL, N, ZERO, ONE, V, LDV);
  } else if (APPLV) {
    MVL = MV;
  }
  RSVEC = RSVEC || APPLV;

  // Initialize SVA[ 1:N ] = ( ||A e_i||_2, i = 1:N )
  // (!)  If necessary, scale A to protect the largest singular value
  // from overflow. It is possible that saving the largest singular
  // value destroys the information about the small ones.
  // This initial scaling is almost minimal in the sense that the
  // goal is to make sure that no column norm overflows, and that
  // sqrt(N)*max_i SVA[i] does not overflow. If INFinite entries
  // in A are detected, the procedure returns with INFO.value=-6.

  SKL = ONE / sqrt(M.toDouble() * N.toDouble());
  NOSCALE = true;
  GOSCALE = true;

  if (LOWER) {
    // the input matrix is M-by-N lower triangular (trapezoidal)
    for (p = 1; p <= N; p++) {
      AAPP.value = ZERO;
      AAQQ.value = ONE;
      dlassq(M - p + 1, A(p, p).asArray(), 1, AAPP, AAQQ);
      if (AAPP.value > BIG) {
        INFO.value = -6;
        xerbla('DGESVJ', -INFO.value);
        return;
      }
      AAQQ.value = sqrt(AAQQ.value);
      if ((AAPP.value < (BIG / AAQQ.value)) && NOSCALE) {
        SVA[p] = AAPP.value * AAQQ.value;
      } else {
        NOSCALE = false;
        SVA[p] = AAPP.value * (AAQQ.value * SKL);
        if (GOSCALE) {
          GOSCALE = false;
          for (q = 1; q <= p - 1; q++) {
            SVA[q] *= SKL;
          }
        }
      }
    }
  } else if (UPPER) {
    // the input matrix is M-by-N upper triangular (trapezoidal)
    for (p = 1; p <= N; p++) {
      AAPP.value = ZERO;
      AAQQ.value = ONE;
      dlassq(p, A(1, p).asArray(), 1, AAPP, AAQQ);
      if (AAPP.value > BIG) {
        INFO.value = -6;
        xerbla('DGESVJ', -INFO.value);
        return;
      }
      AAQQ.value = sqrt(AAQQ.value);
      if ((AAPP.value < (BIG / AAQQ.value)) && NOSCALE) {
        SVA[p] = AAPP.value * AAQQ.value;
      } else {
        NOSCALE = false;
        SVA[p] = AAPP.value * (AAQQ.value * SKL);
        if (GOSCALE) {
          GOSCALE = false;
          for (q = 1; q <= p - 1; q++) {
            SVA[q] *= SKL;
          }
        }
      }
    }
  } else {
    // the input matrix is M-by-N general dense
    for (p = 1; p <= N; p++) {
      AAPP.value = ZERO;
      AAQQ.value = ONE;
      dlassq(M, A(1, p).asArray(), 1, AAPP, AAQQ);
      if (AAPP.value > BIG) {
        INFO.value = -6;
        xerbla('DGESVJ', -INFO.value);
        return;
      }
      AAQQ.value = sqrt(AAQQ.value);
      if ((AAPP.value < (BIG / AAQQ.value)) && NOSCALE) {
        SVA[p] = AAPP.value * AAQQ.value;
      } else {
        NOSCALE = false;
        SVA[p] = AAPP.value * (AAQQ.value * SKL);
        if (GOSCALE) {
          GOSCALE = false;
          for (q = 1; q <= p - 1; q++) {
            SVA[q] *= SKL;
          }
        }
      }
    }
  }

  if (NOSCALE) SKL = ONE;

  // Move the smaller part of the spectrum from the underflow threshold
  // (!)  Start by determining the position of the nonzero entries of the
  // array SVA() relative to ( SFMIN, BIG ).

  AAPP.value = ZERO;
  AAQQ.value = BIG;
  for (p = 1; p <= N; p++) {
    if (SVA[p] != ZERO) AAQQ.value = min(AAQQ.value, SVA[p]);
    AAPP.value = max(AAPP.value, SVA[p]);
  }

  // #:) Quick return for zero matrix

  if (AAPP.value == ZERO) {
    if (LSVEC) dlaset('G', M, N, ZERO, ONE, A, LDA);
    WORK[1] = ONE;
    WORK[2] = ZERO;
    WORK[3] = ZERO;
    WORK[4] = ZERO;
    WORK[5] = ZERO;
    WORK[6] = ZERO;
    return;
  }

  // #:) Quick return for one-column matrix

  if (N == 1) {
    if (LSVEC) dlascl('G', 0, 0, SVA[1], SKL, M, 1, A(1, 1), LDA, IERR);
    WORK[1] = ONE / SKL;
    if (SVA[1] >= SFMIN) {
      WORK[2] = ONE;
    } else {
      WORK[2] = ZERO;
    }
    WORK[3] = ZERO;
    WORK[4] = ZERO;
    WORK[5] = ZERO;
    WORK[6] = ZERO;
    return;
  }

  // Protect small singular values from underflow, and try to
  // avoid underflows/overflows in computing Jacobi rotations.

  SN = sqrt(SFMIN / EPSLN);
  TEMP1.value = sqrt(BIG / N.toDouble());
  if ((AAPP.value <= SN) ||
      (AAQQ.value >= TEMP1.value) ||
      ((SN <= AAQQ.value) && (AAPP.value <= TEMP1.value))) {
    TEMP1.value = min(BIG, TEMP1.value / AAPP.value);
    // AAQQ.value *=TEMP1.value
    // AAPP.value *=TEMP1.value
  } else if ((AAQQ.value <= SN) && (AAPP.value <= TEMP1.value)) {
    TEMP1.value = min(SN / AAQQ.value, BIG / (AAPP.value * sqrt(N.toDouble())));
    // AAQQ.value *=TEMP1.value
    // AAPP.value *=TEMP1.value
  } else if ((AAQQ.value >= SN) && (AAPP.value >= TEMP1.value)) {
    TEMP1.value = max(SN / AAQQ.value, TEMP1.value / AAPP.value);
    // AAQQ.value *=TEMP1.value
    // AAPP.value *=TEMP1.value
  } else if ((AAQQ.value <= SN) && (AAPP.value >= TEMP1.value)) {
    TEMP1.value = min(SN / AAQQ.value, BIG / (sqrt(N.toDouble()) * AAPP.value));
    // AAQQ.value *=TEMP1.value
    // AAPP.value *=TEMP1.value
  } else {
    TEMP1.value = ONE;
  }

  // Scale, if necessary

  if (TEMP1.value != ONE) {
    dlascl('G', 0, 0, ONE, TEMP1.value, N, 1, SVA.asMatrix(N), N, IERR);
  }
  SKL = TEMP1.value * SKL;
  if (SKL != ONE) {
    dlascl(JOBA, 0, 0, ONE, SKL, M, N, A, LDA, IERR);
    SKL = ONE / SKL;
  }

  // Row-cyclic Jacobi SVD algorithm with column pivoting

  EMPTSW = (N * (N - 1)) ~/ 2;
  NOTROT = 0;
  FASTR[1] = ZERO;

  // A is represented in factored form A *= diag(WORK), where diag(WORK)
  // is initialized to identity. WORK is updated during fast scaled
  // rotations.

  for (q = 1; q <= N; q++) {
    WORK[q] = ONE;
  }

  SWBAND = 3;
  // [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
  // if DGESVJ is used as a computational routine in the preconditioned
  // Jacobi SVD algorithm DGESVJ. For sweeps i=1:SWBAND the procedure
  // works on pivots inside a band-like region around the diagonal.
  // The boundaries are determined dynamically, based on the number of
  // pivots above a threshold.

  KBL = min(8, N);
  // [TP] KBL is a tuning parameter that defines the tile size in the
  // tiling of the p-q loops of pivot pairs. In general, an optimal
  // value of KBL depends on the matrix dimensions and on the
  // parameters of the computer's memory.

  NBL = N ~/ KBL;
  if ((NBL * KBL) != N) NBL++;

  BLSKIP = pow(KBL, 2).toInt();
  // [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

  ROWSKIP = min(5, KBL);
  // [TP] ROWSKIP is a tuning parameter.

  LKAHEAD = 1;
  // [TP] LKAHEAD is a tuning parameter.

  // Quasi block transformations, using the lower (upper) triangular
  // structure of the input matrix. The quasi-block-cycling usually
  // invokes cubic convergence. Big part of this cycle is done inside
  // canonical subspaces of dimensions less than M.

  if ((LOWER || UPPER) && (N > max(64, 4 * KBL))) {
    // [TP] The number of partition levels and the actual partition are
    // tuning parameters.
    N4 = N ~/ 4;
    N2 = N ~/ 2;
    N34 = 3 * N4;
    if (APPLV) {
      q = 0;
    } else {
      q = 1;
    }

    if (LOWER) {
      // This works very well on lower triangular matrices, in particular
      // in the framework of the preconditioned Jacobi SVD (xGEJSV).
      // The idea is simple:
      // [+ 0 0 0]   Note that Jacobi transformations of [0 0]
      // [+ + 0 0]                                       [0 0]
      // [+ + x 0]   actually work on [x 0]              [x 0]
      // [+ + x x]                    [x x].             [x x]

      dgsvj0(
          JOBV,
          M - N34,
          N - N34,
          A(N34 + 1, N34 + 1),
          LDA,
          WORK(N34 + 1),
          SVA(N34 + 1),
          MVL,
          V(N34 * q + 1, N34 + 1),
          LDV,
          EPSLN,
          SFMIN,
          TOL,
          2,
          WORK(N + 1),
          LWORK - N,
          IERR);

      dgsvj0(
          JOBV,
          M - N2,
          N34 - N2,
          A(N2 + 1, N2 + 1),
          LDA,
          WORK(N2 + 1),
          SVA(N2 + 1),
          MVL,
          V(N2 * q + 1, N2 + 1),
          LDV,
          EPSLN,
          SFMIN,
          TOL,
          2,
          WORK(N + 1),
          LWORK - N,
          IERR);

      dgsvj1(
          JOBV,
          M - N2,
          N - N2,
          N4,
          A(N2 + 1, N2 + 1),
          LDA,
          WORK(N2 + 1),
          SVA(N2 + 1),
          MVL,
          V(N2 * q + 1, N2 + 1),
          LDV,
          EPSLN,
          SFMIN,
          TOL,
          1,
          WORK(N + 1),
          LWORK - N,
          IERR);

      dgsvj0(
          JOBV,
          M - N4,
          N2 - N4,
          A(N4 + 1, N4 + 1),
          LDA,
          WORK(N4 + 1),
          SVA(N4 + 1),
          MVL,
          V(N4 * q + 1, N4 + 1),
          LDV,
          EPSLN,
          SFMIN,
          TOL,
          1,
          WORK(N + 1),
          LWORK - N,
          IERR);

      dgsvj0(JOBV, M, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1,
          WORK(N + 1), LWORK - N, IERR);

      dgsvj1(JOBV, M, N2, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL,
          1, WORK(N + 1), LWORK - N, IERR);
    } else if (UPPER) {
      dgsvj0(JOBV, N4, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 2,
          WORK(N + 1), LWORK - N, IERR);

      dgsvj0(
          JOBV,
          N2,
          N4,
          A(1, N4 + 1),
          LDA,
          WORK(N4 + 1),
          SVA(N4 + 1),
          MVL,
          V(N4 * q + 1, N4 + 1),
          LDV,
          EPSLN,
          SFMIN,
          TOL,
          1,
          WORK(N + 1),
          LWORK - N,
          IERR);

      dgsvj1(JOBV, N2, N2, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN,
          TOL, 1, WORK(N + 1), LWORK - N, IERR);

      dgsvj0(
          JOBV,
          N2 + N4,
          N4,
          A(1, N2 + 1),
          LDA,
          WORK(N2 + 1),
          SVA(N2 + 1),
          MVL,
          V(N2 * q + 1, N2 + 1),
          LDV,
          EPSLN,
          SFMIN,
          TOL,
          1,
          WORK(N + 1),
          LWORK - N,
          IERR);
    }
  }

  // .. Row-cyclic pivot strategy with de Rijk's pivoting ..

  for (i = 1; i <= NSWEEP; i++) {
    // .. go go go ...

    MXAAPQ = ZERO;
    MXSINJ = ZERO;
    ISWROT = 0;

    NOTROT = 0;
    PSKIPPED = 0;

    // Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
    // 1 <= p < q <= N. This is the first step toward a blocked implementation
    // of the rotations. New implementation, based on block transformations,
    // is under development.

    for (ibr = 1; ibr <= NBL; ibr++) {
      igl = (ibr - 1) * KBL + 1;

      for (ir1 = 0; ir1 <= min(LKAHEAD, NBL - ibr); ir1++) {
        igl += ir1 * KBL;

        for (p = igl; p <= min(igl + KBL - 1, N - 1); p++) {
          // .. de Rijk's pivoting

          q = idamax(N - p + 1, SVA(p), 1) + p - 1;
          if (p != q) {
            dswap(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1);
            if (RSVEC) dswap(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1);
            TEMP1.value = SVA[p];
            SVA[p] = SVA[q];
            SVA[q] = TEMP1.value;
            TEMP1.value = WORK[p];
            WORK[p] = WORK[q];
            WORK[q] = TEMP1.value;
          }

          if (ir1 == 0) {
            // Column norms are periodically updated by explicit
            // norm computation.
            // Caveat:
            // Unfortunately, some BLAS implementations compute dnrm2(M,A[1][p],1)
            // as sqrt(ddot(M,A[1][p],1,A[1][p],1)), which may cause the result to
            // overflow for ||A[:][p]||_2 > sqrt(overflow_threshold), and to
            // underflow for ||A[:][p]||_2 < sqrt(underflow_threshold).
            // Hence, DNRM2 cannot be trusted, not even in the case when
            // the true norm is far from the under(over)flow boundaries.
            // If properly implemented DNRM2 is available, the if-THEN-ELSE
            // below should read "AAPP.value = dnrm2( M, A[1][p], 1 ) * WORK(p)".

            if ((SVA[p] < ROOTBIG) && (SVA[p] > ROOTSFMIN)) {
              SVA[p] = dnrm2(M, A(1, p).asArray(), 1) * WORK[p];
            } else {
              TEMP1.value = ZERO;
              AAPP.value = ONE;
              dlassq(M, A(1, p).asArray(), 1, TEMP1, AAPP);
              SVA[p] = TEMP1.value * sqrt(AAPP.value) * WORK[p];
            }
            AAPP.value = SVA[p];
          } else {
            AAPP.value = SVA[p];
          }

          if (AAPP.value > ZERO) {
            PSKIPPED = 0;

            for (q = p + 1; q <= min(igl + KBL - 1, N); q++) {
              AAQQ.value = SVA[q];

              if (AAQQ.value > ZERO) {
                AAPP0 = AAPP.value;
                if (AAQQ.value >= ONE) {
                  ROTOK = (SMALL * AAPP.value) <= AAQQ.value;
                  if (AAPP.value < (BIG / AAQQ.value)) {
                    AAPQ =
                        (ddot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) *
                                WORK[p] *
                                WORK[q] /
                                AAQQ.value) /
                            AAPP.value;
                  } else {
                    dcopy(M, A(1, p).asArray(), 1, WORK(N + 1), 1);
                    dlascl('G', 0, 0, AAPP.value, WORK[p], M, 1,
                        WORK(N + 1).asMatrix(LDA), LDA, IERR);
                    AAPQ = ddot(M, WORK(N + 1), 1, A(1, q).asArray(), 1) *
                        WORK[q] /
                        AAQQ.value;
                  }
                } else {
                  ROTOK = AAPP.value <= (AAQQ.value / SMALL);
                  if (AAPP.value > (SMALL / AAQQ.value)) {
                    AAPQ =
                        (ddot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) *
                                WORK[p] *
                                WORK[q] /
                                AAQQ.value) /
                            AAPP.value;
                  } else {
                    dcopy(M, A(1, q).asArray(), 1, WORK(N + 1), 1);
                    dlascl('G', 0, 0, AAQQ.value, WORK[q], M, 1,
                        WORK(N + 1).asMatrix(LDA), LDA, IERR);
                    AAPQ = ddot(M, WORK(N + 1), 1, A(1, p).asArray(), 1) *
                        WORK[p] /
                        AAPP.value;
                  }
                }

                MXAAPQ = max(MXAAPQ, AAPQ.abs());

                // TO rotate or NOT to rotate, THAT is the question ...

                if (AAPQ.abs() > TOL) {
                  // .. rotate
                  // [RTD]      ROTATED += ONE

                  if (ir1 == 0) {
                    NOTROT = 0;
                    PSKIPPED = 0;
                    ISWROT++;
                  }

                  if (ROTOK) {
                    AQOAP = AAQQ.value / AAPP.value;
                    APOAQ = AAPP.value / AAQQ.value;
                    THETA = -HALF * (AQOAP - APOAQ).abs() / AAPQ;

                    if (THETA.abs() > BIGTHETA) {
                      T.value = HALF / THETA;
                      FASTR[3] = T.value * WORK[p] / WORK[q];
                      FASTR[4] = -T.value * WORK[q] / WORK[p];
                      drotm(
                          M, A(1, p).asArray(), 1, A(1, q).asArray(), 1, FASTR);
                      if (RSVEC) {
                        drotm(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1,
                            FASTR);
                      }
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ));
                      AAPP.value *=
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ));
                      MXSINJ = max(MXSINJ, T.value.abs());
                    } else {
                      // .. choose correct signum for THETA and rotate

                      THSIGN = -sign(ONE, AAPQ).toDouble();
                      T.value =
                          ONE / (THETA + THSIGN * sqrt(ONE + THETA * THETA));
                      CS = sqrt(ONE / (ONE + T.value * T.value));
                      SN = T.value * CS;

                      MXSINJ = max(MXSINJ, SN.abs());
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ));
                      AAPP.value *=
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ));

                      APOAQ = WORK[p] / WORK[q];
                      AQOAP = WORK[q] / WORK[p];
                      if (WORK[p] >= ONE) {
                        if (WORK[q] >= ONE) {
                          FASTR[3] = T.value * APOAQ;
                          FASTR[4] = -T.value * AQOAP;
                          WORK[p] *= CS;
                          WORK[q] *= CS;
                          drotm(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1,
                              FASTR);
                          if (RSVEC) {
                            drotm(MVL, V(1, p).asArray(), 1, V(1, q).asArray(),
                                1, FASTR);
                          }
                        } else {
                          daxpy(M, -T.value * AQOAP, A(1, q).asArray(), 1,
                              A(1, p).asArray(), 1);
                          daxpy(M, CS * SN * APOAQ, A(1, p).asArray(), 1,
                              A(1, q).asArray(), 1);
                          WORK[p] *= CS;
                          WORK[q] /= CS;
                          if (RSVEC) {
                            daxpy(MVL, -T.value * AQOAP, V(1, q).asArray(), 1,
                                V(1, p).asArray(), 1);
                            daxpy(MVL, CS * SN * APOAQ, V(1, p).asArray(), 1,
                                V(1, q).asArray(), 1);
                          }
                        }
                      } else {
                        if (WORK[q] >= ONE) {
                          daxpy(M, T.value * APOAQ, A(1, p).asArray(), 1,
                              A(1, q).asArray(), 1);
                          daxpy(M, -CS * SN * AQOAP, A(1, q).asArray(), 1,
                              A(1, p).asArray(), 1);
                          WORK[p] /= CS;
                          WORK[q] *= CS;
                          if (RSVEC) {
                            daxpy(MVL, T.value * APOAQ, V(1, p).asArray(), 1,
                                V(1, q).asArray(), 1);
                            daxpy(MVL, -CS * SN * AQOAP, V(1, q).asArray(), 1,
                                V(1, p).asArray(), 1);
                          }
                        } else {
                          if (WORK[p] >= WORK[q]) {
                            daxpy(M, -T.value * AQOAP, A(1, q).asArray(), 1,
                                A(1, p).asArray(), 1);
                            daxpy(M, CS * SN * APOAQ, A(1, p).asArray(), 1,
                                A(1, q).asArray(), 1);
                            WORK[p] *= CS;
                            WORK[q] /= CS;
                            if (RSVEC) {
                              daxpy(MVL, -T.value * AQOAP, V(1, q).asArray(), 1,
                                  V(1, p).asArray(), 1);
                              daxpy(MVL, CS * SN * APOAQ, V(1, p).asArray(), 1,
                                  V(1, q).asArray(), 1);
                            }
                          } else {
                            daxpy(M, T.value * APOAQ, A(1, p).asArray(), 1,
                                A(1, q).asArray(), 1);
                            daxpy(M, -CS * SN * AQOAP, A(1, q).asArray(), 1,
                                A(1, p).asArray(), 1);
                            WORK[p] /= CS;
                            WORK[q] *= CS;
                            if (RSVEC) {
                              daxpy(MVL, T.value * APOAQ, V(1, p).asArray(), 1,
                                  V(1, q).asArray(), 1);
                              daxpy(MVL, -CS * SN * AQOAP, V(1, q).asArray(), 1,
                                  V(1, p).asArray(), 1);
                            }
                          }
                        }
                      }
                    }
                  } else {
                    // .. have to use modified Gram-Schmidt like transformation
                    dcopy(M, A(1, p).asArray(), 1, WORK(N + 1), 1);
                    dlascl('G', 0, 0, AAPP.value, ONE, M, 1,
                        WORK(N + 1).asMatrix(LDA), LDA, IERR);
                    dlascl(
                        'G', 0, 0, AAQQ.value, ONE, M, 1, A(1, q), LDA, IERR);
                    TEMP1.value = -AAPQ * WORK[p] / WORK[q];
                    daxpy(M, TEMP1.value, WORK(N + 1), 1, A(1, q).asArray(), 1);
                    dlascl(
                        'G', 0, 0, ONE, AAQQ.value, M, 1, A(1, q), LDA, IERR);
                    SVA[q] = AAQQ.value * sqrt(max(ZERO, ONE - AAPQ * AAPQ));
                    MXSINJ = max(MXSINJ, SFMIN);
                  }
                  // END if ROTOK THEN ... ELSE

                  // In the case of cancellation in updating SVA[q], SVA[p]
                  // recompute SVA[q], SVA[p].

                  if (pow((SVA[q] / AAQQ.value), 2) <= ROOTEPS) {
                    if ((AAQQ.value < ROOTBIG) && (AAQQ.value > ROOTSFMIN)) {
                      SVA[q] = dnrm2(M, A(1, q).asArray(), 1) * WORK[q];
                    } else {
                      T.value = ZERO;
                      AAQQ.value = ONE;
                      dlassq(M, A(1, q).asArray(), 1, T, AAQQ);
                      SVA[q] = T.value * sqrt(AAQQ.value) * WORK[q];
                    }
                  }
                  if ((AAPP.value / AAPP0) <= ROOTEPS) {
                    if ((AAPP.value < ROOTBIG) && (AAPP.value > ROOTSFMIN)) {
                      AAPP.value = dnrm2(M, A(1, p).asArray(), 1) * WORK[p];
                    } else {
                      T.value = ZERO;
                      AAPP.value = ONE;
                      dlassq(M, A(1, p).asArray(), 1, T, AAPP);
                      AAPP.value = T.value * sqrt(AAPP.value) * WORK[p];
                    }
                    SVA[p] = AAPP.value;
                  }
                } else {
                  // A[:][p] and A[:][q] already numerically orthogonal
                  if (ir1 == 0) NOTROT++;
                  // [RTD]      SKIPPED++
                  PSKIPPED++;
                }
              } else {
                // A[:][q] is zero column
                if (ir1 == 0) NOTROT++;
                PSKIPPED++;
              }

              if ((i <= SWBAND) && (PSKIPPED > ROWSKIP)) {
                if (ir1 == 0) AAPP.value = -AAPP.value;
                NOTROT = 0;
                break;
              }
            }
            // END q-LOOP

            SVA[p] = AAPP.value;
          } else {
            SVA[p] = AAPP.value;
            if ((ir1 == 0) && (AAPP.value == ZERO)) {
              NOTROT += min(igl + KBL - 1, N).toInt() - p;
            }
          }
        }
        // end of the p-loop
        // end of doing the block ( ibr, ibr )
      }
      // end of ir1-loop

      // ... gotothe off diagonal blocks

      igl = (ibr - 1) * KBL + 1;

      jbcLoop:
      for (jbc = ibr + 1; jbc <= NBL; jbc++) {
        jgl = (jbc - 1) * KBL + 1;

        // doing the block at ( ibr, jbc )

        IJBLSK = 0;
        for (p = igl; p <= min(igl + KBL - 1, N); p++) {
          AAPP.value = SVA[p];
          if (AAPP.value > ZERO) {
            PSKIPPED = 0;

            for (q = jgl; q <= min(jgl + KBL - 1, N); q++) {
              AAQQ.value = SVA[q];
              if (AAQQ.value > ZERO) {
                AAPP0 = AAPP.value;

                // .. M x 2 Jacobi SVD ..

                // Safe Gram matrix computation

                if (AAQQ.value >= ONE) {
                  if (AAPP.value >= AAQQ.value) {
                    ROTOK = (SMALL * AAPP.value) <= AAQQ.value;
                  } else {
                    ROTOK = (SMALL * AAQQ.value) <= AAPP.value;
                  }
                  if (AAPP.value < (BIG / AAQQ.value)) {
                    AAPQ =
                        (ddot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) *
                                WORK[p] *
                                WORK[q] /
                                AAQQ.value) /
                            AAPP.value;
                  } else {
                    dcopy(M, A(1, p).asArray(), 1, WORK(N + 1), 1);
                    dlascl('G', 0, 0, AAPP.value, WORK[p], M, 1,
                        WORK(N + 1).asMatrix(LDA), LDA, IERR);
                    AAPQ = ddot(M, WORK(N + 1), 1, A(1, q).asArray(), 1) *
                        WORK[q] /
                        AAQQ.value;
                  }
                } else {
                  if (AAPP.value >= AAQQ.value) {
                    ROTOK = AAPP.value <= (AAQQ.value / SMALL);
                  } else {
                    ROTOK = AAQQ.value <= (AAPP.value / SMALL);
                  }
                  if (AAPP.value > (SMALL / AAQQ.value)) {
                    AAPQ =
                        (ddot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) *
                                WORK[p] *
                                WORK[q] /
                                AAQQ.value) /
                            AAPP.value;
                  } else {
                    dcopy(M, A(1, q).asArray(), 1, WORK(N + 1), 1);
                    dlascl('G', 0, 0, AAQQ.value, WORK[q], M, 1,
                        WORK(N + 1).asMatrix(LDA), LDA, IERR);
                    AAPQ = ddot(M, WORK(N + 1), 1, A(1, p).asArray(), 1) *
                        WORK[p] /
                        AAPP.value;
                  }
                }

                MXAAPQ = max(MXAAPQ, AAPQ.abs());

                // TO rotate or NOT to rotate, THAT is the question ...

                if (AAPQ.abs() > TOL) {
                  NOTROT = 0;
                  // [RTD]      ROTATED++
                  PSKIPPED = 0;
                  ISWROT++;

                  if (ROTOK) {
                    AQOAP = AAQQ.value / AAPP.value;
                    APOAQ = AAPP.value / AAQQ.value;
                    THETA = -HALF * (AQOAP - APOAQ).abs() / AAPQ;
                    if (AAQQ.value > AAPP0) THETA = -THETA;

                    if (THETA.abs() > BIGTHETA) {
                      T.value = HALF / THETA;
                      FASTR[3] = T.value * WORK[p] / WORK[q];
                      FASTR[4] = -T.value * WORK[q] / WORK[p];
                      drotm(
                          M, A(1, p).asArray(), 1, A(1, q).asArray(), 1, FASTR);
                      if (RSVEC) {
                        drotm(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1,
                            FASTR);
                      }
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ));
                      AAPP.value *=
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ));
                      MXSINJ = max(MXSINJ, T.value.abs());
                    } else {
                      // .. choose correct signum for THETA and rotate

                      THSIGN = -sign(ONE, AAPQ).toDouble();
                      if (AAQQ.value > AAPP0) THSIGN = -THSIGN;
                      T.value =
                          ONE / (THETA + THSIGN * sqrt(ONE + THETA * THETA));
                      CS = sqrt(ONE / (ONE + T.value * T.value));
                      SN = T.value * CS;
                      MXSINJ = max(MXSINJ, SN.abs());
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ));
                      AAPP.value *=
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ));

                      APOAQ = WORK[p] / WORK[q];
                      AQOAP = WORK[q] / WORK[p];
                      if (WORK[p] >= ONE) {
                        if (WORK[q] >= ONE) {
                          FASTR[3] = T.value * APOAQ;
                          FASTR[4] = -T.value * AQOAP;
                          WORK[p] *= CS;
                          WORK[q] *= CS;
                          drotm(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1,
                              FASTR);
                          if (RSVEC) {
                            drotm(MVL, V(1, p).asArray(), 1, V(1, q).asArray(),
                                1, FASTR);
                          }
                        } else {
                          daxpy(M, -T.value * AQOAP, A(1, q).asArray(), 1,
                              A(1, p).asArray(), 1);
                          daxpy(M, CS * SN * APOAQ, A(1, p).asArray(), 1,
                              A(1, q).asArray(), 1);
                          if (RSVEC) {
                            daxpy(MVL, -T.value * AQOAP, V(1, q).asArray(), 1,
                                V(1, p).asArray(), 1);
                            daxpy(MVL, CS * SN * APOAQ, V(1, p).asArray(), 1,
                                V(1, q).asArray(), 1);
                          }
                          WORK[p] *= CS;
                          WORK[q] /= CS;
                        }
                      } else {
                        if (WORK[q] >= ONE) {
                          daxpy(M, T.value * APOAQ, A(1, p).asArray(), 1,
                              A(1, q).asArray(), 1);
                          daxpy(M, -CS * SN * AQOAP, A(1, q).asArray(), 1,
                              A(1, p).asArray(), 1);
                          if (RSVEC) {
                            daxpy(MVL, T.value * APOAQ, V(1, p).asArray(), 1,
                                V(1, q).asArray(), 1);
                            daxpy(MVL, -CS * SN * AQOAP, V(1, q).asArray(), 1,
                                V(1, p).asArray(), 1);
                          }
                          WORK[p] /= CS;
                          WORK[q] *= CS;
                        } else {
                          if (WORK[p] >= WORK[q]) {
                            daxpy(M, -T.value * AQOAP, A(1, q).asArray(), 1,
                                A(1, p).asArray(), 1);
                            daxpy(M, CS * SN * APOAQ, A(1, p).asArray(), 1,
                                A(1, q).asArray(), 1);
                            WORK[p] *= CS;
                            WORK[q] /= CS;
                            if (RSVEC) {
                              daxpy(MVL, -T.value * AQOAP, V(1, q).asArray(), 1,
                                  V(1, p).asArray(), 1);
                              daxpy(MVL, CS * SN * APOAQ, V(1, p).asArray(), 1,
                                  V(1, q).asArray(), 1);
                            }
                          } else {
                            daxpy(M, T.value * APOAQ, A(1, p).asArray(), 1,
                                A(1, q).asArray(), 1);
                            daxpy(M, -CS * SN * AQOAP, A(1, q).asArray(), 1,
                                A(1, p).asArray(), 1);
                            WORK[p] /= CS;
                            WORK[q] *= CS;
                            if (RSVEC) {
                              daxpy(MVL, T.value * APOAQ, V(1, p).asArray(), 1,
                                  V(1, q).asArray(), 1);
                              daxpy(MVL, -CS * SN * AQOAP, V(1, q).asArray(), 1,
                                  V(1, p).asArray(), 1);
                            }
                          }
                        }
                      }
                    }
                  } else {
                    if (AAPP.value > AAQQ.value) {
                      dcopy(M, A(1, p).asArray(), 1, WORK(N + 1), 1);
                      dlascl('G', 0, 0, AAPP.value, ONE, M, 1,
                          WORK(N + 1).asMatrix(LDA), LDA, IERR);
                      dlascl(
                          'G', 0, 0, AAQQ.value, ONE, M, 1, A(1, q), LDA, IERR);
                      TEMP1.value = -AAPQ * WORK[p] / WORK[q];
                      daxpy(
                          M, TEMP1.value, WORK(N + 1), 1, A(1, q).asArray(), 1);
                      dlascl(
                          'G', 0, 0, ONE, AAQQ.value, M, 1, A(1, q), LDA, IERR);
                      SVA[q] = AAQQ.value * sqrt(max(ZERO, ONE - AAPQ * AAPQ));
                      MXSINJ = max(MXSINJ, SFMIN);
                    } else {
                      dcopy(M, A(1, q).asArray(), 1, WORK(N + 1), 1);
                      dlascl('G', 0, 0, AAQQ.value, ONE, M, 1,
                          WORK(N + 1).asMatrix(LDA), LDA, IERR);
                      dlascl(
                          'G', 0, 0, AAPP.value, ONE, M, 1, A(1, p), LDA, IERR);
                      TEMP1.value = -AAPQ * WORK[q] / WORK[p];
                      daxpy(
                          M, TEMP1.value, WORK(N + 1), 1, A(1, p).asArray(), 1);
                      dlascl(
                          'G', 0, 0, ONE, AAPP.value, M, 1, A(1, p), LDA, IERR);
                      SVA[p] = AAPP.value * sqrt(max(ZERO, ONE - AAPQ * AAPQ));
                      MXSINJ = max(MXSINJ, SFMIN);
                    }
                  }
                  // END if ROTOK THEN ... ELSE

                  // In the case of cancellation in updating SVA[q]
                  // .. recompute SVA[q]
                  if (pow((SVA[q] / AAQQ.value), 2) <= ROOTEPS) {
                    if ((AAQQ.value < ROOTBIG) && (AAQQ.value > ROOTSFMIN)) {
                      SVA[q] = dnrm2(M, A(1, q).asArray(), 1) * WORK[q];
                    } else {
                      T.value = ZERO;
                      AAQQ.value = ONE;
                      dlassq(M, A(1, q).asArray(), 1, T, AAQQ);
                      SVA[q] = T.value * sqrt(AAQQ.value) * WORK[q];
                    }
                  }
                  if (pow((AAPP.value / AAPP0), 2) <= ROOTEPS) {
                    if ((AAPP.value < ROOTBIG) && (AAPP.value > ROOTSFMIN)) {
                      AAPP.value = dnrm2(M, A(1, p).asArray(), 1) * WORK[p];
                    } else {
                      T.value = ZERO;
                      AAPP.value = ONE;
                      dlassq(M, A(1, p).asArray(), 1, T, AAPP);
                      AAPP.value = T.value * sqrt(AAPP.value) * WORK[p];
                    }
                    SVA[p] = AAPP.value;
                  }
                  // end of OK rotation
                } else {
                  NOTROT++;
                  // [RTD]      SKIPPED++
                  PSKIPPED++;
                  IJBLSK++;
                }
              } else {
                NOTROT++;
                PSKIPPED++;
                IJBLSK++;
              }

              if ((i <= SWBAND) && (IJBLSK >= BLSKIP)) {
                SVA[p] = AAPP.value;
                NOTROT = 0;
                break jbcLoop;
              }
              if ((i <= SWBAND) && (PSKIPPED > ROWSKIP)) {
                AAPP.value = -AAPP.value;
                NOTROT = 0;
                break;
              }
            }
            // end of the q-loop

            SVA[p] = AAPP.value;
          } else {
            if (AAPP.value == ZERO) {
              NOTROT += min(jgl + KBL - 1, N).toInt() - jgl + 1;
            }
            if (AAPP.value < ZERO) NOTROT = 0;
          }
        }
        // end of the p-loop
      }
      // end of the jbc-loop
      // }
      // 2011 bailed out of the jbc-loop
      for (p = igl; p <= min(igl + KBL - 1, N); p++) {
        SVA[p] = SVA[p].abs();
      }
      // **
    }
    // 2000 :: end of the ibr-loop

    // .. update SVA[N]
    if ((SVA[N] < ROOTBIG) && (SVA[N] > ROOTSFMIN)) {
      SVA[N] = dnrm2(M, A(1, N).asArray(), 1) * WORK[N];
    } else {
      T.value = ZERO;
      AAPP.value = ONE;
      dlassq(M, A(1, N).asArray(), 1, T, AAPP);
      SVA[N] = T.value * sqrt(AAPP.value) * WORK[N];
    }

    // Additional steering devices

    if ((i < SWBAND) && ((MXAAPQ <= ROOTTOL) || (ISWROT <= N))) SWBAND = i;

    if ((i > SWBAND + 1) &&
            (MXAAPQ < sqrt(N.toDouble()) * TOL) &&
            (N.toDouble() * MXAAPQ * MXSINJ < TOL) ||
        NOTROT >= EMPTSW) {
      // #:) Reaching this point means numerical convergence after the i-th
      // sweep.

      // #:) INFO.value = 0 confirms successful iterations.
      INFO.value = 0;
      continue;
    }

    // #:( Reaching this point means that the procedure has not converged.
    INFO.value = NSWEEP - 1;
  }

  // Sort the singular values and find how many are above
  // the underflow threshold.

  N2 = 0;
  N4 = 0;
  for (p = 1; p <= N - 1; p++) {
    q = idamax(N - p + 1, SVA(p), 1) + p - 1;
    if (p != q) {
      TEMP1.value = SVA[p];
      SVA[p] = SVA[q];
      SVA[q] = TEMP1.value;
      TEMP1.value = WORK[p];
      WORK[p] = WORK[q];
      WORK[q] = TEMP1.value;
      dswap(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1);
      if (RSVEC) dswap(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1);
    }
    if (SVA[p] != ZERO) {
      N4++;
      if (SVA[p] * SKL > SFMIN) N2++;
    }
  }
  if (SVA[N] != ZERO) {
    N4++;
    if (SVA[N] * SKL > SFMIN) N2++;
  }

  // Normalize the left singular vectors.

  if (LSVEC || UCTOL) {
    for (p = 1; p <= N2; p++) {
      dscal(M, WORK[p] / SVA[p], A(1, p).asArray(), 1);
    }
  }

  // Scale the product of Jacobi rotations (assemble the fast rotations).

  if (RSVEC) {
    if (APPLV) {
      for (p = 1; p <= N; p++) {
        dscal(MVL, WORK[p], V(1, p).asArray(), 1);
      }
    } else {
      for (p = 1; p <= N; p++) {
        TEMP1.value = ONE / dnrm2(MVL, V(1, p).asArray(), 1);
        dscal(MVL, TEMP1.value, V(1, p).asArray(), 1);
      }
    }
  }

  // Undo scaling, if necessary (and possible).
  if (((SKL > ONE) && (SVA[1] < (BIG / SKL))) ||
      ((SKL < ONE) && (SVA[max(N2, 1)] > (SFMIN / SKL)))) {
    for (p = 1; p <= N; p++) {
      SVA[p] = SKL * SVA[p];
    }
    SKL = ONE;
  }

  WORK[1] = SKL;
  // The singular values of A are SKL*SVA[1:N]. If SKL != ONE
  // then some of the singular values may overflow or underflow and
  // the spectrum is given in this factored representation.

  WORK[2] = N4.toDouble();
  // N4 is the number of computed nonzero singular values of A.

  WORK[3] = N2.toDouble();
  // N2 is the number of singular values of A greater than SFMIN.
  // If N2<N, SVA[N2:N] contains ZEROS and/or denormalized numbers
  // that may carry some information.

  WORK[4] = i.toDouble();
  // i is the index of the last sweep before declaring convergence.

  WORK[5] = MXAAPQ;
  // MXAAPQ is the largest absolute value of scaled pivots in the
  // last sweep

  WORK[6] = MXSINJ;
  // MXSINJ is the largest absolute value of the sines of Jacobi angles
  // in the last sweep
}
