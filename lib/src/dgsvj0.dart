import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/drotm.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgsvj0(
  final String JOBV,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> SVA_,
  final int MV,
  final Matrix<double> V_,
  final int LDV,
  final double EPS,
  final double SFMIN,
  final double TOL,
  final int NSWEEP,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final SVA = SVA_.having();
  final V = V_.having(ld: LDV);
  final WORK = WORK_.having();
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  double AAPP0,
      AAPQ,
      APOAQ,
      AQOAP,
      BIG,
      BIGTHETA,
      CS,
      MXAAPQ,
      MXSINJ,
      ROOTBIG,
      ROOTEPS,
      ROOTSFMIN,
      ROOTTOL,
      SMALL,
      SN,
      THETA,
      THSIGN;
  int BLSKIP,
      EMPTSW,
      i,
      ibr,
      igl,
      IJBLSK,
      ir1,
      ISWROT,
      jbc,
      jgl,
      KBL,
      LKAHEAD,
      MVL = 0,
      NBL,
      NOTROT,
      p,
      PSKIPPED,
      q,
      ROWSKIP,
      SWBAND;
  bool APPLV, ROTOK, RSVEC;
  final FASTR = Array<double>(5);
  final IERR = Box(0);
  final AAPP = Box(0.0), TEMP1 = Box(0.0), T = Box(0.0), AAQQ = Box(0.0);

  // Test the input parameters.

  APPLV = lsame(JOBV, 'A');
  RSVEC = lsame(JOBV, 'V');
  if (!(RSVEC || APPLV || lsame(JOBV, 'N'))) {
    INFO.value = -1;
  } else if (M < 0) {
    INFO.value = -2;
  } else if ((N < 0) || (N > M)) {
    INFO.value = -3;
  } else if (LDA < M) {
    INFO.value = -5;
  } else if ((RSVEC || APPLV) && (MV < 0)) {
    INFO.value = -8;
  } else if ((RSVEC && (LDV < N)) || (APPLV && (LDV < MV))) {
    INFO.value = -10;
  } else if (TOL <= EPS) {
    INFO.value = -13;
  } else if (NSWEEP < 0) {
    INFO.value = -14;
  } else if (LWORK < M) {
    INFO.value = -16;
  } else {
    INFO.value = 0;
  }

  if (INFO.value != 0) {
    xerbla('DGSVJ0', -INFO.value);
    return;
  }

  if (RSVEC) {
    MVL = N;
  } else if (APPLV) {
    MVL = MV;
  }
  RSVEC = RSVEC || APPLV;

  ROOTEPS = sqrt(EPS);
  ROOTSFMIN = sqrt(SFMIN);
  SMALL = SFMIN / EPS;
  BIG = ONE / SFMIN;
  ROOTBIG = ONE / ROOTSFMIN;
  BIGTHETA = ONE / ROOTEPS;
  ROOTTOL = sqrt(TOL);

  // -#- Row-cyclic Jacobi SVD algorithm with column pivoting -#-

  EMPTSW = (N * (N - 1)) ~/ 2;
  NOTROT = 0;
  FASTR[1] = ZERO;

  // -#- Row-cyclic pivot strategy with de Rijk's pivoting -#-

  SWBAND = 0;
  // [TP] SWBAND is a tuning parameter. It is meaningful and effective
  // if SGESVJ is used as a computational routine in the preconditioned
  // Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure
  // ......

  KBL = min(8, N);
  // [TP] KBL is a tuning parameter that defines the tile size in the
  // tiling of the p-q loops of pivot pairs. In general, an optimal
  // value of KBL depends on the matrix dimensions and on the
  // parameters of the computer's memory.

  NBL = N ~/ KBL;
  if ((NBL * KBL) != N) NBL++;

  BLSKIP = pow(KBL, 2).toInt() + 1;
  // [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

  ROWSKIP = min(5, KBL);
  // [TP] ROWSKIP is a tuning parameter.

  LKAHEAD = 1;
  // [TP] LKAHEAD is a tuning parameter.
  SWBAND = 0;
  PSKIPPED = 0;

  var isBelowTolerance = false;
  for (i = 1; i <= NSWEEP; i++) {
    // .. go go go ...

    MXAAPQ = ZERO;
    MXSINJ = ZERO;
    ISWROT = 0;

    NOTROT = 0;
    PSKIPPED = 0;

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
            TEMP1.value = D[p];
            D[p] = D[q];
            D[q] = TEMP1.value;
          }

          if (ir1 == 0) {
            // Column norms are periodically updated by explicit
            // norm computation.
            // Caveat:
            // Some BLAS implementations compute dnrm2(M,A[1][p],1)
            // as sqrt(ddot(M,A[1][p],1,A[1][p],1)), which may result in
            // overflow for ||A[:][p]||_2 > sqrt(overflow_threshold), and
            // underflow for ||A[:][p]||_2 < sqrt(underflow_threshold).
            // Hence, DNRM2 cannot be trusted, not even in the case when
            // the true norm is far from the under(over)flow boundaries.
            // If properly implemented DNRM2 is available, the if-THEN-ELSE
            // below should read "AAPP.value = dnrm2( M, A[1][p], 1 ) * D[p]".

            if ((SVA[p] < ROOTBIG) && (SVA[p] > ROOTSFMIN)) {
              SVA[p] = dnrm2(M, A(1, p).asArray(), 1) * D[p];
            } else {
              TEMP1.value = ZERO;
              AAPP.value = ONE;
              dlassq(M, A(1, p).asArray(), 1, TEMP1, AAPP);
              SVA[p] = TEMP1.value * sqrt(AAPP.value) * D[p];
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
                                D[p] *
                                D[q] /
                                AAQQ.value) /
                            AAPP.value;
                  } else {
                    dcopy(M, A(1, p).asArray(), 1, WORK, 1);
                    dlascl('G', 0, 0, AAPP.value, D[p], M, 1,
                        WORK.asMatrix(LDA), LDA, IERR);
                    AAPQ = ddot(M, WORK, 1, A(1, q).asArray(), 1) *
                        D[q] /
                        AAQQ.value;
                  }
                } else {
                  ROTOK = AAPP.value <= (AAQQ.value / SMALL);
                  if (AAPP.value > (SMALL / AAQQ.value)) {
                    AAPQ =
                        (ddot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) *
                                D[p] *
                                D[q] /
                                AAQQ.value) /
                            AAPP.value;
                  } else {
                    dcopy(M, A(1, q).asArray(), 1, WORK, 1);
                    dlascl('G', 0, 0, AAQQ.value, D[q], M, 1,
                        WORK.asMatrix(LDA), LDA, IERR);
                    AAPQ = ddot(M, WORK, 1, A(1, p).asArray(), 1) *
                        D[p] /
                        AAPP.value;
                  }
                }

                MXAAPQ = max(MXAAPQ, AAPQ.abs());

                // TO rotate or NOT to rotate, THAT is the question ...

                if (AAPQ.abs() > TOL) {
                  // .. rotate
                  // ROTATED = ROTATED + ONE

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
                      FASTR[3] = T.value * D[p] / D[q];
                      FASTR[4] = -T.value * D[q] / D[p];
                      drotm(
                          M, A(1, p).asArray(), 1, A(1, q).asArray(), 1, FASTR);
                      if (RSVEC) {
                        drotm(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1,
                            FASTR);
                      }
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ));
                      AAPP.value = AAPP.value *
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
                      AAPP.value = AAPP.value *
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ));

                      APOAQ = D[p] / D[q];
                      AQOAP = D[q] / D[p];
                      if (D[p] >= ONE) {
                        if (D[q] >= ONE) {
                          FASTR[3] = T.value * APOAQ;
                          FASTR[4] = -T.value * AQOAP;
                          D[p] = D[p] * CS;
                          D[q] = D[q] * CS;
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
                          D[p] = D[p] * CS;
                          D[q] = D[q] / CS;
                          if (RSVEC) {
                            daxpy(MVL, -T.value * AQOAP, V(1, q).asArray(), 1,
                                V(1, p).asArray(), 1);
                            daxpy(MVL, CS * SN * APOAQ, V(1, p).asArray(), 1,
                                V(1, q).asArray(), 1);
                          }
                        }
                      } else {
                        if (D[q] >= ONE) {
                          daxpy(M, T.value * APOAQ, A(1, p).asArray(), 1,
                              A(1, q).asArray(), 1);
                          daxpy(M, -CS * SN * AQOAP, A(1, q).asArray(), 1,
                              A(1, p).asArray(), 1);
                          D[p] = D[p] / CS;
                          D[q] = D[q] * CS;
                          if (RSVEC) {
                            daxpy(MVL, T.value * APOAQ, V(1, p).asArray(), 1,
                                V(1, q).asArray(), 1);
                            daxpy(MVL, -CS * SN * AQOAP, V(1, q).asArray(), 1,
                                V(1, p).asArray(), 1);
                          }
                        } else {
                          if (D[p] >= D[q]) {
                            daxpy(M, -T.value * AQOAP, A(1, q).asArray(), 1,
                                A(1, p).asArray(), 1);
                            daxpy(M, CS * SN * APOAQ, A(1, p).asArray(), 1,
                                A(1, q).asArray(), 1);
                            D[p] = D[p] * CS;
                            D[q] = D[q] / CS;
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
                            D[p] = D[p] / CS;
                            D[q] = D[q] * CS;
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
                    dcopy(M, A(1, p).asArray(), 1, WORK, 1);
                    dlascl('G', 0, 0, AAPP.value, ONE, M, 1, WORK.asMatrix(LDA),
                        LDA, IERR);
                    dlascl(
                        'G', 0, 0, AAQQ.value, ONE, M, 1, A(1, q), LDA, IERR);
                    TEMP1.value = -AAPQ * D[p] / D[q];
                    daxpy(M, TEMP1.value, WORK, 1, A(1, q).asArray(), 1);
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
                      SVA[q] = dnrm2(M, A(1, q).asArray(), 1) * D[q];
                    } else {
                      T.value = ZERO;
                      AAQQ.value = ONE;
                      dlassq(M, A(1, q).asArray(), 1, T, AAQQ);
                      SVA[q] = T.value * sqrt(AAQQ.value) * D[q];
                    }
                  }
                  if ((AAPP.value / AAPP0) <= ROOTEPS) {
                    if ((AAPP.value < ROOTBIG) && (AAPP.value > ROOTSFMIN)) {
                      AAPP.value = dnrm2(M, A(1, p).asArray(), 1) * D[p];
                    } else {
                      T.value = ZERO;
                      AAPP.value = ONE;
                      dlassq(M, A(1, p).asArray(), 1, T, AAPP);
                      AAPP.value = T.value * sqrt(AAPP.value) * D[p];
                    }
                    SVA[p] = AAPP.value;
                  }
                } else {
                  // A[:][p] and A[:][q] already numerically orthogonal
                  if (ir1 == 0) NOTROT++;
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
                //  GO TO 2103;
                break;
              }
            }
            // END q-LOOP

            //  }
            // bailed out of q-loop

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

      // ........................................................
      // ... go to the off diagonal blocks

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

                // -#- M x 2 Jacobi SVD -#-

                // -#- Safe Gram matrix computation -#-

                if (AAQQ.value >= ONE) {
                  if (AAPP.value >= AAQQ.value) {
                    ROTOK = (SMALL * AAPP.value) <= AAQQ.value;
                  } else {
                    ROTOK = (SMALL * AAQQ.value) <= AAPP.value;
                  }
                  if (AAPP.value < (BIG / AAQQ.value)) {
                    AAPQ =
                        (ddot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) *
                                D[p] *
                                D[q] /
                                AAQQ.value) /
                            AAPP.value;
                  } else {
                    dcopy(M, A(1, p).asArray(), 1, WORK, 1);
                    dlascl('G', 0, 0, AAPP.value, D[p], M, 1,
                        WORK.asMatrix(LDA), LDA, IERR);
                    AAPQ = ddot(M, WORK, 1, A(1, q).asArray(), 1) *
                        D[q] /
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
                                D[p] *
                                D[q] /
                                AAQQ.value) /
                            AAPP.value;
                  } else {
                    dcopy(M, A(1, q).asArray(), 1, WORK, 1);
                    dlascl('G', 0, 0, AAQQ.value, D[q], M, 1,
                        WORK.asMatrix(LDA), LDA, IERR);
                    AAPQ = ddot(M, WORK, 1, A(1, p).asArray(), 1) *
                        D[p] /
                        AAPP.value;
                  }
                }

                MXAAPQ = max(MXAAPQ, AAPQ.abs());

                // TO rotate or NOT to rotate, THAT is the question ...

                if (AAPQ.abs() > TOL) {
                  NOTROT = 0;
                  // ROTATED++
                  PSKIPPED = 0;
                  ISWROT++;

                  if (ROTOK) {
                    AQOAP = AAQQ.value / AAPP.value;
                    APOAQ = AAPP.value / AAQQ.value;
                    THETA = -HALF * (AQOAP - APOAQ).abs() / AAPQ;
                    if (AAQQ.value > AAPP0) THETA = -THETA;

                    if (THETA.abs() > BIGTHETA) {
                      T.value = HALF / THETA;
                      FASTR[3] = T.value * D[p] / D[q];
                      FASTR[4] = -T.value * D[q] / D[p];
                      drotm(
                          M, A(1, p).asArray(), 1, A(1, q).asArray(), 1, FASTR);
                      if (RSVEC) {
                        drotm(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1,
                            FASTR);
                      }
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ));
                      AAPP.value = AAPP.value *
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
                      AAPP.value = AAPP.value *
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ));

                      APOAQ = D[p] / D[q];
                      AQOAP = D[q] / D[p];
                      if (D[p] >= ONE) {
                        if (D[q] >= ONE) {
                          FASTR[3] = T.value * APOAQ;
                          FASTR[4] = -T.value * AQOAP;
                          D[p] = D[p] * CS;
                          D[q] = D[q] * CS;
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
                          D[p] = D[p] * CS;
                          D[q] = D[q] / CS;
                        }
                      } else {
                        if (D[q] >= ONE) {
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
                          D[p] = D[p] / CS;
                          D[q] = D[q] * CS;
                        } else {
                          if (D[p] >= D[q]) {
                            daxpy(M, -T.value * AQOAP, A(1, q).asArray(), 1,
                                A(1, p).asArray(), 1);
                            daxpy(M, CS * SN * APOAQ, A(1, p).asArray(), 1,
                                A(1, q).asArray(), 1);
                            D[p] = D[p] * CS;
                            D[q] = D[q] / CS;
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
                            D[p] = D[p] / CS;
                            D[q] = D[q] * CS;
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
                      dcopy(M, A(1, p).asArray(), 1, WORK, 1);
                      dlascl('G', 0, 0, AAPP.value, ONE, M, 1,
                          WORK.asMatrix(LDA), LDA, IERR);
                      dlascl(
                          'G', 0, 0, AAQQ.value, ONE, M, 1, A(1, q), LDA, IERR);
                      TEMP1.value = -AAPQ * D[p] / D[q];
                      daxpy(M, TEMP1.value, WORK, 1, A(1, q).asArray(), 1);
                      dlascl(
                          'G', 0, 0, ONE, AAQQ.value, M, 1, A(1, q), LDA, IERR);
                      SVA[q] = AAQQ.value * sqrt(max(ZERO, ONE - AAPQ * AAPQ));
                      MXSINJ = max(MXSINJ, SFMIN);
                    } else {
                      dcopy(M, A(1, q).asArray(), 1, WORK, 1);
                      dlascl('G', 0, 0, AAQQ.value, ONE, M, 1,
                          WORK.asMatrix(LDA), LDA, IERR);
                      dlascl(
                          'G', 0, 0, AAPP.value, ONE, M, 1, A(1, p), LDA, IERR);
                      TEMP1.value = -AAPQ * D[q] / D[p];
                      daxpy(M, TEMP1.value, WORK, 1, A(1, p).asArray(), 1);
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
                      SVA[q] = dnrm2(M, A(1, q).asArray(), 1) * D[q];
                    } else {
                      T.value = ZERO;
                      AAQQ.value = ONE;
                      dlassq(M, A(1, q).asArray(), 1, T, AAQQ);
                      SVA[q] = T.value * sqrt(AAQQ.value) * D[q];
                    }
                  }
                  if (pow((AAPP.value / AAPP0), 2) <= ROOTEPS) {
                    if ((AAPP.value < ROOTBIG) && (AAPP.value > ROOTSFMIN)) {
                      AAPP.value = dnrm2(M, A(1, p).asArray(), 1) * D[p];
                    } else {
                      T.value = ZERO;
                      AAPP.value = ONE;
                      dlassq(M, A(1, p).asArray(), 1, T, AAPP);
                      AAPP.value = T.value * sqrt(AAPP.value) * D[p];
                    }
                    SVA[p] = AAPP.value;
                  }
                  // end of OK rotation
                } else {
                  NOTROT++;
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
      for (p = igl; p <= min(igl + KBL - 1, N); p++) {
        SVA[p] = SVA[p].abs();
      }
    }
    // end of the ibr-loop

    // .. update SVA[N]
    if ((SVA[N] < ROOTBIG) && (SVA[N] > ROOTSFMIN)) {
      SVA[N] = dnrm2(M, A(1, N).asArray(), 1) * D[N];
    } else {
      T.value = ZERO;
      AAPP.value = ONE;
      dlassq(M, A(1, N).asArray(), 1, T, AAPP);
      SVA[N] = T.value * sqrt(AAPP.value) * D[N];
    }

    // Additional steering devices

    if ((i < SWBAND) && ((MXAAPQ <= ROOTTOL) || (ISWROT <= N))) SWBAND = i;

    if ((i > SWBAND + 1) &&
        (MXAAPQ < (N).toDouble() * TOL) &&
        (N.toDouble() * MXAAPQ * MXSINJ < TOL)) {
      isBelowTolerance = true;
      break;
    }

    if (NOTROT >= EMPTSW) {
      isBelowTolerance = true;
      break;
    }
  }
  // end i=1:NSWEEP loop
  if (!isBelowTolerance) {
    // #:) Reaching this point means that the procedure has completed the given
    // number of iterations.
    INFO.value = NSWEEP - 1;
  } else {
    // #:) Reaching this point means that during the i-th sweep all pivots were
    // below the given tolerance, causing early exit.
    INFO.value = 0;
    // #:) INFO.value = 0 confirms successful iterations.
  }
  // Sort the vector D.
  for (p = 1; p <= N - 1; p++) {
    q = idamax(N - p + 1, SVA(p), 1) + p - 1;
    if (p != q) {
      TEMP1.value = SVA[p];
      SVA[p] = SVA[q];
      SVA[q] = TEMP1.value;
      TEMP1.value = D[p];
      D[p] = D[q];
      D[q] = TEMP1.value;
      dswap(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1);
      if (RSVEC) dswap(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1);
    }
  }
}
