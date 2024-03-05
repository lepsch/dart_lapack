import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlassq.dart';
import 'package:lapack/src/zrot.dart';

void zgsvj0(
  final String JOBV,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> D_,
  final Array<double> SVA_,
  final int MV,
  final Matrix<Complex> V_,
  final int LDV,
  final double EPS,
  final double SFMIN,
  final double TOL,
  final int NSWEEP,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having(length: N);
  final A = A_.having(ld: LDA);
  final V = V_.having(ld: LDV);
  final SVA = SVA_.having(length: N);
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  Complex AAPQ, OMPQ;
  double AAPP0,
      AAPQ1,
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
  final AAPP = Box(0.0), AAQQ = Box(0.0), TEMP1 = Box(0.0), T = Box(0.0);
  final IERR = Box(0);

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

  // #:(
  if (INFO.value != 0) {
    xerbla('ZGSVJ0', -INFO.value);
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

  // .. Row-cyclic Jacobi SVD algorithm with column pivoting ..

  EMPTSW = (N * (N - 1)) ~/ 2;
  NOTROT = 0;

  // .. Row-cyclic pivot strategy with de Rijk's pivoting ..

  SWBAND = 0;
  // [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
  // if ZGESVJ is used as a computational routine in the preconditioned
  // Jacobi SVD algorithm ZGEJSV. For sweeps i=1:SWBAND the procedure
  // works on pivots inside a band-like region around the diagonal.
  // The boundaries are determined dynamically, based on the number of
  // pivots above a threshold.

  KBL = min(8, N);
  // [TP] KBL is a tuning parameter that defines the tile size in the
  // tiling of the p-q loops of pivot pairs. In general, an optimal
  // value of KBL depends on the matrix dimensions and on the
  // parameters of the computer's memory.

  NBL = N ~/ KBL;
  if ((NBL * KBL) != N) NBL = NBL + 1;

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

  // .. Row-cyclic pivot strategy with de Rijk's pivoting ..
  var exhausted = true;
  for (i = 1; i <= NSWEEP; i++) {
    // 1993

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
      // 2000

      igl = (ibr - 1) * KBL + 1;

      for (ir1 = 0; ir1 <= min(LKAHEAD, NBL - ibr); ir1++) {
        // 1002

        igl = igl + ir1 * KBL;

        for (p = igl; p <= min(igl + KBL - 1, N - 1); p++) {
          // 2001

          // .. de Rijk's pivoting

          q = idamax(N - p + 1, SVA(p), 1) + p - 1;
          if (p != q) {
            zswap(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1);
            if (RSVEC) zswap(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1);
            TEMP1.value = SVA[p];
            SVA[p] = SVA[q];
            SVA[q] = TEMP1.value;
            AAPQ = D[p];
            D[p] = D[q];
            D[q] = AAPQ;
          }

          if (ir1 == 0) {
            // Column norms are periodically updated by explicit
            // norm computation.
            // Caveat:
            // Unfortunately, some BLAS implementations compute SNCRM2(M,A(1,p),1)
            // as sqrt(S=zdotc(M,A(1,p),1,A(1,p),1)), which may cause the result to
            // overflow for ||A(:,p)||_2 > sqrt(overflow_threshold), and to
            // underflow for ||A(:,p)||_2 < sqrt(underflow_threshold).
            // Hence, dznrm2 cannot be trusted, not even in the case when
            // the true norm is far from the under(over)flow boundaries.
            // If properly implemented dznrm2 is available, the IF-THEN-ELSE-END IF
            // below should be replaced with "AAPP.value = dznrm2( M, A(1,p), 1 )".

            if ((SVA[p] < ROOTBIG) && (SVA[p] > ROOTSFMIN)) {
              SVA[p] = dznrm2(M, A(1, p).asArray(), 1);
            } else {
              TEMP1.value = ZERO;
              AAPP.value = ONE;
              zlassq(M, A(1, p).asArray(), 1, TEMP1, AAPP);
              SVA[p] = TEMP1.value * sqrt(AAPP.value);
            }
            AAPP.value = SVA[p];
          } else {
            AAPP.value = SVA[p];
          }

          if (AAPP.value > ZERO) {
            PSKIPPED = 0;

            for (q = p + 1; q <= min(igl + KBL - 1, N); q++) {
              // 2002

              AAQQ.value = SVA[q];

              if (AAQQ.value > ZERO) {
                AAPP0 = AAPP.value;
                if (AAQQ.value >= ONE) {
                  ROTOK = (SMALL * AAPP.value) <= AAQQ.value;
                  if (AAPP.value < (BIG / AAQQ.value)) {
                    AAPQ =
                        (zdotc(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) /
                                AAQQ.value.toComplex()) /
                            AAPP.value.toComplex();
                  } else {
                    zcopy(M, A(1, p).asArray(), 1, WORK, 1);
                    zlascl('G', 0, 0, AAPP.value, ONE, M, 1, WORK.asMatrix(LDA),
                        LDA, IERR);
                    AAPQ = zdotc(M, WORK, 1, A(1, q).asArray(), 1) /
                        AAQQ.value.toComplex();
                  }
                } else {
                  ROTOK = AAPP.value <= (AAQQ.value / SMALL);
                  if (AAPP.value > (SMALL / AAQQ.value)) {
                    AAPQ =
                        (zdotc(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) /
                                AAPP.value.toComplex()) /
                            AAQQ.value.toComplex();
                  } else {
                    zcopy(M, A(1, q).asArray(), 1, WORK, 1);
                    zlascl('G', 0, 0, AAQQ.value, ONE, M, 1, WORK.asMatrix(LDA),
                        LDA, IERR);
                    AAPQ = zdotc(M, A(1, p).asArray(), 1, WORK, 1) /
                        AAPP.value.toComplex();
                  }
                }

                // AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q)
                AAPQ1 = -(AAPQ).abs();
                MXAAPQ = max(MXAAPQ, -AAPQ1);

                // TO rotate or NOT to rotate, THAT is the question ...

                if ((AAPQ1).abs() > TOL) {
                  OMPQ = AAPQ / AAPQ.abs().toComplex();

                  // .. rotate
                  // [RTD]      ROTATED = ROTATED + ONE

                  if (ir1 == 0) {
                    NOTROT = 0;
                    PSKIPPED = 0;
                    ISWROT = ISWROT + 1;
                  }

                  if (ROTOK) {
                    AQOAP = AAQQ.value / AAPP.value;
                    APOAQ = AAPP.value / AAQQ.value;
                    THETA = -HALF * (AQOAP - APOAQ).abs() / AAPQ1;

                    if ((THETA).abs() > BIGTHETA) {
                      T.value = HALF / THETA;
                      CS = ONE;
                      zrot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1, CS,
                          OMPQ.conjugate() * T.value.toComplex());
                      if (RSVEC) {
                        zrot(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1,
                            CS, OMPQ.conjugate() * T.value.toComplex());
                      }
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ1));
                      AAPP.value = AAPP.value *
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ1));
                      MXSINJ = max(MXSINJ, (T.value).abs());
                    } else {
                      // .. choose correct signum for THETA and rotate

                      THSIGN = -sign(ONE, AAPQ1).toDouble();
                      T.value =
                          ONE / (THETA + THSIGN * sqrt(ONE + THETA * THETA));
                      CS = sqrt(ONE / (ONE + T.value * T.value));
                      SN = T.value * CS;

                      MXSINJ = max(MXSINJ, (SN).abs());
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ1));
                      AAPP.value = AAPP.value *
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ1));

                      zrot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1, CS,
                          OMPQ.conjugate() * SN.toComplex());
                      if (RSVEC) {
                        zrot(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1,
                            CS, OMPQ.conjugate() * SN.toComplex());
                      }
                    }
                    D[p] = -D[q] * OMPQ;
                  } else {
                    // .. have to use modified Gram-Schmidt like transformation
                    zcopy(M, A(1, p).asArray(), 1, WORK, 1);
                    zlascl('G', 0, 0, AAPP.value, ONE, M, 1, WORK.asMatrix(LDA),
                        LDA, IERR);
                    zlascl(
                        'G', 0, 0, AAQQ.value, ONE, M, 1, A(1, q), LDA, IERR);
                    zaxpy(M, -AAPQ, WORK, 1, A(1, q).asArray(), 1);
                    zlascl(
                        'G', 0, 0, ONE, AAQQ.value, M, 1, A(1, q), LDA, IERR);
                    SVA[q] = AAQQ.value * sqrt(max(ZERO, ONE - AAPQ1 * AAPQ1));
                    MXSINJ = max(MXSINJ, SFMIN);
                  }
                  // END IF ROTOK THEN ... ELSE

                  // In the case of cancellation in updating SVA(q), SVA(p)
                  // recompute SVA(q), SVA(p).

                  if (pow(SVA[q] / AAQQ.value, 2) <= ROOTEPS) {
                    if ((AAQQ.value < ROOTBIG) && (AAQQ.value > ROOTSFMIN)) {
                      SVA[q] = dznrm2(M, A(1, q).asArray(), 1);
                    } else {
                      T.value = ZERO;
                      AAQQ.value = ONE;
                      zlassq(M, A(1, q).asArray(), 1, T, AAQQ);
                      SVA[q] = T.value * sqrt(AAQQ.value);
                    }
                  }
                  if ((AAPP.value / AAPP0) <= ROOTEPS) {
                    if ((AAPP.value < ROOTBIG) && (AAPP.value > ROOTSFMIN)) {
                      AAPP.value = dznrm2(M, A(1, p).asArray(), 1);
                    } else {
                      T.value = ZERO;
                      AAPP.value = ONE;
                      zlassq(M, A(1, p).asArray(), 1, T, AAPP);
                      AAPP.value = T.value * sqrt(AAPP.value);
                    }
                    SVA[p] = AAPP.value;
                  }
                } else {
                  // A(:,p) and A(:,q) already numerically orthogonal
                  if (ir1 == 0) NOTROT = NOTROT + 1;
                  // [RTD]      SKIPPED  = SKIPPED  + 1
                  PSKIPPED = PSKIPPED + 1;
                }
              } else {
                // A(:,q) is zero column
                if (ir1 == 0) NOTROT = NOTROT + 1;
                PSKIPPED = PSKIPPED + 1;
              }

              if ((i <= SWBAND) && (PSKIPPED > ROWSKIP)) {
                if (ir1 == 0) AAPP.value = -AAPP.value;
                NOTROT = 0;
                break;
              }
            } // 2002
            // END q-LOOP

            //  } // 2103
            // bailed out of q-loop

            SVA[p] = AAPP.value;
          } else {
            SVA[p] = AAPP.value;
            if ((ir1 == 0) && (AAPP.value == ZERO)) {
              NOTROT = NOTROT + min(igl + KBL - 1, N).toInt() - p;
            }
          }
        } // 2001
        // end of the p-loop
        // end of doing the block ( ibr, ibr )
      } // 1002
      // end of ir1-loop

      // ... go to the off diagonal blocks

      igl = (ibr - 1) * KBL + 1;
      jbcLoop:
      for (jbc = ibr + 1; jbc <= NBL; jbc++) {
        // 2010

        jgl = (jbc - 1) * KBL + 1;

        // doing the block at ( ibr, jbc )

        IJBLSK = 0;
        for (p = igl; p <= min(igl + KBL - 1, N); p++) {
          // 2100

          AAPP.value = SVA[p];
          if (AAPP.value > ZERO) {
            PSKIPPED = 0;

            for (q = jgl; q <= min(jgl + KBL - 1, N); q++) {
              // 2200

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
                        (zdotc(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) /
                                AAQQ.value.toComplex()) /
                            AAPP.value.toComplex();
                  } else {
                    zcopy(M, A(1, p).asArray(), 1, WORK, 1);
                    zlascl('G', 0, 0, AAPP.value, ONE, M, 1, WORK.asMatrix(LDA),
                        LDA, IERR);
                    AAPQ = zdotc(M, WORK, 1, A(1, q).asArray(), 1) /
                        AAQQ.value.toComplex();
                  }
                } else {
                  if (AAPP.value >= AAQQ.value) {
                    ROTOK = AAPP.value <= (AAQQ.value / SMALL);
                  } else {
                    ROTOK = AAQQ.value <= (AAPP.value / SMALL);
                  }
                  if (AAPP.value > (SMALL / AAQQ.value)) {
                    AAPQ =
                        (zdotc(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1) /
                                max(AAQQ.value, AAPP.value).toComplex()) /
                            min(AAQQ.value, AAPP.value).toComplex();
                  } else {
                    zcopy(M, A(1, q).asArray(), 1, WORK, 1);
                    zlascl('G', 0, 0, AAQQ.value, ONE, M, 1, WORK.asMatrix(LDA),
                        LDA, IERR);
                    AAPQ = zdotc(M, A(1, p).asArray(), 1, WORK, 1) /
                        AAPP.value.toComplex();
                  }
                }

                // AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                AAPQ1 = -(AAPQ).abs();
                MXAAPQ = max(MXAAPQ, -AAPQ1);

                // TO rotate or NOT to rotate, THAT is the question ...

                if ((AAPQ1).abs() > TOL) {
                  OMPQ = AAPQ / AAPQ.abs().toComplex();
                  NOTROT = 0;
                  // [RTD]      ROTATED  = ROTATED + 1
                  PSKIPPED = 0;
                  ISWROT = ISWROT + 1;

                  if (ROTOK) {
                    AQOAP = AAQQ.value / AAPP.value;
                    APOAQ = AAPP.value / AAQQ.value;
                    THETA = -HALF * (AQOAP - APOAQ).abs() / AAPQ1;
                    if (AAQQ.value > AAPP0) THETA = -THETA;

                    if ((THETA).abs() > BIGTHETA) {
                      T.value = HALF / THETA;
                      CS = ONE;
                      zrot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1, CS,
                          OMPQ.conjugate() * T.value.toComplex());
                      if (RSVEC) {
                        zrot(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1,
                            CS, OMPQ.conjugate() * T.value.toComplex());
                      }
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ1));
                      AAPP.value = AAPP.value *
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ1));
                      MXSINJ = max(MXSINJ, (T.value).abs());
                    } else {
                      // .. choose correct signum for THETA and rotate

                      THSIGN = -sign(ONE, AAPQ1).toDouble();
                      if (AAQQ.value > AAPP0) THSIGN = -THSIGN;
                      T.value =
                          ONE / (THETA + THSIGN * sqrt(ONE + THETA * THETA));
                      CS = sqrt(ONE / (ONE + T.value * T.value));
                      SN = T.value * CS;
                      MXSINJ = max(MXSINJ, (SN).abs());
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ1));
                      AAPP.value = AAPP.value *
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ1));

                      zrot(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1, CS,
                          OMPQ.conjugate() * SN.toComplex());
                      if (RSVEC) {
                        zrot(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1,
                            CS, OMPQ.conjugate() * SN.toComplex());
                      }
                    }
                    D[p] = -D[q] * OMPQ;
                  } else {
                    // .. have to use modified Gram-Schmidt like transformation
                    if (AAPP.value > AAQQ.value) {
                      zcopy(M, A(1, p).asArray(), 1, WORK, 1);
                      zlascl('G', 0, 0, AAPP.value, ONE, M, 1,
                          WORK.asMatrix(LDA), LDA, IERR);
                      zlascl(
                          'G', 0, 0, AAQQ.value, ONE, M, 1, A(1, q), LDA, IERR);
                      zaxpy(M, -AAPQ, WORK, 1, A(1, q).asArray(), 1);
                      zlascl(
                          'G', 0, 0, ONE, AAQQ.value, M, 1, A(1, q), LDA, IERR);
                      SVA[q] =
                          AAQQ.value * sqrt(max(ZERO, ONE - AAPQ1 * AAPQ1));
                      MXSINJ = max(MXSINJ, SFMIN);
                    } else {
                      zcopy(M, A(1, q).asArray(), 1, WORK, 1);
                      zlascl('G', 0, 0, AAQQ.value, ONE, M, 1,
                          WORK.asMatrix(LDA), LDA, IERR);
                      zlascl(
                          'G', 0, 0, AAPP.value, ONE, M, 1, A(1, p), LDA, IERR);
                      zaxpy(
                          M, -AAPQ.conjugate(), WORK, 1, A(1, p).asArray(), 1);
                      zlascl(
                          'G', 0, 0, ONE, AAPP.value, M, 1, A(1, p), LDA, IERR);
                      SVA[p] =
                          AAPP.value * sqrt(max(ZERO, ONE - AAPQ1 * AAPQ1));
                      MXSINJ = max(MXSINJ, SFMIN);
                    }
                  }
                  // END IF ROTOK THEN ... ELSE

                  // In the case of cancellation in updating SVA(q), SVA(p)
                  // .. recompute SVA(q), SVA(p)
                  if (pow(SVA[q] / AAQQ.value, 2) <= ROOTEPS) {
                    if ((AAQQ.value < ROOTBIG) && (AAQQ.value > ROOTSFMIN)) {
                      SVA[q] = dznrm2(M, A(1, q).asArray(), 1);
                    } else {
                      T.value = ZERO;
                      AAQQ.value = ONE;
                      zlassq(M, A(1, q).asArray(), 1, T, AAQQ);
                      SVA[q] = T.value * sqrt(AAQQ.value);
                    }
                  }
                  if (pow(AAPP.value / AAPP0, 2) <= ROOTEPS) {
                    if ((AAPP.value < ROOTBIG) && (AAPP.value > ROOTSFMIN)) {
                      AAPP.value = dznrm2(M, A(1, p).asArray(), 1);
                    } else {
                      T.value = ZERO;
                      AAPP.value = ONE;
                      zlassq(M, A(1, p).asArray(), 1, T, AAPP);
                      AAPP.value = T.value * sqrt(AAPP.value);
                    }
                    SVA[p] = AAPP.value;
                  }
                  // end of OK rotation
                } else {
                  NOTROT = NOTROT + 1;
                  // [RTD]      SKIPPED  = SKIPPED  + 1
                  PSKIPPED = PSKIPPED + 1;
                  IJBLSK = IJBLSK + 1;
                }
              } else {
                NOTROT = NOTROT + 1;
                PSKIPPED = PSKIPPED + 1;
                IJBLSK = IJBLSK + 1;
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
            } // 2200
            // end of the q-loop
            //  } // 2203

            SVA[p] = AAPP.value;
          } else {
            if (AAPP.value == ZERO) {
              NOTROT = NOTROT + min(jgl + KBL - 1, N).toInt() - jgl + 1;
            }
            if (AAPP.value < ZERO) NOTROT = 0;
          }
        } // 2100
        // end of the p-loop
      } // 2010
      // end of the jbc-loop
      //} // 2011
      // 2011 bailed out of the jbc-loop
      for (p = igl; p <= min(igl + KBL - 1, N); p++) {
        // 2012
        SVA[p] = (SVA[p]).abs();
      } // 2012
    } // 2000
    // 2000 :: end of the ibr-loop

    // .. update SVA(N)
    if ((SVA[N] < ROOTBIG) && (SVA[N] > ROOTSFMIN)) {
      SVA[N] = dznrm2(M, A(1, N).asArray(), 1);
    } else {
      T.value = ZERO;
      AAPP.value = ONE;
      zlassq(M, A(1, N).asArray(), 1, T, AAPP);
      SVA[N] = T.value * sqrt(AAPP.value);
    }

    // Additional steering devices

    if ((i < SWBAND) && ((MXAAPQ <= ROOTTOL) || (ISWROT <= N))) SWBAND = i;

    if ((i > SWBAND + 1) &&
        (MXAAPQ < sqrt((N).toDouble()) * TOL) &&
        (N.toDouble() * MXAAPQ * MXSINJ < TOL)) {
      exhausted = false;
      break;
    }

    if (NOTROT >= EMPTSW) {
      exhausted = false;
      break;
    }
  } // 1993
  // end i=1:NSWEEP loop
  if (exhausted) {
    // #:( Reaching this point means that the procedure has not converged.
    INFO.value = NSWEEP - 1;
  } else {
    // #:) Reaching this point means numerical convergence after the i-th
    // sweep.

    INFO.value = 0;
    // #:) INFO.value = 0 confirms successful iterations.
  } // 1995

  // Sort the vector SVA() of column norms.
  for (p = 1; p <= N - 1; p++) {
    // 5991
    q = idamax(N - p + 1, SVA(p), 1) + p - 1;
    if (p != q) {
      TEMP1.value = SVA[p];
      SVA[p] = SVA[q];
      SVA[q] = TEMP1.value;
      AAPQ = D[p];
      D[p] = D[q];
      D[q] = AAPQ;
      zswap(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1);
      if (RSVEC) zswap(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1);
    }
  } // 5991
}
