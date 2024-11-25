// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dznrm2.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zdotc.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zlassq.dart';
import 'package:dart_lapack/src/zrot.dart';

void zgsvj1(
  final String JOBV,
  final int M,
  final int N,
  final int N1,
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
  final A = A_.having(ld: LDA);
  final V = V_.having(ld: LDV);
  final D = D_.having(length: N);
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
      MXAAPQ = 0,
      MXSINJ = 0,
      ROOTBIG,
      ROOTEPS,
      ROOTSFMIN,
      ROOTTOL,
      SMALL,
      SN,
      TEMP1,
      THETA,
      THSIGN;
  int BLSKIP,
      EMPTSW,
      i,
      ibr,
      igl = 0,
      IJBLSK,
      ISWROT = 0,
      jbc,
      jgl,
      KBL,
      MVL = 0,
      NOTROT,
      NBLC,
      NBLR,
      p,
      PSKIPPED,
      q,
      ROWSKIP,
      SWBAND;
  bool APPLV, ROTOK, RSVEC;
  final IERR = Box(0);
  final AAPP = Box(0.0), AAQQ = Box(0.0), T = Box(0.0);

  // Test the input parameters.
  APPLV = lsame(JOBV, 'A');
  RSVEC = lsame(JOBV, 'V');
  if (!(RSVEC || APPLV || lsame(JOBV, 'N'))) {
    INFO.value = -1;
  } else if (M < 0) {
    INFO.value = -2;
  } else if ((N < 0) || (N > M)) {
    INFO.value = -3;
  } else if (N1 < 0) {
    INFO.value = -4;
  } else if (LDA < M) {
    INFO.value = -6;
  } else if ((RSVEC || APPLV) && (MV < 0)) {
    INFO.value = -9;
  } else if ((RSVEC && (LDV < N)) || (APPLV && (LDV < MV))) {
    INFO.value = -11;
  } else if (TOL <= EPS) {
    INFO.value = -14;
  } else if (NSWEEP < 0) {
    INFO.value = -15;
  } else if (LWORK < M) {
    INFO.value = -17;
  } else {
    INFO.value = 0;
  }

  if (INFO.value != 0) {
    xerbla('ZGSVJ1', -INFO.value);
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
  // LARGE = BIG / sqrt( (M*N) )
  BIGTHETA = ONE / ROOTEPS;
  ROOTTOL = sqrt(TOL);

  // Initialize the right singular vector matrix
  // RSVEC = lsame( JOBV, 'Y' )

  EMPTSW = N1 * (N - N1);
  NOTROT = 0;

  // Row-cyclic pivot strategy with de Rijk's pivoting
  KBL = min(8, N);
  NBLR = N1 ~/ KBL;
  if ((NBLR * KBL) != N1) NBLR++;

  // the tiling is NBLR-by-NBLC [tiles]

  NBLC = (N - N1) ~/ KBL;
  if ((NBLC * KBL) != (N - N1)) NBLC++;
  BLSKIP = pow(KBL, 2).toInt() + 1;
  // [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

  ROWSKIP = min(5, KBL);
  // [TP] ROWSKIP is a tuning parameter.
  SWBAND = 0;
  // [TP] SWBAND is a tuning parameter. It is meaningful and effective
  // if ZGESVJ is used as a computational routine in the preconditioned
  // Jacobi SVD algorithm ZGEJSV.

  // | *   *   * [x] [x] [x]|
  // | *   *   * [x] [x] [x]|    Row-cycling in the NBLR-by-NBLC [x] blocks.
  // | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
  // |[x] [x] [x] *   *   * |
  // |[x] [x] [x] *   *   * |
  // |[x] [x] [x] *   *   * |
  var exhausted = true;
  for (i = 1; i <= NSWEEP; i++) {
    MXAAPQ = ZERO;
    MXSINJ = ZERO;
    ISWROT = 0;

    NOTROT = 0;
    PSKIPPED = 0;

    // Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
    // 1 <= p < q <= N. This is the first step toward a blocked implementation
    // of the rotations. New implementation, based on block transformations,
    // is under development.
    for (ibr = 1; ibr <= NBLR; ibr++) {
      igl = (ibr - 1) * KBL + 1;

      // go to the off diagonal blocks
      igl = (ibr - 1) * KBL + 1;

      jbcLoop:
      for (jbc = 1; jbc <= NBLC; jbc++) {
        jgl = (jbc - 1) * KBL + N1 + 1;

        // doing the block at ( ibr, jbc )
        IJBLSK = 0;
        for (p = igl; p <= min(igl + KBL - 1, N1); p++) {
          AAPP.value = SVA[p];
          if (AAPP.value > ZERO) {
            PSKIPPED = 0;

            for (q = jgl; q <= min(jgl + KBL - 1, N); q++) {
              AAQQ.value = SVA[q];
              if (AAQQ.value > ZERO) {
                AAPP0 = AAPP.value;

                // M x 2 Jacobi SVD

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

                // AAPQ *= CONJG(CWORK(p))*CWORK(q)
                AAPQ1 = -AAPQ.abs();
                MXAAPQ = max(MXAAPQ, -AAPQ1);

                // TO rotate or NOT to rotate, THAT is the question .
                if (AAPQ1.abs() > TOL) {
                  OMPQ = AAPQ / AAPQ.abs().toComplex();
                  NOTROT = 0;
                  // [RTD]      ROTATED++
                  PSKIPPED = 0;
                  ISWROT++;

                  if (ROTOK) {
                    AQOAP = AAQQ.value / AAPP.value;
                    APOAQ = AAPP.value / AAQQ.value;
                    THETA = -HALF * (AQOAP - APOAQ).abs() / AAPQ1;
                    if (AAQQ.value > AAPP0) THETA = -THETA;

                    if (THETA.abs() > BIGTHETA) {
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
                      AAPP.value *=
                          sqrt(max(ZERO, ONE - T.value * AQOAP * AAPQ1));
                      MXSINJ = max(MXSINJ, T.value.abs());
                    } else {
                      // choose correct signum for THETA and rotate
                      THSIGN = -sign(ONE, AAPQ1);
                      if (AAQQ.value > AAPP0) THSIGN = -THSIGN;
                      T.value =
                          ONE / (THETA + THSIGN * sqrt(ONE + THETA * THETA));
                      CS = sqrt(ONE / (ONE + T.value * T.value));
                      SN = T.value * CS;
                      MXSINJ = max(MXSINJ, SN.abs());
                      SVA[q] = AAQQ.value *
                          sqrt(max(ZERO, ONE + T.value * APOAQ * AAPQ1));
                      AAPP.value *=
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
                    // have to use modified Gram-Schmidt like transformation
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
                  // recompute SVA(q), SVA(p)
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
      // 2011 bailed out of the jbc-loop
      for (p = igl; p <= min(igl + KBL - 1, N); p++) {
        SVA[p] = SVA[p].abs();
      }
    }
    // 2000 :: end of the ibr-loop

    // update SVA(N)
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
        (MXAAPQ < sqrt(N) * TOL) &&
        (N * MXAAPQ * MXSINJ < TOL)) {
      exhausted = false;
      break;
    }

    if (NOTROT >= EMPTSW) {
      exhausted = false;
      break;
    }
  }
  // end i=1:NSWEEP loop
  if (exhausted) {
    // Reaching this point means that the procedure has not converged.
    INFO.value = NSWEEP - 1;
  } else {
    // Reaching this point means numerical convergence after the i-th
    // sweep.

    INFO.value = 0;
    // INFO = 0 confirms successful iterations.
  }

  // Sort the vector SVA() of column norms.
  for (p = 1; p <= N - 1; p++) {
    q = idamax(N - p + 1, SVA(p), 1) + p - 1;
    if (p != q) {
      TEMP1 = SVA[p];
      SVA[p] = SVA[q];
      SVA[q] = TEMP1;
      AAPQ = D[p];
      D[p] = D[q];
      D[q] = AAPQ;
      zswap(M, A(1, p).asArray(), 1, A(1, q).asArray(), 1);
      if (RSVEC) zswap(MVL, V(1, p).asArray(), 1, V(1, q).asArray(), 1);
    }
  }
}
