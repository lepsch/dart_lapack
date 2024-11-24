// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:collection/collection.dart';
import 'package:lapack/lapack.dart';

import 'alaerh.dart';
import 'alahd.dart';
import 'alasvm.dart';
import 'common.dart';
import 'xlaenv.dart';
import 'zerrls.dart';
import 'zqrt12.dart';
import 'zqrt13.dart';
import 'zqrt14.dart';
import 'zqrt15.dart';
import 'zqrt16.dart';
import 'zqrt17.dart';

void zdrvls(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final Array<int> NXVAL_,
  final double THRESH,
  final bool TSTERR,
  final Array<Complex> A_,
  final Array<Complex> COPYA_,
  final Array<Complex> B_,
  final Array<Complex> COPYB_,
  final Array<Complex> C_,
  final Array<double> S_,
  final Array<double> COPYS_,
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final NSVAL = NSVAL_.having();
  final NBVAL = NBVAL_.having();
  final NXVAL = NXVAL_.having();
  final A = A_.having();
  final COPYA = COPYA_.having();
  final B = B_.having();
  final COPYB = COPYB_.having();
  final C = C_.having();
  final S = S_.having();
  final COPYS = COPYS_.having();
  const NTESTS = 18, SMLSIZ = 25;
  const ONE = 1.0, ZERO = 0.0;
  final IWQ = Array<int>(1);
  final RESULT = Array<double>(NTESTS), RWQ = Array<double>(1);
  final WQ = Array<Complex>(1);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}LS';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);
  final EPS = dlamch('Epsilon');

  // Threshold for rank estimation

  final RCOND = sqrt(EPS) - (sqrt(EPS) - EPS) / 2;

  // Test the error exits

  xlaenv(9, SMLSIZ);
  if (TSTERR) zerrls(PATH, NOUT);

  // Print the header if NM = 0 or NN = 0 and THRESH = 0.
  if ((NM == 0 || NN == 0) && THRESH == ZERO) alahd(NOUT, PATH);
  infoc.INFOT = 0;

  // Compute maximal workspace needed for all routines
  var NMAX = 0;
  var MMAX = 0;
  var NSMAX = 0;
  for (var I = 1; I <= NM; I++) {
    if (MVAL[I] > MMAX) {
      MMAX = MVAL[I];
    }
  }
  for (var I = 1; I <= NN; I++) {
    if (NVAL[I] > NMAX) {
      NMAX = NVAL[I];
    }
  }
  for (var I = 1; I <= NNS; I++) {
    if (NSVAL[I] > NSMAX) {
      NSMAX = NSVAL[I];
    }
  }
  final M = MMAX;
  final N = NMAX;
  final NRHS = NSMAX;
  final MNMIN = max(min(M, N), 1);

  // Compute workspace needed for routines
  // ZQRT14, ZQRT17 (two side cases), ZQRT15 and ZQRT12 
  var LWORK = <int>[
    1,
    (M + N) * NRHS,
    (N + NRHS) * (M + 2),
    (M + NRHS) * (N + 2),
    max(M + MNMIN, max(NRHS * MNMIN, 2 * N + M)),
    max(M * N + 4 * MNMIN + max(M, N), M * N + 2 * MNMIN + 4 * N)
  ].max;
  var LRWORK = 1;
  var LIWORK = 1;

  // Iterate through all test cases and compute necessary workspace
  // sizes for ?GELS, ?GELST, ?GETSLS, ?GELSY, ?GELSS and ?GELSD
  // routines.

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];
    final LDA = max(1, M);
    for (var IN = 1; IN <= NN; IN++) {
      final N = NVAL[IN];
      final MNMIN = max(min(M, N), 1);
      final LDB = max(1, max(M, N));
      for (var INS = 1; INS <= NNS; INS++) {
        final NRHS = NSVAL[INS];
        var LWORK_ZGELS = 0;
        var LWORK_ZGELST = 0;
        var LWORK_ZGETSLS = 0;
        for (var IRANK = 1; IRANK <= 2; IRANK++) {
          for (var ISCALE = 1; ISCALE <= 3; ISCALE++) {
            final ITYPE = (IRANK - 1) * 3 + ISCALE;
            if (DOTYPE[ITYPE]) {
              if (IRANK == 1) {
                for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
                  final TRANS = ITRAN == 1 ? 'N' : 'C';

                  // Compute workspace needed for ZGELS
                  zgels(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB,
                      WQ, -1, INFO);
                  LWORK_ZGELS = WQ[1].toInt();
                  // Compute workspace needed for ZGELST
                  zgelst(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(),
                      LDB, WQ, -1, INFO);
                  LWORK_ZGELST = WQ[1].toInt();
                  // Compute workspace needed for ZGETSLS
                  zgetsls(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(),
                      LDB, WQ, -1, INFO);
                  LWORK_ZGETSLS = WQ[1].toInt();
                }
              }
              final CRANK = Box(0);
              // Compute workspace needed for ZGELSY
              zgelsy(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, IWQ,
                  RCOND, CRANK, WQ, -1, RWQ, INFO);
              final LWORK_ZGELSY = WQ[1].toInt();
              final LRWORK_ZGELSY = 2 * N;
              // Compute workspace needed for ZGELSS
              zgelss(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, S, RCOND,
                  CRANK, WQ, -1, RWQ, INFO);
              final LWORK_ZGELSS = WQ[1].toInt();
              final LRWORK_ZGELSS = 5 * MNMIN;
              // Compute workspace needed for ZGELSD
              zgelsd(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, S, RCOND,
                  CRANK, WQ, -1, RWQ, IWQ, INFO);
              final LWORK_ZGELSD = WQ[1].toInt();
              final LRWORK_ZGELSD = RWQ[1].toInt();
              // Compute LIWORK workspace needed for ZGELSY and ZGELSD
              LIWORK = [LIWORK, N, IWQ[1]].max;
              // Compute LRWORK workspace needed for ZGELSY, ZGELSS and ZGELSD
              LRWORK =
                  [LRWORK, LRWORK_ZGELSY, LRWORK_ZGELSS, LRWORK_ZGELSD].max;
              // Compute LWORK workspace needed for all functions
              LWORK = [
                LWORK,
                LWORK_ZGELS,
                LWORK_ZGELST,
                LWORK_ZGETSLS,
                LWORK_ZGELSY,
                LWORK_ZGELSS,
                LWORK_ZGELSD
              ].max;
            }
          }
        }
      }
    }
  }

  final LWLSY = LWORK;

  final WORK = Array<Complex>(LWORK);
  final WORK2 = Array<double>(2 * LWORK);
  final IWORK = Array<int>(LIWORK);
  final RWORK = Array<double>(LRWORK);

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];
    final LDA = max(1, M);

    for (var IN = 1; IN <= NN; IN++) {
      final N = NVAL[IN];
      final MNMIN = max(min(M, N), 1);
      final LDB = max(1, max(M, N));
      // final MB = (MNMIN+1);

      for (var INS = 1; INS <= NNS; INS++) {
        final NRHS = NSVAL[INS];

        for (var IRANK = 1; IRANK <= 2; IRANK++) {
          for (var ISCALE = 1; ISCALE <= 3; ISCALE++) {
            final ITYPE = (IRANK - 1) * 3 + ISCALE;
            if (!DOTYPE[ITYPE]) continue;
            // =====================================================
            //       Begin test ZGELS
            // =====================================================
            if (IRANK == 1) {
              // Generate a matrix of scaling type ISCALE

              final NORMA = Box(ZERO);
              zqrt13(ISCALE, M, N, COPYA.asMatrix(), LDA, NORMA, ISEED);

              // Loop for testing different block sizes.

              for (var INB = 1; INB <= NNB; INB++) {
                final NB = NBVAL[INB];
                xlaenv(1, NB);
                xlaenv(3, NXVAL[INB]);

                // Loop for testing non-transposed and transposed.

                for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
                  final (TRANS, NROWS, NCOLS) =
                      ITRAN == 1 ? ('N', M, N) : ('C', N, M);
                  final LDWORK = max(1, NCOLS);

                  // Set up a consistent rhs

                  if (NCOLS > 0) {
                    zlarnv(2, ISEED, NCOLS * NRHS, WORK);
                    zdscal(NCOLS * NRHS, ONE / NCOLS, WORK, 1);
                  }
                  zgemm(
                      TRANS,
                      'No transpose',
                      NROWS,
                      NRHS,
                      NCOLS,
                      Complex.one,
                      COPYA.asMatrix(),
                      LDA,
                      WORK.asMatrix(),
                      LDWORK,
                      Complex.zero,
                      B.asMatrix(),
                      LDB);
                  zlacpy('Full', NROWS, NRHS, B.asMatrix(), LDB,
                      COPYB.asMatrix(), LDB);

                  // Solve LS or overdetermined system

                  if (M > 0 && N > 0) {
                    zlacpy(
                        'Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
                    zlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                        B.asMatrix(), LDB);
                  }
                  srnamc.SRNAMT = 'ZGELS';
                  zgels(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB,
                      WORK, LWORK, INFO);

                  if (INFO.value != 0) {
                    alaerh(PATH, 'ZGELS ', INFO.value, 0, TRANS, M, N, NRHS, -1,
                        NB, ITYPE, NFAIL, NERRS, NOUT);
                  }

                  // Test 1: Check correctness of results
                  // for ZGELS, compute the residual:
                  // RESID = norm(B - A*X) /
                  // / ( max(m,n) * norm(A) * norm(X) * EPS )

                  if (NROWS > 0 && NRHS > 0) {
                    zlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                        C.asMatrix(), LDB);
                  }
                  zqrt16(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA, B.asMatrix(),
                      LDB, C.asMatrix(), LDB, RWORK, RESULT(1));

                  // Test 2: Check correctness of results
                  // for ZGELS.

                  if ((ITRAN == 1 && M >= N) || (ITRAN == 2 && M < N)) {
                    // Solving LS system

                    RESULT[2] = zqrt17(
                        TRANS,
                        1,
                        M,
                        N,
                        NRHS,
                        COPYA.asMatrix(),
                        LDA,
                        B.asMatrix(),
                        LDB,
                        COPYB.asMatrix(),
                        LDB,
                        C.asMatrix(),
                        WORK,
                        LWORK);
                  } else {
                    // Solving overdetermined system

                    RESULT[2] = zqrt14(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA,
                        B.asMatrix(), LDB, WORK, LWORK);
                  }

                  // Print information about the tests that
                  // did not pass the threshold.

                  for (var K = 1; K <= 2; K++) {
                    if (RESULT[K] >= THRESH) {
                      if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                      NOUT.print9999(
                          TRANS, M, N, NRHS, NB, ITYPE, K, RESULT[K]);
                      NFAIL++;
                    }
                  }
                  NRUN += 2;
                }
              }
            }
            // =====================================================
            //       End test ZGELS
            // =====================================================
            // =====================================================
            //       Begin test ZGELST
            // =====================================================
            if (IRANK == 1) {
              // Generate a matrix of scaling type ISCALE

              final NORMA = Box(ZERO);
              zqrt13(ISCALE, M, N, COPYA.asMatrix(), LDA, NORMA, ISEED);

              // Loop for testing different block sizes.

              for (var INB = 1; INB <= NNB; INB++) {
                final NB = NBVAL[INB];
                xlaenv(1, NB);
                xlaenv(3, NXVAL[INB]);

                // Loop for testing non-transposed and transposed.

                for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
                  final (TRANS, NROWS, NCOLS) =
                      ITRAN == 1 ? ('N', M, N) : ('C', N, M);
                  final LDWORK = max(1, NCOLS);

                  // Set up a consistent rhs

                  if (NCOLS > 0) {
                    zlarnv(2, ISEED, NCOLS * NRHS, WORK);
                    zdscal(NCOLS * NRHS, ONE / NCOLS, WORK, 1);
                  }
                  zgemm(
                      TRANS,
                      'No transpose',
                      NROWS,
                      NRHS,
                      NCOLS,
                      Complex.one,
                      COPYA.asMatrix(),
                      LDA,
                      WORK.asMatrix(),
                      LDWORK,
                      Complex.zero,
                      B.asMatrix(),
                      LDB);
                  zlacpy('Full', NROWS, NRHS, B.asMatrix(), LDB,
                      COPYB.asMatrix(), LDB);

                  // Solve LS or overdetermined system

                  if (M > 0 && N > 0) {
                    zlacpy(
                        'Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
                    zlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                        B.asMatrix(), LDB);
                  }
                  srnamc.SRNAMT = 'ZGELST';
                  zgelst(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(),
                      LDB, WORK, LWORK, INFO);

                  if (INFO.value != 0) {
                    alaerh(PATH, 'ZGELST', INFO.value, 0, TRANS, M, N, NRHS, -1,
                        NB, ITYPE, NFAIL, NERRS, NOUT);
                  }

                  // Test 3: Check correctness of results
                  // for ZGELST, compute the residual:
                  // RESID = norm(B - A*X) /
                  // / ( max(m,n) * norm(A) * norm(X) * EPS )

                  if (NROWS > 0 && NRHS > 0) {
                    zlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                        C.asMatrix(), LDB);
                  }
                  zqrt16(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA, B.asMatrix(),
                      LDB, C.asMatrix(), LDB, RWORK, RESULT(3));

                  // Test 4: Check correctness of results
                  // for ZGELST.

                  if ((ITRAN == 1 && M >= N) || (ITRAN == 2 && M < N)) {
                    // Solving LS system

                    RESULT[4] = zqrt17(
                        TRANS,
                        1,
                        M,
                        N,
                        NRHS,
                        COPYA.asMatrix(),
                        LDA,
                        B.asMatrix(),
                        LDB,
                        COPYB.asMatrix(),
                        LDB,
                        C.asMatrix(),
                        WORK,
                        LWORK);
                  } else {
                    // Solving overdetermined system

                    RESULT[4] = zqrt14(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA,
                        B.asMatrix(), LDB, WORK, LWORK);
                  }

                  // Print information about the tests that
                  // did not pass the threshold.

                  for (var K = 3; K <= 4; K++) {
                    if (RESULT[K] >= THRESH) {
                      if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                      NOUT.print9999(
                          TRANS, M, N, NRHS, NB, ITYPE, K, RESULT[K]);
                      NFAIL++;
                    }
                  }
                  NRUN += 2;
                }
              }
            }
            // =====================================================
            //       End test ZGELST
            // =====================================================
            // =====================================================
            //       Begin test ZGELSTSLS
            // =====================================================
            if (IRANK == 1) {
              // Generate a matrix of scaling type ISCALE

              final NORMA = Box(ZERO);
              zqrt13(ISCALE, M, N, COPYA.asMatrix(), LDA, NORMA, ISEED);

              // Loop for testing different block sizes MB.

              for (var INB = 1; INB <= NNB; INB++) {
                final MB = NBVAL[INB];
                xlaenv(1, MB);

                // Loop for testing different block sizes NB.

                for (var IMB = 1; IMB <= NNB; IMB++) {
                  final NB = NBVAL[IMB];
                  xlaenv(2, NB);

                  // Loop for testing non-transposed
                  // and transposed.

                  for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
                    final (TRANS, NROWS, NCOLS) =
                        ITRAN == 1 ? ('N', M, N) : ('C', N, M);
                    final LDWORK = max(1, NCOLS);

                    // Set up a consistent rhs

                    if (NCOLS > 0) {
                      zlarnv(2, ISEED, NCOLS * NRHS, WORK);
                      zscal(NCOLS * NRHS, Complex.one / NCOLS.toComplex(), WORK,
                          1);
                    }
                    zgemm(
                        TRANS,
                        'No transpose',
                        NROWS,
                        NRHS,
                        NCOLS,
                        Complex.one,
                        COPYA.asMatrix(),
                        LDA,
                        WORK.asMatrix(),
                        LDWORK,
                        Complex.zero,
                        B.asMatrix(),
                        LDB);
                    zlacpy('Full', NROWS, NRHS, B.asMatrix(), LDB,
                        COPYB.asMatrix(), LDB);

                    // Solve LS or overdetermined system

                    if (M > 0 && N > 0) {
                      zlacpy('Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(),
                          LDA);
                      zlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                          B.asMatrix(), LDB);
                    }
                    srnamc.SRNAMT = 'ZGETSLS';
                    zgetsls(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(),
                        LDB, WORK, LWORK, INFO);
                    if (INFO.value != 0) {
                      alaerh(PATH, 'ZGETSLS ', INFO.value, 0, TRANS, M, N, NRHS,
                          -1, NB, ITYPE, NFAIL, NERRS, NOUT);
                    }

                    // Test 5: Check correctness of results
                    // for ZGETSLS, compute the residual:
                    // RESID = norm(B - A*X) /
                    // / ( max(m,n) * norm(A) * norm(X) * EPS )

                    if (NROWS > 0 && NRHS > 0) {
                      zlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                          C.asMatrix(), LDB);
                    }
                    zqrt16(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA,
                        B.asMatrix(), LDB, C.asMatrix(), LDB, WORK2, RESULT(5));

                    // Test 6: Check correctness of results
                    // for ZGETSLS.

                    if ((ITRAN == 1 && M >= N) || (ITRAN == 2 && M < N)) {
                      // Solving LS system, compute:
                      // r = norm((B- A*X)**T * A) /
                      // / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                      RESULT[6] = zqrt17(
                          TRANS,
                          1,
                          M,
                          N,
                          NRHS,
                          COPYA.asMatrix(),
                          LDA,
                          B.asMatrix(),
                          LDB,
                          COPYB.asMatrix(),
                          LDB,
                          C.asMatrix(),
                          WORK,
                          LWORK);
                    } else {
                      // Solving overdetermined system

                      RESULT[6] = zqrt14(TRANS, M, N, NRHS, COPYA.asMatrix(),
                          LDA, B.asMatrix(), LDB, WORK, LWORK);
                    }

                    // Print information about the tests that
                    // did not pass the threshold.

                    for (var K = 5; K <= 6; K++) {
                      if (RESULT[K] >= THRESH) {
                        if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                        NOUT.println(
                            ' TRANS=\'${TRANS.a1} M=${M.i5}, N=${N.i5}, NRHS=${NRHS.i4}, MB=${MB.i4}, NB=${NB.i4}, type${ITYPE.i2}, test(${K.i2})=${RESULT[K].g12_5}');
                        NFAIL++;
                      }
                    }
                    NRUN += 2;
                  }
                }
              }
            }
            // =====================================================
            //       End test ZGELSTSLS
            // =====================================================

            // Generate a matrix of scaling type ISCALE and rank
            // type IRANK.

            final RANK = Box(0);
            final NORMA = Box(ZERO), NORMB = Box(ZERO);
            zqrt15(
                ISCALE,
                IRANK,
                M,
                N,
                NRHS,
                COPYA.asMatrix(),
                LDA,
                COPYB.asMatrix(),
                LDB,
                COPYS,
                RANK,
                NORMA,
                NORMB,
                ISEED,
                WORK,
                LWORK);

            // workspace used: max(M+min(M,N),NRHS*min(M,N),2*N+M)

            final LDWORK = max(1, M);

            // Loop for testing different block sizes.

            for (var INB = 1; INB <= NNB; INB++) {
              final NB = NBVAL[INB];
              xlaenv(1, NB);
              xlaenv(3, NXVAL[INB]);

              // Test ZGELSY

              // ZGELSY:  Compute the minimum-norm solution
              // X to min( norm( A * X - B ) )
              // using the rank-revealing orthogonal
              // factorization.

              zlacpy('Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
              zlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, B.asMatrix(), LDB);

              // Initialize vector IWORK.

              for (var J = 1; J <= N; J++) {
                IWORK[J] = 0;
              }

              final CRANK = Box(0);
              srnamc.SRNAMT = 'ZGELSY';
              zgelsy(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, IWORK,
                  RCOND, CRANK, WORK, LWLSY, RWORK, INFO);
              if (INFO.value != 0) {
                alaerh(PATH, 'ZGELSY', INFO.value, 0, ' ', M, N, NRHS, -1, NB,
                    ITYPE, NFAIL, NERRS, NOUT);
              }

              // workspace used: 2*MNMIN+NB*NB+NB*max(N,NRHS)

              // Test 7:  Compute relative error in svd
              //          workspace: M*N + 4*min(M,N) + max(M,N)

              RESULT[7] = zqrt12(CRANK.value, CRANK.value, A.asMatrix(), LDA,
                  COPYS, WORK, LWORK, RWORK);

              // Test 8:  Compute error in solution
              //          workspace:  M*NRHS + M

              zlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, WORK.asMatrix(),
                  LDWORK);
              zqrt16('No transpose', M, N, NRHS, COPYA.asMatrix(), LDA,
                  B.asMatrix(), LDB, WORK.asMatrix(), LDWORK, RWORK, RESULT(8));

              // Test 9:  Check norm of r'*A
              //          workspace: NRHS*(M+N)

              RESULT[9] = ZERO;
              if (M > CRANK.value) {
                RESULT[9] = zqrt17(
                    'No transpose',
                    1,
                    M,
                    N,
                    NRHS,
                    COPYA.asMatrix(),
                    LDA,
                    B.asMatrix(),
                    LDB,
                    COPYB.asMatrix(),
                    LDB,
                    C.asMatrix(),
                    WORK,
                    LWORK);
              }

              // Test 10:  Check if x is in the rowspace of A
              //          workspace: (M+NRHS)*(N+2)

              RESULT[10] = ZERO;

              if (N > CRANK.value) {
                RESULT[10] = zqrt14('No transpose', M, N, NRHS,
                    COPYA.asMatrix(), LDA, B.asMatrix(), LDB, WORK, LWORK);
              }

              // Test ZGELSS

              // ZGELSS:  Compute the minimum-norm solution
              // X to min( norm( A * X - B ) )
              // using the SVD.

              zlacpy('Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
              zlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, B.asMatrix(), LDB);
              srnamc.SRNAMT = 'ZGELSS';
              zgelss(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, S, RCOND,
                  CRANK, WORK, LWORK, RWORK, INFO);

              if (INFO.value != 0) {
                alaerh(PATH, 'ZGELSS', INFO.value, 0, ' ', M, N, NRHS, -1, NB,
                    ITYPE, NFAIL, NERRS, NOUT);
              }

              // workspace used: 3*min(m,n) +
              //                 max(2*min(m,n),nrhs,max(m,n))

              // Test 11:  Compute relative error in svd

              if (RANK.value > 0) {
                daxpy(MNMIN, -ONE, COPYS, 1, S, 1);
                RESULT[11] =
                    dasum(MNMIN, S, 1) / dasum(MNMIN, COPYS, 1) / (EPS * MNMIN);
              } else {
                RESULT[11] = ZERO;
              }

              // Test 12:  Compute error in solution

              zlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, WORK.asMatrix(),
                  LDWORK);
              zqrt16(
                  'No transpose',
                  M,
                  N,
                  NRHS,
                  COPYA.asMatrix(),
                  LDA,
                  B.asMatrix(),
                  LDB,
                  WORK.asMatrix(),
                  LDWORK,
                  RWORK,
                  RESULT(12));

              // Test 13:  Check norm of r'*A

              RESULT[13] = ZERO;
              if (M > CRANK.value) {
                RESULT[13] = zqrt17(
                    'No transpose',
                    1,
                    M,
                    N,
                    NRHS,
                    COPYA.asMatrix(),
                    LDA,
                    B.asMatrix(),
                    LDB,
                    COPYB.asMatrix(),
                    LDB,
                    C.asMatrix(),
                    WORK,
                    LWORK);
              }

              // Test 14:  Check if x is in the rowspace of A

              RESULT[14] = ZERO;
              if (N > CRANK.value) {
                RESULT[14] = zqrt14('No transpose', M, N, NRHS,
                    COPYA.asMatrix(), LDA, B.asMatrix(), LDB, WORK, LWORK);
              }

              // Test ZGELSD

              // ZGELSD:  Compute the minimum-norm solution X
              // to min( norm( A * X - B ) ) using a
              // divide and conquer SVD.

              xlaenv(9, 25);

              zlacpy('Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
              zlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, B.asMatrix(), LDB);

              srnamc.SRNAMT = 'ZGELSD';
              zgelsd(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, S, RCOND,
                  CRANK, WORK, LWORK, RWORK, IWORK, INFO);
              if (INFO.value != 0) {
                alaerh(PATH, 'ZGELSD', INFO.value, 0, ' ', M, N, NRHS, -1, NB,
                    ITYPE, NFAIL, NERRS, NOUT);
              }

              // Test 15:  Compute relative error in svd

              if (RANK.value > 0) {
                daxpy(MNMIN, -ONE, COPYS, 1, S, 1);
                RESULT[15] =
                    dasum(MNMIN, S, 1) / dasum(MNMIN, COPYS, 1) / (EPS * MNMIN);
              } else {
                RESULT[15] = ZERO;
              }

              // Test 16:  Compute error in solution

              zlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, WORK.asMatrix(),
                  LDWORK);
              zqrt16(
                  'No transpose',
                  M,
                  N,
                  NRHS,
                  COPYA.asMatrix(),
                  LDA,
                  B.asMatrix(),
                  LDB,
                  WORK.asMatrix(),
                  LDWORK,
                  RWORK,
                  RESULT(16));

              // Test 17:  Check norm of r'*A

              RESULT[17] = ZERO;
              if (M > CRANK.value) {
                RESULT[17] = zqrt17(
                    'No transpose',
                    1,
                    M,
                    N,
                    NRHS,
                    COPYA.asMatrix(),
                    LDA,
                    B.asMatrix(),
                    LDB,
                    COPYB.asMatrix(),
                    LDB,
                    C.asMatrix(),
                    WORK,
                    LWORK);
              }

              // Test 18:  Check if x is in the rowspace of A

              RESULT[18] = ZERO;
              if (N > CRANK.value) {
                RESULT[18] = zqrt14('No transpose', M, N, NRHS,
                    COPYA.asMatrix(), LDA, B.asMatrix(), LDB, WORK, LWORK);
              }

              // Print information about the tests that did not
              // pass the threshold.

              for (var K = 7; K <= 18; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                  NOUT.println(
                      ' M=${M.i5}, N=${N.i5}, NRHS=${NRHS.i4}, NB=${NB.i4}, type${ITYPE.i2}, test(${K.i2})=${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += 12;
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}

extension on Nout {
  void print9999(String trans, int m, int n, int nrhs, int nb, int type,
      int test, double ratio) {
    println(
        ' TRANS=\'${trans.a1}\', M=${m.i5}, N=${n.i5}, NRHS=${nrhs.i4}, NB=${nb.i4}, type${type.i2}, test(${test.i2})=${ratio.g12_5}');
  }
}
