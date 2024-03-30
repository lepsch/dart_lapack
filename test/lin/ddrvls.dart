import 'dart:math';

import 'package:collection/collection.dart';
import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgels.dart';
import 'package:lapack/src/dgelsd.dart';
import 'package:lapack/src/dgelss.dart';
import 'package:lapack/src/dgelst.dart';
import 'package:lapack/src/dgelsy.dart';
import 'package:lapack/src/dgetsls.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alaerh.dart';
import 'alahd.dart';
import 'alasvm.dart';
import 'common.dart';
import 'derrls.dart';
import 'dqrt12.dart';
import 'dqrt13.dart';
import 'dqrt14.dart';
import 'dqrt15.dart';
import 'dqrt16.dart';
import 'dqrt17.dart';
import 'xlaenv.dart';

void ddrvls(
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
  final Array<double> A_,
  final Array<double> COPYA_,
  final Array<double> B_,
  final Array<double> COPYB_,
  final Array<double> C_,
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
  const NTESTS = 18;
  const SMLSIZ = 25;
  const ONE = 1.0, ZERO = 0.0;
  final ISEED = Array<int>(4), IWQ = Array<int>(1);
  final RESULT = Array<double>(NTESTS), WQ = Array<double>(1);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}LS';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }
  final EPS = dlamch('Epsilon');

  // Threshold for rank estimation

  final RCOND = Box(sqrt(EPS) - (sqrt(EPS) - EPS) / 2);

  // Test the error exits

  xlaenv(2, 2);
  xlaenv(9, SMLSIZ);
  if (TSTERR) derrls(PATH, NOUT);

  // Print the header if NM = 0 or NN = 0 and THRESH = 0.

  if ((NM == 0 || NN == 0) && THRESH == ZERO) alahd(NOUT, PATH);
  infoc.INFOT = 0;
  xlaenv(2, 2);
  xlaenv(9, SMLSIZ);

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
  // dqrt14, dqrt17 (two side cases), DQRT15 and dqrt12

  var LWORK = [
    1,
    (M + N) * NRHS,
    (N + NRHS) * (M + 2),
    (M + NRHS) * (N + 2),
    M + MNMIN,
    NRHS * MNMIN,
    2 * N + M,
    M * N + 4 * MNMIN + max(M, N),
    M * N + 2 * MNMIN + 4 * N
  ].max.toInt();
  var LIWORK = 1;

  // Iterate through all test cases and compute necessary workspace
  // sizes for ?GELS, ?GELST, ?GETSLS, ?GELSY, ?GELSS and ?GELSD
  // routines.

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];
    final LDA = max(1, M);
    for (var IN = 1; IN <= NN; IN++) {
      final N = NVAL[IN];
      //  final MNMIN = max(min( M, N ),1);
      final LDB = max(1, max(M, N));
      for (var INS = 1; INS <= NNS; INS++) {
        final NRHS = NSVAL[INS];
        var LWORK_DGELS = 0, LWORK_DGELST = 0, LWORK_DGETSLS = 0;
        for (var IRANK = 1; IRANK <= 2; IRANK++) {
          for (var ISCALE = 1; ISCALE <= 3; ISCALE++) {
            final ITYPE = (IRANK - 1) * 3 + ISCALE;
            if (DOTYPE[ITYPE]) {
              if (IRANK == 1) {
                for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
                  final TRANS = ITRAN == 1 ? 'N' : 'T';

                  // Compute workspace needed for DGELS
                  dgels(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB,
                      WQ, -1, INFO);
                  LWORK_DGELS = WQ[1].toInt();
                  // Compute workspace needed for DGELST
                  dgelst(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(),
                      LDB, WQ, -1, INFO);
                  LWORK_DGELST = WQ[1].toInt();
                  // Compute workspace needed for DGETSLS
                  dgetsls(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(),
                      LDB, WQ, -1, INFO);
                  LWORK_DGETSLS = WQ[1].toInt();
                }
              }
              final CRANK = Box(0);
              // Compute workspace needed for DGELSY
              dgelsy(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, IWQ,
                  RCOND.value, CRANK, WQ, -1, INFO);
              final LWORK_DGELSY = WQ[1].toInt();
              // Compute workspace needed for DGELSS
              dgelss(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, S, RCOND,
                  CRANK, WQ, -1, INFO);
              final LWORK_DGELSS = WQ[1].toInt();
              // Compute workspace needed for DGELSD
              dgelsd(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, S,
                  RCOND.value, CRANK, WQ, -1, IWQ, INFO);
              final LWORK_DGELSD = WQ[1].toInt();
              // Compute LIWORK workspace needed for DGELSY and DGELSD
              LIWORK = max(LIWORK, max(N, IWQ[1]));
              // Compute LWORK workspace needed for all functions
              LWORK = [
                LWORK,
                LWORK_DGELS,
                LWORK_DGELST,
                LWORK_DGETSLS,
                LWORK_DGELSY,
                LWORK_DGELSS,
                LWORK_DGELSD
              ].max.toInt();
            }
          }
        }
      }
    }
  }

  final LWLSY = LWORK;

  final WORK = Array<double>(LWORK);
  final IWORK = Array<int>(LIWORK);

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];
    final LDA = max(1, M);

    for (var IN = 1; IN <= NN; IN++) {
      final N = NVAL[IN];
      final MNMIN = max(min(M, N), 1);
      final LDB = max(1, max(M, N));
      // final MB = MNMIN + 1;

      for (var INS = 1; INS <= NNS; INS++) {
        final NRHS = NSVAL[INS];

        for (var IRANK = 1; IRANK <= 2; IRANK++) {
          for (var ISCALE = 1; ISCALE <= 3; ISCALE++) {
            final ITYPE = (IRANK - 1) * 3 + ISCALE;
            if (!DOTYPE[ITYPE]) continue;
            // =====================================================
            //       Begin test DGELS
            // =====================================================
            final NORMA = Box(ZERO);
            if (IRANK == 1) {
              // Generate a matrix of scaling type ISCALE

              dqrt13(ISCALE, M, N, COPYA.asMatrix(), LDA, NORMA, ISEED);

              // Loop for testing different block sizes.

              for (var INB = 1; INB <= NNB; INB++) {
                final NB = NBVAL[INB];
                xlaenv(1, NB);
                xlaenv(3, NXVAL[INB]);

                // Loop for testing non-transposed and transposed.

                for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
                  final (TRANS, NROWS, NCOLS) =
                      ITRAN == 1 ? ('N', M, N) : ('T', N, M);
                  final LDWORK = max(1, NCOLS);

                  // Set up a consistent rhs

                  if (NCOLS > 0) {
                    dlarnv(2, ISEED, NCOLS * NRHS, WORK);
                    dscal(NCOLS * NRHS, ONE / NCOLS, WORK, 1);
                  }
                  dgemm(
                      TRANS,
                      'No transpose',
                      NROWS,
                      NRHS,
                      NCOLS,
                      ONE,
                      COPYA.asMatrix(),
                      LDA,
                      WORK.asMatrix(),
                      LDWORK,
                      ZERO,
                      B.asMatrix(),
                      LDB);
                  dlacpy('Full', NROWS, NRHS, B.asMatrix(), LDB,
                      COPYB.asMatrix(), LDB);

                  // Solve LS or overdetermined system

                  if (M > 0 && N > 0) {
                    dlacpy(
                        'Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
                    dlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                        B.asMatrix(), LDB);
                  }
                  srnamc.SRNAMT = 'DGELS ';
                  dgels(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB,
                      WORK, LWORK, INFO);
                  if (INFO.value != 0) {
                    alaerh(PATH, 'DGELS ', INFO.value, 0, TRANS, M, N, NRHS, -1,
                        NB, ITYPE, NFAIL, NERRS, NOUT);
                  }

                  // Test 1: Check correctness of results
                  // for DGELS, compute the residual:
                  // RESID = norm(B - A*X) /
                  // / ( max(m,n) * norm(A) * norm(X) * EPS )

                  if (NROWS > 0 && NRHS > 0) {
                    dlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                        C.asMatrix(), LDB);
                  }
                  dqrt16(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA, B.asMatrix(),
                      LDB, C.asMatrix(), LDB, WORK, RESULT(1));

                  // Test 2: Check correctness of results
                  // for DGELS.

                  if ((ITRAN == 1 && M >= N) || (ITRAN == 2 && M < N)) {
                    // Solving LS system, compute:
                    // r = norm((B- A*X)**T * A) /
                    //  / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                    RESULT[2] = dqrt17(
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

                    RESULT[2] = dqrt14(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA,
                        B.asMatrix(), LDB, WORK, LWORK);
                  }

                  // Print information about the tests that
                  // did not pass the threshold.

                  for (var K = 1; K <= 2; K++) {
                    if (RESULT[K] >= THRESH) {
                      if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                      NOUT.println(
                          ' TRANS=\'${TRANS.a1}\', M=${M.i5}, N=${N.i5}, NRHS=${NRHS.i4}, NB=${NB.i4}, type${ITYPE.i2}, test(${K.i2})=${RESULT[K].g12_5}');
                      NFAIL++;
                    }
                  }
                  NRUN += 2;
                }
              }
            }
            // =====================================================
            //       End test DGELS
            // =====================================================
            // =====================================================
            //       Begin test DGELST
            // =====================================================
            if (IRANK == 1) {
              // Generate a matrix of scaling type ISCALE

              dqrt13(ISCALE, M, N, COPYA.asMatrix(), LDA, NORMA, ISEED);

              // Loop for testing different block sizes.

              for (var INB = 1; INB <= NNB; INB++) {
                final NB = NBVAL[INB];
                xlaenv(1, NB);

                // Loop for testing non-transposed and transposed.

                for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
                  final (TRANS, NROWS, NCOLS) =
                      ITRAN == 1 ? ('N', M, N) : ('T', N, M);
                  final LDWORK = max(1, NCOLS);

                  // Set up a consistent rhs

                  if (NCOLS > 0) {
                    dlarnv(2, ISEED, NCOLS * NRHS, WORK);
                    dscal(NCOLS * NRHS, ONE / NCOLS, WORK, 1);
                  }
                  dgemm(
                      TRANS,
                      'No transpose',
                      NROWS,
                      NRHS,
                      NCOLS,
                      ONE,
                      COPYA.asMatrix(),
                      LDA,
                      WORK.asMatrix(),
                      LDWORK,
                      ZERO,
                      B.asMatrix(),
                      LDB);
                  dlacpy('Full', NROWS, NRHS, B.asMatrix(), LDB,
                      COPYB.asMatrix(), LDB);

                  // Solve LS or overdetermined system

                  if (M > 0 && N > 0) {
                    dlacpy(
                        'Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
                    dlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                        B.asMatrix(), LDB);
                  }
                  srnamc.SRNAMT = 'DGELST';
                  dgelst(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(),
                      LDB, WORK, LWORK, INFO);
                  if (INFO.value != 0) {
                    alaerh(PATH, 'DGELST', INFO.value, 0, TRANS, M, N, NRHS, -1,
                        NB, ITYPE, NFAIL, NERRS, NOUT);
                  }

                  // Test 3: Check correctness of results
                  // for DGELST, compute the residual:
                  // RESID = norm(B - A*X) /
                  // / ( max(m,n) * norm(A) * norm(X) * EPS )

                  if (NROWS > 0 && NRHS > 0) {
                    dlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                        C.asMatrix(), LDB);
                  }
                  dqrt16(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA, B.asMatrix(),
                      LDB, C.asMatrix(), LDB, WORK, RESULT(3));

                  // Test 4: Check correctness of results
                  // for DGELST.

                  if ((ITRAN == 1 && M >= N) || (ITRAN == 2 && M < N)) {
                    // Solving LS system, compute:
                    // r = norm((B- A*X)**T * A) /
                    //  / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                    RESULT[4] = dqrt17(
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

                    RESULT[4] = dqrt14(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA,
                        B.asMatrix(), LDB, WORK, LWORK);
                  }

                  // Print information about the tests that
                  // did not pass the threshold.

                  for (var K = 3; K <= 4; K++) {
                    if (RESULT[K] >= THRESH) {
                      if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                      NOUT.println(
                          ' TRANS=\'${TRANS.a1}\', M=${M.i5}, N=${N.i5}, NRHS=${NRHS.i4}, NB=${NB.i4}, type${ITYPE.i2}, test(${K.i2})=${RESULT[K].g12_5}');
                      NFAIL++;
                    }
                  }
                  NRUN += 2;
                }
              }
            }
            // =====================================================
            //       End test DGELST
            // =====================================================
            // =====================================================
            //       Begin test DGETSLS
            // =====================================================
            if (IRANK == 1) {
              // Generate a matrix of scaling type ISCALE

              dqrt13(ISCALE, M, N, COPYA.asMatrix(), LDA, NORMA, ISEED);

              // Loop for testing different block sizes MB.

              for (var IMB = 1; IMB <= NNB; IMB++) {
                final MB = NBVAL[IMB];
                xlaenv(1, MB);

                // Loop for testing different block sizes NB.

                for (var INB = 1; INB <= NNB; INB++) {
                  final NB = NBVAL[INB];
                  xlaenv(2, NB);

                  // Loop for testing non-transposed
                  // and transposed.

                  for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
                    final (TRANS, NROWS, NCOLS) =
                        ITRAN == 1 ? ('N', M, N) : ('T', N, M);
                    final LDWORK = max(1, NCOLS);

                    // Set up a consistent rhs

                    if (NCOLS > 0) {
                      dlarnv(2, ISEED, NCOLS * NRHS, WORK);
                      dscal(NCOLS * NRHS, ONE / NCOLS, WORK, 1);
                    }
                    dgemm(
                        TRANS,
                        'No transpose',
                        NROWS,
                        NRHS,
                        NCOLS,
                        ONE,
                        COPYA.asMatrix(),
                        LDA,
                        WORK.asMatrix(),
                        LDWORK,
                        ZERO,
                        B.asMatrix(),
                        LDB);
                    dlacpy('Full', NROWS, NRHS, B.asMatrix(), LDB,
                        COPYB.asMatrix(), LDB);

                    // Solve LS or overdetermined system

                    if (M > 0 && N > 0) {
                      dlacpy('Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(),
                          LDA);
                      dlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                          B.asMatrix(), LDB);
                    }
                    srnamc.SRNAMT = 'DGETSLS';
                    dgetsls(TRANS, M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(),
                        LDB, WORK, LWORK, INFO);
                    if (INFO.value != 0) {
                      alaerh(PATH, 'DGETSLS', INFO.value, 0, TRANS, M, N, NRHS,
                          -1, NB, ITYPE, NFAIL, NERRS, NOUT);
                    }

                    // Test 5: Check correctness of results
                    // for DGETSLS, compute the residual:
                    // RESID = norm(B - A*X) /
                    // / ( max(m,n) * norm(A) * norm(X) * EPS )

                    if (NROWS > 0 && NRHS > 0) {
                      dlacpy('Full', NROWS, NRHS, COPYB.asMatrix(), LDB,
                          C.asMatrix(), LDB);
                    }
                    dqrt16(TRANS, M, N, NRHS, COPYA.asMatrix(), LDA,
                        B.asMatrix(), LDB, C.asMatrix(), LDB, WORK, RESULT(5));

                    // Test 6: Check correctness of results
                    // for DGETSLS.

                    if ((ITRAN == 1 && M >= N) || (ITRAN == 2 && M < N)) {
                      // Solving LS system, compute:
                      // r = norm((B- A*X)**T * A) /
                      // / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                      RESULT[6] = dqrt17(
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

                      RESULT[6] = dqrt14(TRANS, M, N, NRHS, COPYA.asMatrix(),
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
            //       End test DGETSLS
            // =====================================================

            // Generate a matrix of scaling type ISCALE and rank
            // type IRANK.

            final RANK = Box(0);
            final NORMB = Box(ZERO);
            dqrt15(
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

              // Test DGELSY

              // DGELSY:  Compute the minimum-norm solution X
              // to min( norm( A * X - B ) )
              // using the rank-revealing orthogonal
              // factorization.

              // Initialize vector IWORK.

              for (var J = 1; J <= N; J++) {
                IWORK[J] = 0;
              }

              dlacpy('Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
              dlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, B.asMatrix(), LDB);

              final CRANK = Box(0);
              srnamc.SRNAMT = 'DGELSY';
              dgelsy(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, IWORK,
                  RCOND.value, CRANK, WORK, LWLSY, INFO);
              if (INFO.value != 0) {
                alaerh(PATH, 'DGELSY', INFO.value, 0, ' ', M, N, NRHS, -1, NB,
                    ITYPE, NFAIL, NERRS, NOUT);
              }

              // Test 7:  Compute relative error in svd
              //          workspace: M*N + 4*min(M,N) + max(M,N)

              RESULT[7] = dqrt12(CRANK.value, CRANK.value, A.asMatrix(), LDA,
                  COPYS, WORK, LWORK);

              // Test 8:  Compute error in solution
              //          workspace:  M*NRHS + M

              dlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, WORK.asMatrix(),
                  LDWORK);
              dqrt16(
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
                  WORK(M * NRHS + 1),
                  RESULT(8));

              // Test 9:  Check norm of r'*A
              //          workspace: NRHS*(M+N)

              RESULT[9] = ZERO;
              if (M > CRANK.value) {
                RESULT[9] = dqrt17(
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
                RESULT[10] = dqrt14('No transpose', M, N, NRHS,
                    COPYA.asMatrix(), LDA, B.asMatrix(), LDB, WORK, LWORK);
              }

              // Test DGELSS

              // DGELSS:  Compute the minimum-norm solution X
              // to min( norm( A * X - B ) )
              // using the SVD.

              dlacpy('Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
              dlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, B.asMatrix(), LDB);
              srnamc.SRNAMT = 'DGELSS';
              dgelss(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, S, RCOND,
                  CRANK, WORK, LWORK, INFO);
              if (INFO.value != 0) {
                alaerh(PATH, 'DGELSS', INFO.value, 0, ' ', M, N, NRHS, -1, NB,
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

              dlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, WORK.asMatrix(),
                  LDWORK);
              dqrt16(
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
                  WORK(M * NRHS + 1),
                  RESULT(12));

              // Test 13:  Check norm of r'*A

              RESULT[13] = ZERO;
              if (M > CRANK.value) {
                RESULT[13] = dqrt17(
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
                RESULT[14] = dqrt14('No transpose', M, N, NRHS,
                    COPYA.asMatrix(), LDA, B.asMatrix(), LDB, WORK, LWORK);
              }

              // Test DGELSD

              // DGELSD:  Compute the minimum-norm solution X
              // to min( norm( A * X - B ) ) using a
              // divide and conquer SVD.

              // Initialize vector IWORK.

              for (var J = 1; J <= N; J++) {
                IWORK[J] = 0;
              }

              dlacpy('Full', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
              dlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, B.asMatrix(), LDB);

              srnamc.SRNAMT = 'DGELSD';
              dgelsd(M, N, NRHS, A.asMatrix(), LDA, B.asMatrix(), LDB, S,
                  RCOND.value, CRANK, WORK, LWORK, IWORK, INFO);
              if (INFO.value != 0) {
                alaerh(PATH, 'DGELSD', INFO.value, 0, ' ', M, N, NRHS, -1, NB,
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

              dlacpy('Full', M, NRHS, COPYB.asMatrix(), LDB, WORK.asMatrix(),
                  LDWORK);
              dqrt16(
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
                  WORK(M * NRHS + 1),
                  RESULT(16));

              // Test 17:  Check norm of r'*A

              RESULT[17] = ZERO;
              if (M > CRANK.value) {
                RESULT[17] = dqrt17(
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
                RESULT[18] = dqrt14('No transpose', M, N, NRHS,
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
