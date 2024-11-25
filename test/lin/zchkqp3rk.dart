// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgeqp3rk.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zunmqr.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dlaord.dart';
import 'icopy.dart';
import 'xlaenv.dart';
import 'zlatb4.dart';
import 'zqpt01.dart';
import 'zqrt11.dart';
import 'zqrt12.dart';

void zchkqp3rk(
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
  final Array<Complex> A_,
  final Array<Complex> COPYA_,
  final Array<Complex> B_,
  final Array<Complex> COPYB_,
  final Array<double> S_,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
) {
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
  final S = S_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const NTYPES = 19, NTESTS = 5;
  const ZERO = 0.0, BIGNUM = 1.0e+38;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS), RDUMMY = Array<double>(1);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}QK';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }
  // final EPS = dlamch('Epsilon');
  infoc.INFOT = 0;

  for (var IM = 1; IM <= NM; IM++) {
    // Do for each value of M in MVAL.

    final M = MVAL[IM];
    final LDA = max(1, M).toInt();

    for (var IN = 1; IN <= NN; IN++) {
      // Do for each value of N in NVAL.

      final N = NVAL[IN];
      final MINMN = min(M, N);
      final LWORK = max(
          1,
          max(M * max(M, N).toInt() + 4 * MINMN + max(M, N).toInt(),
              M * N + 2 * MINMN + 4 * N));

      for (var INS = 1; INS <= NNS; INS++) {
        final NRHS = NSVAL[INS];

        // Set up parameters with ZLATB4 and generate
        // M-by-NRHS B matrix with ZLATMS.
        // IMAT = 14:
        // Random matrix, CNDNUM = 2, NORM = ONE,
        // MODE = 3 (geometric distribution of singular values).

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
            zlatb4(PATH, 14, M, NRHS);

        srnamc.SRNAMT = 'ZLATMS';
        zlatms(M, NRHS, DIST, ISEED, TYPE, S, MODE, CNDNUM, ANORM, KL, KU,
            'No packing', COPYB.asMatrix(), LDA, WORK, INFO);

        // Check error code from ZLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZLATMS', INFO.value, 0, ' ', M, NRHS, -1, -1, -1, 6,
              NFAIL, NERRS, NOUT);
          continue;
        }

        for (var IMAT = 1; IMAT <= NTYPES; IMAT++) {
          // Do the tests only if DOTYPE( IMAT ) is true.

          if (!DOTYPE[IMAT]) continue;

          // The type of distribution used to generate the random
          // eigen-/singular values:
          // ( 'S' for symmetric distribution ) => UNIFORM( -1, 1 )

          // Do for each type of NON-SYMMETRIC matrix:                               CNDNUM                     NORM                                     MODE
          //  1. Zero matrix
          //  2. Random, Diagonal, CNDNUM = 2                                        CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          //  3. Random, Upper triangular, CNDNUM = 2                                CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          //  4. Random, Lower triangular, CNDNUM = 2                                CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          //  5. Random, First column is zero, CNDNUM = 2                            CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          //  6. Random, Last MINMN column is zero, CNDNUM = 2                       CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          //  7. Random, Last N column is zero, CNDNUM = 2                           CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          //  8. Random, Middle column in MINMN is zero, CNDNUM = 2                  CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          //  9. Random, First half of MINMN columns are zero, CNDNUM = 2            CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          // 10. Random, Last columns are zero starting from MINMN/2+1, CNDNUM = 2   CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          // 11. Random, Half MINMN columns in the middle are zero starting
          //        from  MINMN/2-(MINMN/2)/2+1, CNDNUM = 2                          CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          // 12. Random, Odd columns are ZERO, CNDNUM = 2                            CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          // 13. Random, Even columns are ZERO, CNDNUM = 2                           CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          // 14. Random, CNDNUM = 2                                                  CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
          // 15. Random, CNDNUM = sqrt(0.1/EPS)                                      CNDNUM = BADC1 = sqrt(0.1/EPS)  ONE                                      3 ( geometric distribution of singular values )
          // 16. Random, CNDNUM = 0.1/EPS                                            CNDNUM = BADC2 = 0.1/EPS        ONE                                      3 ( geometric distribution of singular values )
          // 17. Random, CNDNUM = 0.1/EPS,                                           CNDNUM = BADC2 = 0.1/EPS        ONE                                      2 ( one small singular value, S(N)=1/CNDNUM )
          //       one small singular value S(N)=1/CNDNUM
          // 18. Random, CNDNUM = 2, scaled near underflow                           CNDNUM = 2                      SMALL = SAFMIN
          // 19. Random, CNDNUM = 2, scaled near overflow                            CNDNUM = 2                      LARGE = 1.0/( 0.25 * ( SAFMIN / EPS ) )  3 ( geometric distribution of singular values )

          if (IMAT == 1) {
            // Matrix 1: Zero matrix

            zlaset('Full', M, N, Complex.zero, Complex.zero, COPYA.asMatrix(),
                LDA);
            for (var I = 1; I <= MINMN; I++) {
              S[I] = ZERO;
            }
          } else if ((IMAT >= 2 && IMAT <= 4) || (IMAT >= 14 && IMAT <= 19)) {
            // Matrices 2-5.

            // Set up parameters with DLATB4 and generate a test
            // matrix with ZLATMS.

            final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
                zlatb4(PATH, IMAT, M, N);

            srnamc.SRNAMT = 'ZLATMS';
            zlatms(M, N, DIST, ISEED, TYPE, S, MODE, CNDNUM, ANORM, KL, KU,
                'No packing', COPYA.asMatrix(), LDA, WORK, INFO);

            // Check error code from ZLATMS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZLATMS', INFO.value, 0, ' ', M, N, -1, -1, -1, IMAT,
                  NFAIL, NERRS, NOUT);
              continue;
            }

            dlaord('Decreasing', MINMN, S, 1);
          } else if (MINMN >= 2 && IMAT >= 5 && IMAT <= 13) {
            // Rectangular matrices 5-13 that contain zero columns,
            // only for matrices MINMN >=2.

            // JB_ZERO is the column index of ZERO block.
            // NB_ZERO is the column block size of ZERO block.
            // NB_GEN is the column blcok size of the
            // generated block.
            // J_INC in the non_zero column index increment
            // for matrix 12 and 13.
            // J_FIRS_NZ is the index of the first non-zero
            // column.

            var JB_ZERO = 0, NB_ZERO = 0, NB_GEN = 0, J_INC = 0, J_FIRST_NZ = 0;

            if (IMAT == 5) {
              // First column is zero.

              JB_ZERO = 1;
              NB_ZERO = 1;
              NB_GEN = N - NB_ZERO;
            } else if (IMAT == 6) {
              // Last column MINMN is zero.

              JB_ZERO = MINMN;
              NB_ZERO = 1;
              NB_GEN = N - NB_ZERO;
            } else if (IMAT == 7) {
              // Last column N is zero.

              JB_ZERO = N;
              NB_ZERO = 1;
              NB_GEN = N - NB_ZERO;
            } else if (IMAT == 8) {
              // Middle column in MINMN is zero.

              JB_ZERO = MINMN ~/ 2 + 1;
              NB_ZERO = 1;
              NB_GEN = N - NB_ZERO;
            } else if (IMAT == 9) {
              // First half of MINMN columns is zero.

              JB_ZERO = 1;
              NB_ZERO = MINMN ~/ 2;
              NB_GEN = N - NB_ZERO;
            } else if (IMAT == 10) {
              // Last columns are zero columns,
              // starting from (MINMN / 2 + 1) column.

              JB_ZERO = MINMN ~/ 2 + 1;
              NB_ZERO = N - JB_ZERO + 1;
              NB_GEN = N - NB_ZERO;
            } else if (IMAT == 11) {
              // Half of the columns in the middle of MINMN
              // columns is zero, starting from
              // MINMN/2 - (MINMN/2)/2 + 1 column.

              JB_ZERO = MINMN ~/ 2 - (MINMN ~/ 2) ~/ 2 + 1;
              NB_ZERO = MINMN ~/ 2;
              NB_GEN = N - NB_ZERO;
            } else if (IMAT == 12) {
              // Odd-numbered columns are zero,

              NB_GEN = N ~/ 2;
              NB_ZERO = N - NB_GEN;
              J_INC = 2;
              J_FIRST_NZ = 2;
            } else if (IMAT == 13) {
              // Even-numbered columns are zero.

              NB_ZERO = N ~/ 2;
              NB_GEN = N - NB_ZERO;
              J_INC = 2;
              J_FIRST_NZ = 1;
            } else {
              throw UnimplementedError();
            }

            // 1) Set the first NB_ZERO columns in COPYA(1:M,1:N)
            //    to zero.

            zlaset('Full', M, NB_ZERO, Complex.zero, Complex.zero,
                COPYA.asMatrix(), LDA);

            // 2) Generate an M-by-(N-NB_ZERO) matrix with the
            //    chosen singular value distribution
            //    in COPYA(1:M,NB_ZERO+1:N).

            final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
                zlatb4(PATH, IMAT, M, NB_GEN);

            srnamc.SRNAMT = 'ZLATMS';

            final IND_OFFSET_GEN = NB_ZERO * LDA;

            zlatms(
                M,
                NB_GEN,
                DIST,
                ISEED,
                TYPE,
                S,
                MODE,
                CNDNUM,
                ANORM,
                KL,
                KU,
                'No packing',
                COPYA(IND_OFFSET_GEN + 1).asMatrix(),
                LDA,
                WORK,
                INFO);

            // Check error code from ZLATMS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZLATMS', INFO.value, 0, ' ', M, NB_GEN, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
              continue;
            }

            // 3) Swap the gererated colums from the right side
            // NB_GEN-size block in COPYA into correct column
            // positions.

            if (IMAT == 6 ||
                IMAT == 7 ||
                IMAT == 8 ||
                IMAT == 10 ||
                IMAT == 11) {
              // Move by swapping the generated columns
              // from the right NB_GEN-size block from
              // (NB_ZERO+1:NB_ZERO+JB_ZERO)
              // into columns (1:JB_ZERO-1).

              for (var J = 1; J <= JB_ZERO - 1; J++) {
                zswap(M, COPYA((NB_ZERO + J - 1) * LDA + 1), 1,
                    COPYA((J - 1) * LDA + 1), 1);
              }
            } else if (IMAT == 12 || IMAT == 13) {
              // ( IMAT = 12, Odd-numbered ZERO columns. )
              // Swap the generated columns from the right
              // NB_GEN-size block into the even zero colums in the
              // left NB_ZERO-size block.

              // ( IMAT = 13, Even-numbered ZERO columns. )
              // Swap the generated columns from the right
              // NB_GEN-size block into the odd zero colums in the
              // left NB_ZERO-size block.

              for (var J = 1; J <= NB_GEN; J += 1) {
                final IND_OUT = (NB_ZERO + J - 1) * LDA + 1;
                final IND_IN = (J_INC * (J - 1) + (J_FIRST_NZ - 1)) * LDA + 1;
                zswap(M, COPYA(IND_OUT), 1, COPYA(IND_IN), 1);
              }
            }

            // 5) Order the singular values generated by
            //    DLAMTS in decreasing order and add trailing zeros
            //    that correspond to zero columns.
            //    The total number of singular values is MINMN.

            final MINMNB_GEN = min(M, NB_GEN);

            dlaord('Decreasing', MINMNB_GEN, S, 1);

            for (var I = MINMNB_GEN + 1; I <= MINMN; I++) {
              S[I] = ZERO;
            }
          } else {
            // IF(MINMN < 2) skip this size for this matrix type.

            continue;
          }

          // Initialize a copy array for a pivot array for DGEQP3RK.

          for (var I = 1; I <= N; I++) {
            IWORK[I] = 0;
          }

          for (var INB = 1; INB <= NNB; INB++) {
            // Do for each pair of values (NB,NX) in NBVAL and NXVAL.

            final NB = NBVAL[INB];
            xlaenv(1, NB);
            final NX = NXVAL[INB];
            xlaenv(3, NX);

            // We do min(M,N)+1 because we need a test for KMAX > N,
            // when KMAX is larger than min(M,N), KMAX should be
            // KMAX = min(M,N)

            for (var KMAX = 0; KMAX <= min(M, N) + 1; KMAX++) {
              // Get a working copy of COPYA into A( 1:M,1:N ).
              // Get a working copy of COPYB into A( 1:M, (N+1):NRHS ).
              // Get a working copy of COPYB into into B( 1:M, 1:NRHS ).
              // Get a working copy of IWORK(1:N) awith zeroes into
              // which is going to be used as pivot array IWORK( N+1:2N ).
              // NOTE: IWORK(2N+1:3N) is going to be used as a WORK array
              // for the routine.

              zlacpy('All', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
              zlacpy('All', M, NRHS, COPYB.asMatrix(), LDA,
                  A(LDA * N + 1).asMatrix(), LDA);
              zlacpy('All', M, NRHS, COPYB.asMatrix(), LDA, B.asMatrix(), LDA);
              icopy(N, IWORK(1), 1, IWORK(N + 1), 1);

              // Compute the QR factorization with pivoting of A

              final LW =
                  max(1, max(2 * N + NB * (N + NRHS + 1), 3 * N + NRHS - 1));

              // Compute ZGEQP3RK factorization of A.

              final ABSTOL = Box(-1.0),
                  RELTOL = Box(-1.0),
                  MAXC2NRMK = Box(ZERO),
                  RELMAXC2NRMK = Box(ZERO);
              final KFACT = Box(0);
              srnamc.SRNAMT = 'ZGEQP3RK';
              zgeqp3rk(
                  M,
                  N,
                  NRHS,
                  KMAX,
                  ABSTOL,
                  RELTOL,
                  A.asMatrix(),
                  LDA,
                  KFACT,
                  MAXC2NRMK,
                  RELMAXC2NRMK,
                  IWORK(N + 1),
                  TAU,
                  WORK,
                  LW,
                  RWORK,
                  IWORK(2 * N + 1),
                  INFO);

              // Check error code from ZGEQP3RK.

              if (INFO.value < 0) {
                alaerh(PATH, 'ZGEQP3RK', INFO.value, 0, ' ', M, N, NX, -1, NB,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              if (KFACT.value == MINMN) {
                // Compute test 1:

                // This test in only for the full rank factorization of
                // the matrix A.

                // Array S(1:min(M,N)) contains svd(A) the sigular values
                // of the original matrix A in decreasing absolute value
                // order. The test computes svd(R), the vector sigular
                // values of the upper trapezoid of A(1:M,1:N) that
                // contains the factor R, in decreasing order. The test
                // returns the ratio:

                // 2-norm(svd(R) - svd(A)) / ( max(M,N) * 2-norm(svd(A)) * EPS )

                RESULT[1] =
                    zqrt12(M, N, A.asMatrix(), LDA, S, WORK, LWORK, RWORK);

                for (var T = 1; T <= 1; T++) {
                  if (RESULT[T] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                    NOUT.print9999('ZGEQP3RK', M, N, NRHS, KMAX, ABSTOL.value,
                        RELTOL.value, NB, NX, IMAT, T, RESULT[T]);
                    NFAIL++;
                  }
                }
                NRUN++;

                // End test 1
              }

              // Compute test 2:

              // The test returns the ratio:

              // 1-norm( A*P - Q*R ) / ( max(M,N) * 1-norm(A) * EPS )

              RESULT[2] = zqpt01(M, N, KFACT.value, COPYA.asMatrix(),
                  A.asMatrix(), LDA, TAU, IWORK(N + 1), WORK, LWORK);

              // Compute test 3:

              // The test returns the ratio:

              // 1-norm( Q**T * Q - I ) / ( M * EPS )

              RESULT[3] =
                  zqrt11(M, KFACT.value, A.asMatrix(), LDA, TAU, WORK, LWORK);

              // Print information about the tests that did not pass
              // the threshold.

              for (var T = 2; T <= 3; T++) {
                if (RESULT[T] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                  NOUT.print9999('ZGEQP3RK', M, N, NRHS, KMAX, ABSTOL.value,
                      RELTOL.value, NB, NX, IMAT, T, RESULT[T]);
                  NFAIL++;
                }
              }
              NRUN += 2;

              // Compute test 4:

              // This test is only for the factorizations with the
              // rank greater than 2.
              // The elements on the diagonal of R should be non-
              // increasing.

              // The test returns the ratio:

              // Returns 1.0e+100 if abs(R(K+1,K+1)) > abs(R(K,K)),
              // K=1:KFACT-1

              if (min(KFACT.value, MINMN) >= 2) {
                for (var J = 1; J <= KFACT.value - 1; J++) {
                  final DTEMP =
                      ((A[(J - 1) * M + J].abs() - A[J * M + J + 1].abs()) /
                          A[1].abs());

                  if (DTEMP < ZERO) {
                    RESULT[4] = BIGNUM;
                  }
                }

                // Print information about the tests that did not
                // pass the threshold.

                for (var T = 4; T <= 4; T++) {
                  if (RESULT[T] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                    NOUT.print9999('ZGEQP3RK', M, N, NRHS, KMAX, ABSTOL.value,
                        RELTOL.value, NB, NX, IMAT, T, RESULT[T]);
                    NFAIL++;
                  }
                }
                NRUN++;

                // End test 4.
              }

              // Compute test 5:

              // This test in only for matrix A with min(M,N) > 0.

              // The test returns the ratio:

              // 1-norm(Q**T * B - Q**T * B ) /
              //       ( M * EPS )

              // (1) Compute B:=Q**T * B in the matrix B.

              if (MINMN > 0) {
                final LWORK_MQR = max(1, NRHS);
                zunmqr(
                    'Left',
                    'Conjugate transpose',
                    M,
                    NRHS,
                    KFACT.value,
                    A.asMatrix(),
                    LDA,
                    TAU,
                    B.asMatrix(),
                    LDA,
                    WORK,
                    LWORK_MQR,
                    INFO);

                for (var I = 1; I <= NRHS; I++) {
                  // Compare N+J-th column of A and J-column of B.

                  zaxpy(M, -Complex.one, A((N + I - 1) * LDA + 1), 1,
                      B((I - 1) * LDA + 1), 1);
                }

                RESULT[5] =
                    (zlange('One-norm', M, NRHS, B.asMatrix(), LDA, RDUMMY) /
                            (M * dlamch('Epsilon')))
                        .abs();

                // Print information about the tests that did not pass
                // the threshold.

                for (var T = 5; T <= 5; T++) {
                  if (RESULT[T] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                    NOUT.print9999('ZGEQP3RK', M, N, NRHS, KMAX, ABSTOL.value,
                        RELTOL.value, NB, NX, IMAT, T, RESULT[T]);
                    NFAIL++;
                  }
                }
                NRUN++;

                // End compute test 5.
              }

              // END DO KMAX = 1, min(M,N)+1
            }

            // END DO for INB = 1, NNB
          }

          // END DO  for IMAT = 1, NTYPES
        }

        // END DO for INS = 1, NNS
      }

      // END DO for IN = 1, NN
    }

    // END DO for IM = 1, NM
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}

extension on Nout {
  void print9999(String s, int m, int n, int nrhs, int kmax, double abstol,
      double reltol, int nb, int nx, int type, int test, double ratio) {
    println(
        ' $s M =${m.i5}, N =${n.i5}, NRHS =${nrhs.i5}, KMAX =${kmax.i5}, ABSTOL =${abstol.g12_5}, RELTOL =${reltol.g12_5}, NB =${nb.i4}, NX =${nx.i4}, type ${type.i2}, test ${test.i2}, ratio =${ratio.g12_5}');
  }
}
