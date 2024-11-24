// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:io';
import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:lapack/src/zgedmd.dart';
import 'package:lapack/src/zgedmdq.dart';
import 'package:lapack/src/zgeev.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlarnv.dart';
import 'package:lapack/src/zlascl.dart';

import '../matgen/zlatmr.dart';

void main() async {
  final NIN = Nin(stdin), NOUT = Nout(stdout);
  const ONE = 1.0, ZERO = 0.0;
  final WDUMMY = Array<double>(2);
  final IDUMMY = Array<int>(4), ISEED = Array<int>(4);
  double ANORM,
      COND,
      CONDL,
      CONDR,
      EPS,
      TOL,
      TOL2,
      SVDIFF = 0,
      TMP,
      TMP_AU = 0,
      TMP_FQR = 0,
      TMP_REZ = 0,
      TMP_REZQ = 0,
      TMP_ZXW = 0,
      TMP_EX;
  Complex ZMAX;
  int LZWORK;
  final ZDUMMY = Array<Complex>(22), ZDUM2X2 = Matrix<Complex>(2, 2);
  int LDF, LDS, LDA, LDAU, LDW, LDX, LDY, LDZ, LIWORK, LWORK, NRNK;
  int NFAIL_AU,
      NFAIL_F_QR,
      NFAIL_REZ,
      NFAIL_REZQ,
      NFAIL_SVDIFF,
      NFAIL_TOTAL,
      NFAILQ_TOTAL,
      NFAIL_Z_XV,
      MODEL,
      MODER,
      WHTSVD;
  String GRADE, PIVTNG, RSIGN, WANTQ, WANTR;
  bool TEST_QRDMD;
  final INFO = Box(0), K = Box(0), KQ = Box(0);

  // The test is always in pairs : ( zgedmd and zgedmdq )
  // because the test includes comparing the results (in pairs).
  //
  TEST_QRDMD = true; // This code by default performs tests on zgedmdq
  // Since the QR factorizations based algorithm is designed for
  // single trajectory data, only single trajectory tests will
  // be performed with xGEDMDQ.
  WANTQ = 'Q';
  WANTR = 'R';

  EPS = dlamch('P'); // machine precision DP

  // Global counters of failures of some particular tests
  NFAIL_REZ = 0;
  NFAIL_REZQ = 0;
  NFAIL_Z_XV = 0;
  NFAIL_F_QR = 0;
  NFAIL_AU = 0;
  NFAIL_SVDIFF = 0;
  NFAIL_TOTAL = 0;
  NFAILQ_TOTAL = 0;

  for (final LLOOP in 1.through(4)) {
    NOUT.println('L Loop Index = $LLOOP');

    // Set the dimensions of the problem ...
    NOUT.println('M = ');
    final M = await NIN.readInt();
    NOUT.println('$M');
    // and the number of snapshots.
    NOUT.println('N = ');
    final N = await NIN.readInt();
    NOUT.println('$N');

    // Test the dimensions
    if ((min(M, N) == 0) || (M < N)) {
      NOUT.println('Bad dimensions. Required: M >= N > 0.');
      return;
    }
    //
    // The seed inside the LLOOP so that each pass can be reproduced easily.
    ISEED[1] = 4;
    ISEED[2] = 3;
    ISEED[3] = 2;
    ISEED[4] = 1;

    LDA = M;
    LDF = M;
    LDX = M;
    LDY = M;
    LDW = N;
    LDZ = M;
    LDAU = M;
    LDS = N;

    TMP_ZXW = ZERO;
    TMP_AU = ZERO;
    TMP_REZ = ZERO;
    TMP_REZQ = ZERO;
    SVDIFF = ZERO;
    TMP_EX = ZERO;

    // REAL(KIND=WP), ALLOCATABLE, DIMENSION(:)   ::   WORK
    // int      , ALLOCATABLE, DIMENSION(:)   :: IWORK
    // COMPLEX(KIND=WP), ALLOCATABLE, DIMENSION(:)   ::  ZEIGSA, ZWORK

    final ZA = Matrix<Complex>(LDA, M),
        ZAC = Matrix<Complex>(LDA, M),
        ZF = Matrix<Complex>(LDF, N + 1),
        ZF0 = Matrix<Complex>(LDF, N + 1),
        ZF1 = Matrix<Complex>(LDF, N + 1),
        ZX = Matrix<Complex>(LDX, N),
        ZX0 = Matrix<Complex>(LDX, N),
        ZY = Matrix<Complex>(LDY, N + 1),
        ZY0 = Matrix<Complex>(LDY, N + 1),
        ZY1 = Matrix<Complex>(LDY, N + 1),
        ZAU = Matrix<Complex>(LDAU, N),
        ZW = Matrix<Complex>(LDW, N),
        ZS = Matrix<Complex>(LDS, N),
        ZZ = Matrix<Complex>(LDZ, N),
        ZZ1 = Matrix<Complex>(LDZ, N);
    final ZEIGS = Array<Complex>(N);

    final RES = Array<double>(N),
        RES1 = Array<double>(N),
        RESEX = Array<double>(N),
        SINGVX = Array<double>(N),
        SINGVQX = Array<double>(N);

    TOL = M * EPS;
    // This mimics O(M*N)*EPS bound for accumulated roundoff error.
    // The factor 10 is somewhat arbitrary.
    TOL2 = 10 * M * N * EPS;

    //

    for (final K_TRAJ in 1.through(2)) {
      // Number of intial conditions in the simulation/trajectories (1 or 2)

      COND = 1.0e4;
      ZMAX = Complex(1.0e1, 1.0e1);
      RSIGN = 'F';
      GRADE = 'N';
      MODEL = 6;
      CONDL = 1.0e1;
      MODER = 6;
      CONDR = 1.0e1;
      PIVTNG = 'N';

      // Loop over all parameter MODE values for zlatmr (+1,..,+6)
      for (final MODE in 1.through(6)) {
        final IWORK = Array<int>(2 * M);
        final ZDA = Array<Complex>(M);
        final ZDL = Array<Complex>(M);
        final ZDR = Array<Complex>(M);

        zlatmr(
            M,
            M,
            'N',
            ISEED,
            'N',
            ZDA,
            MODE,
            COND,
            ZMAX,
            RSIGN,
            GRADE,
            ZDL,
            MODEL,
            CONDL,
            ZDR,
            MODER,
            CONDR,
            PIVTNG,
            IWORK,
            M,
            M,
            ZERO,
            -ONE,
            'N',
            ZA,
            LDA,
            IWORK(M + 1),
            INFO);
        // DEALLOCATE( ZDR )
        // DEALLOCATE( ZDL )
        // DEALLOCATE( ZDA )
        // DEALLOCATE( IWORK )

        LZWORK = max(1, 2 * M);
        final ZEIGSA = Array<Complex>(M);
        final ZWORK = Array<Complex>(LZWORK);
        final WORK = Array<double>(2 * M);

        for (final i in 1.through(M)) {
          for (final j in 1.through(M)) {
            ZAC[i][j] = ZA[i][j];
          }
        }

        zgeev('N', 'N', M, ZAC, LDA, ZEIGSA, ZDUM2X2, 2, ZDUM2X2, 2, ZWORK,
            LZWORK, WORK, INFO);
        // DEALLOCATE(WORK)
        // DEALLOCATE(ZWORK)

        TMP = ZEIGSA[izamax(M, ZEIGSA, 1)].abs(); // The spectral radius of ZA
        // Scale the matrix ZA to have unit spectral radius.
        zlascl('G', 0, 0, TMP, ONE, M, M, ZA, LDA, INFO);
        zlascl('G', 0, 0, TMP, ONE, M, 1, ZEIGSA.asMatrix(), M, INFO);
        ANORM = zlange('F', M, M, ZA, LDA, WDUMMY);

        if (K_TRAJ == 2) {
          // generate data as two trajectories
          // with two inital conditions
          zlarnv(2, ISEED, M, ZF(1, 1).asArray());
          for (final i in 1.through(N ~/ 2)) {
            zgemv('N', M, M, Complex.one, ZA, LDA, ZF(1, i).asArray(), 1,
                Complex.zero, ZF(1, i + 1).asArray(), 1);
          }

          for (final i in 1.through(M)) {
            for (final j in 1.through(N ~/ 2)) {
              ZX0[i][j] = ZF[i][j];
              ZY0[i][j] = ZF[i][j + 1];
            }
          }

          zlarnv(2, ISEED, M, ZF(1, 1).asArray());
          for (final i in 1.through(N - N ~/ 2)) {
            zgemv('N', M, M, Complex.one, ZA, LDA, ZF(1, i).asArray(), 1,
                Complex.zero, ZF(1, i + 1).asArray(), 1);
          }

          for (final i in 1.through(M)) {
            for (final j in 1.through(N - N ~/ 2)) {
              ZX0[i][j + N ~/ 2 + 1] = ZF[i][j];
              ZY0[i][j + N ~/ 2 + 1] = ZF[i][j + 1];
            }
          }
        } else {
          zlarnv(2, ISEED, M, ZF(1, 1).asArray());
          for (final i in 1.through(N)) {
            zgemv('N', M, M, Complex.one, ZA, M, ZF(1, i).asArray(), 1,
                Complex.zero, ZF(1, i + 1).asArray(), 1);
          }

          for (final i in 1.through(M)) {
            for (final j in 1.through(N + 1)) {
              ZF0[i][j] = ZF[i][j];
            }
          }

          for (final i in 1.through(M)) {
            for (final j in 1.through(N)) {
              ZX0[i][j] = ZF0[i][j];
              ZY0[i][j] = ZF0[i][j + 1];
            }
          }
        }

        // DEALLOCATE( ZEIGSA )
        //

        for (final iJOBZ in 1.through(4)) {
          final (JOBZ, RESIDS) = switch (iJOBZ) {
            1 => ('V', 'R'),
            2 => ('V', 'N'),
            3 => ('F', 'N'),
            4 => ('N', 'N'),
            _ => throw UnimplementedError(),
          };

          for (final iJOBREF in 1.through(3)) {
            final JOBREF = switch (iJOBREF) {
              1 => 'R',
              2 => 'E',
              3 => 'N',
              _ => throw UnimplementedError(),
            };

            for (final iSCALE in 1.through(4)) {
              final SCALE = switch (iSCALE) {
                1 => 'S',
                2 => 'C',
                3 => 'Y',
                4 => 'N',
                _ => throw UnimplementedError(),
              };

              for (final iNRNK in [-1 - 2]) {
                NRNK = iNRNK;

                for (final iWHTSVD in 1.through(3)) {
                  // Check all four options to compute the POD basis
                  // via the SVD.
                  WHTSVD = iWHTSVD;

                  for (final LWMINOPT in 1.through(2)) {
                    // Workspace query for the minimal (1) and for the optimal
                    // (2) workspace lengths determined by workspace query.

                    // zgedmd is always tested and its results are also used for
                    // comparisons with zgedmdq.

                    for (final i in 1.through(M)) {
                      for (final j in 1.through(N)) {
                        ZX[i][j] = ZX0[i][j];
                        ZY[i][j] = ZY0[i][j];
                      }
                    }

                    zgedmd(
                        SCALE,
                        JOBZ,
                        RESIDS,
                        JOBREF,
                        WHTSVD,
                        M,
                        N,
                        ZX,
                        LDX,
                        ZY,
                        LDY,
                        NRNK,
                        TOL,
                        K,
                        ZEIGS,
                        ZZ,
                        LDZ,
                        RES,
                        ZAU,
                        LDAU,
                        ZW,
                        LDW,
                        ZS,
                        LDS,
                        ZDUMMY,
                        -1,
                        WDUMMY,
                        -1,
                        IDUMMY,
                        -1,
                        INFO);
                    if ((INFO.value == 2) ||
                        (INFO.value == 3) ||
                        (INFO.value < 0)) {
                      NOUT.println(
                          'to zgedmd workspace query failed. Check the calling sequence and the code.');
                      NOUT.println('The error code is ${INFO.value}');
                      NOUT.println(
                          'The input parameters were $SCALE $JOBZ $RESIDS $JOBREF $WHTSVD $M $N $LDX $LDY $NRNK $TOL $LDZ $LDAU $LDW $LDS');
                      return;
                    }

                    LZWORK = ZDUMMY[LWMINOPT].toInt();
                    LWORK = WDUMMY[1].toInt();
                    LIWORK = IDUMMY[1];

                    final ZWORK = Array<Complex>(LZWORK);
                    final WORK = Array<double>(LWORK);
                    final IWORK = Array<int>(LIWORK);

                    zgedmd(
                        SCALE,
                        JOBZ,
                        RESIDS,
                        JOBREF,
                        WHTSVD,
                        M,
                        N,
                        ZX,
                        LDX,
                        ZY,
                        LDY,
                        NRNK,
                        TOL,
                        K,
                        ZEIGS,
                        ZZ,
                        LDZ,
                        RES,
                        ZAU,
                        LDAU,
                        ZW,
                        LDW,
                        ZS,
                        LDS,
                        ZWORK,
                        LZWORK,
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        INFO);

                    if (INFO.value != 0) {
                      NOUT.println(
                          'to zgedmd failed.            &Check the calling sequence and the code.');
                      NOUT.println('The error code is ${INFO.value}');
                      NOUT.println(
                          'The input parameters were $SCALE $JOBZ $RESIDS $JOBREF $WHTSVD $M $N $LDX $LDY $NRNK $TOL');
                      return;
                    }

                    for (final i in 1.through(N)) {
                      SINGVX[i] = WORK[i];
                    }

                    // zgedmd check point
                    if (lsame(JOBZ, 'V')) {
                      // Check that Z = X*W, on return from zgedmd
                      // This checks that the returned eigenvectors in Z are
                      // the product of the SVD'POD basis returned in X
                      // and the eigenvectors of the rayleigh quotient
                      // returned in W
                      zgemm('N', 'N', M, K.value, K.value, Complex.one, ZX, LDX,
                          ZW, LDW, Complex.zero, ZZ1, LDZ);
                      TMP = ZERO;
                      for (final i in 1.through(K.value)) {
                        zaxpy(M, -Complex.one, ZZ(1, i).asArray(), 1,
                            ZZ1(1, i).asArray(), 1);
                        TMP = max(TMP, dznrm2(M, ZZ1(1, i).asArray(), 1));
                      }
                      TMP_ZXW = max(TMP_ZXW, TMP);
                      if (TMP_ZXW <= 10 * M * EPS) {
                        // NOUT.println( ' :) .... OK .........zgedmd PASSED.');
                      } else {
                        NFAIL_Z_XV = NFAIL_Z_XV + 1;
                        NOUT.println(
                            ':( .................zgedmd FAILED!Check the code for implementation errors.');
                        NOUT.println(
                            'The input parameters were $SCALE $JOBZ $RESIDS $JOBREF $WHTSVD $M $N $LDX $LDY $NRNK $TOL');
                      }
                    }

                    // zgedmd check point
                    if (lsame(JOBREF, 'R')) {
                      // The matrix A*U is returned for computing refined Ritz vectors.
                      // Check that A*U is computed correctly using the formula
                      // A*U = Y * V * inv(SIGMA). This depends on the
                      // accuracy in the computed singular values and vectors of X.
                      // See the paper for an error analysis.
                      // Note that the left singular vectors of the input matrix X
                      // are returned in the array X.
                      zgemm('N', 'N', M, K.value, M, Complex.one, ZA, LDA, ZX,
                          LDX, Complex.zero, ZZ1, LDZ);
                      TMP = ZERO;
                      for (final i in 1.through(K.value)) {
                        zaxpy(M, -Complex.one, ZAU(1, i).asArray(), 1,
                            ZZ1(1, i).asArray(), 1);
                        TMP = max(
                            TMP,
                            dznrm2(M, ZZ1(1, i).asArray(), 1) *
                                SINGVX[K.value] /
                                (ANORM * SINGVX[1]));
                      }
                      TMP_AU = max(TMP_AU, TMP);
                      if (TMP <= TOL2) {
                        // NOUT.println( ':) .... OK .........zgedmd PASSED.');
                      } else {
                        NFAIL_AU = NFAIL_AU + 1;
                        NOUT.println(
                            ':( .................zgedmd FAILED!Check the code for implementation errors.');
                        NOUT.println(
                            'The input parameters were $SCALE $JOBZ $RESIDS $JOBREF $WHTSVD $M $N $LDX $LDY $NRNK $TOL');
                      }
                    } else if (lsame(JOBREF, 'E')) {
                      // The unscaled vectors of the Exact DMD are computed.
                      // This option is included for the sake of completeness,
                      // for users who prefer the Exact DMD vectors. The
                      // returned vectors are in the real form, in the same way
                      // as the Ritz vectors. Here we just save the vectors
                      // and test them separately using a Matlab script.

                      zgemm('N', 'N', M, K.value, M, Complex.one, ZA, LDA, ZAU,
                          LDAU, Complex.zero, ZY1, LDY);

                      for (final i in 1.through(K.value)) {
                        // have a real eigenvalue with real eigenvector
                        zaxpy(M, -ZEIGS[i], ZAU(1, i).asArray(), 1,
                            ZY1(1, i).asArray(), 1);
                        RESEX[i] = dznrm2(M, ZY1(1, i).asArray(), 1) /
                            dznrm2(M, ZAU(1, i).asArray(), 1);
                      }
                    }
                    // zgedmd check point

                    if (lsame(RESIDS, 'R')) {
                      // Compare the residuals returned by zgedmd with the
                      // explicitly computed residuals using the matrix A.
                      // Compute explicitly Y1 = A*Z
                      zgemm('N', 'N', M, K.value, M, Complex.one, ZA, LDA, ZZ,
                          LDZ, Complex.zero, ZY1, LDY);
                      // and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
                      // of the invariant subspaces that correspond to complex conjugate
                      // pairs of eigencalues. (See the description of Z in zgedmd,)

                      for (final i in 1.through(K.value)) {
                        // have a real eigenvalue with real eigenvector
                        zaxpy(M, -ZEIGS[i], ZZ(1, i).asArray(), 1,
                            ZY1(1, i).asArray(), 1);
                        RES1[i] = dznrm2(M, ZY1(1, i).asArray(), 1);
                      }
                      TMP = ZERO;
                      for (final i in 1.through(K.value)) {
                        TMP = max(
                            TMP,
                            (RES[i] - RES1[i]).abs() *
                                SINGVX[K.value] /
                                (ANORM * SINGVX[1]));
                      }
                      TMP_REZ = max(TMP_REZ, TMP);
                      if (TMP <= TOL2) {
                        // NOUT.println( ':) .... OK ..........zgedmd PASSED.');
                      } else {
                        NFAIL_REZ = NFAIL_REZ + 1;
                        NOUT.println(
                            ':( ..................zgedmd FAILED!Check the code for implementation errors.');
                        NOUT.println(
                            'The input parameters were $SCALE $JOBZ $RESIDS $JOBREF $WHTSVD $M $N $LDX $LDY $NRNK $TOL');
                      }

                      if (lsame(JOBREF, 'E')) {
                        TMP = ZERO;
                        for (final i in 1.through(K.value)) {
                          TMP = max(
                              TMP,
                              (RES1[i] - RESEX[i]).abs() /
                                  (RES1[i] + RESEX[i]));
                        }
                        TMP_EX = max(TMP_EX, TMP);
                      }
                    }

                    // DEALLOCATE(ZWORK)
                    // DEALLOCATE(WORK)
                    // DEALLOCATE(IWORK)

                    if (TEST_QRDMD && (K_TRAJ == 1)) {
                      for (final i in 1.through(M)) {
                        for (final j in 1.through(N + 1)) {
                          ZF[i][j] = ZF0[i][j];
                        }
                      }

                      zgedmdq(
                          SCALE,
                          JOBZ,
                          RESIDS,
                          WANTQ,
                          WANTR,
                          JOBREF,
                          WHTSVD,
                          M,
                          N + 1,
                          ZF,
                          LDF,
                          ZX,
                          LDX,
                          ZY,
                          LDY,
                          NRNK,
                          TOL,
                          K,
                          ZEIGS,
                          ZZ,
                          LDZ,
                          RES,
                          ZAU,
                          LDAU,
                          ZW,
                          LDW,
                          ZS,
                          LDS,
                          ZDUMMY,
                          -1,
                          WDUMMY,
                          -1,
                          IDUMMY,
                          -1,
                          INFO);

                      LZWORK = ZDUMMY[LWMINOPT].toInt();
                      final ZWORK = Array<Complex>(LZWORK);
                      LIWORK = IDUMMY[1];
                      final IWORK = Array<int>(LIWORK);
                      LWORK = WDUMMY[1].toInt();
                      final WORK = Array<double>(LWORK);

                      zgedmdq(
                          SCALE,
                          JOBZ,
                          RESIDS,
                          WANTQ,
                          WANTR,
                          JOBREF,
                          WHTSVD,
                          M,
                          N + 1,
                          ZF,
                          LDF,
                          ZX,
                          LDX,
                          ZY,
                          LDY,
                          NRNK,
                          TOL,
                          KQ,
                          ZEIGS,
                          ZZ,
                          LDZ,
                          RES,
                          ZAU,
                          LDAU,
                          ZW,
                          LDW,
                          ZS,
                          LDS,
                          ZWORK,
                          LZWORK,
                          WORK,
                          LWORK,
                          IWORK,
                          LIWORK,
                          INFO);

                      if (INFO.value != 0) {
                        NOUT.println(
                            'to zgedmdq failed.              &Check the calling sequence and the code.');
                        NOUT.println('The error code is ${INFO.value}');
                        NOUT.println(
                            'The input parameters were $SCALE $JOBZ $RESIDS $WANTQ $WANTR $WHTSVD $M $N $LDX $LDY $NRNK $TOL');
                        return;
                      }

                      for (final i in 1.through(N)) {
                        SINGVQX[i] = WORK[i];
                      }

                      // zgedmdq check point

                      if (1 == 0) {
                        // Comparison of zgedmd and zgedmdq singular values disabled
                        TMP = ZERO;
                        for (final i in 1.through(min(K.value, KQ.value))) {
                          TMP = max(
                              TMP, (SINGVX[i] - SINGVQX[i]).abs() / SINGVX[1]);
                        }
                        SVDIFF = max(SVDIFF, TMP);
                        if (TMP > M * N * EPS) {
                          NOUT.println(
                              'FAILED! Something was wrong with the run.');
                          NFAIL_SVDIFF = NFAIL_SVDIFF + 1;
                          for (final j in 1.through(3)) {
                            NOUT.println('$j ${SINGVX[j]} ${SINGVQX[j]}');
                            //  await NIN.read(); // pause
                          }
                        }
                      }

                      // zgedmdq check point
                      if (lsame(WANTQ, 'Q') && lsame(WANTR, 'R')) {
                        // Check that the QR factors are computed and returned
                        // as requested. The residual ||F-Q*R||_F / ||F||_F
                        // is compared to M*N*EPS.

                        for (final i in 1.through(M)) {
                          for (final j in 1.through(N + 1)) {
                            ZF1[i][j] = ZF0[i][j];
                          }
                        }

                        zgemm('N', 'N', M, N + 1, min(M, N + 1), -Complex.one,
                            ZF, LDF, ZY, LDY, Complex.one, ZF1, LDF);
                        TMP_FQR = zlange('F', M, N + 1, ZF1, LDF, WORK) /
                            zlange('F', M, N + 1, ZF0, LDF, WORK);
                        if (TMP_FQR > TOL2) {
                          NOUT.println(
                              'FAILED! Something was wrong with the run.');
                          NFAIL_F_QR = NFAIL_F_QR + 1;
                        } else {
                          // NOUT.println( '........ PASSED.');
                        }
                      }

                      // zgedmdq check point
                      if (lsame(RESIDS, 'R')) {
                        // Compare the residuals returned by zgedmdq with the
                        // explicitly computed residuals using the matrix A.
                        // Compute explicitly Y1 = A*Z
                        zgemm('N', 'N', M, KQ.value, M, Complex.one, ZA, LDA,
                            ZZ, LDZ, Complex.zero, ZY1, LDY);
                        // and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
                        // of the invariant subspaces that correspond to complex conjugate
                        // pairs of eigencalues. (See the description of Z in zgedmdq)

                        for (final i in 1.through(KQ.value)) {
                          // have a real eigenvalue with real eigenvector
                          zaxpy(M, -ZEIGS[i], ZZ(1, i).asArray(), 1,
                              ZY1(1, i).asArray(), 1);
                          // Y(1:M,i) = Y(1:M,i) - REIG[i]*Z(1:M,i)
                          RES1[i] = dznrm2(M, ZY1(1, i).asArray(), 1);
                        }
                        TMP = ZERO;
                        for (final i in 1.through(KQ.value)) {
                          TMP = max(
                              TMP,
                              (RES[i] - RES1[i]).abs() *
                                  SINGVQX[KQ.value] /
                                  (ANORM * SINGVQX[1]));
                        }
                        TMP_REZQ = max(TMP_REZQ, TMP);
                        if (TMP <= TOL2) {
                          // NOUT.println( '.... OK ........ zgedmdq PASSED.');
                        } else {
                          NFAIL_REZQ = NFAIL_REZQ + 1;
                          NOUT.println(
                              '................ zgedmdq FAILED!Check the code for implementation errors.');
                          return;
                        }
                      }

                      // DEALLOCATE( ZWORK )
                      // DEALLOCATE( WORK  )
                      // DEALLOCATE( IWORK )
                    } // zgedmdq

                    //
                  } // LWMINOPT
                  // NOUT.println( 'LWMINOPT loop completed');
                } // iWHTSVD
                // NOUT.println( 'WHTSVD loop completed');
              } // iNRNK  -2:-1
              // NOUT.println( 'NRNK loop completed');
            } // iSCALE  1:4
            // NOUT.println( 'SCALE loop completed');
          }
          // NOUT.println( 'JOBREF loop completed');
        } // iJOBZ
        // NOUT.println( 'JOBZ loop completed');
      } // MODE -6:6
      // NOUT.println( 'MODE loop completed');
    } // 1 or 2 trajectories
    // NOUT.println( 'trajectories  loop completed');
  } // LLOOP

  NOUT.println('>>>>>>>>>>>>>>>>>>>>>>>>>>');
  NOUT.println(' Test summary for zgedmd :');
  NOUT.println('>>>>>>>>>>>>>>>>>>>>>>>>>>');
  NOUT.println();
  if (NFAIL_Z_XV == 0) {
    NOUT.println('>>>> Z - U*V test PASSED.');
  } else {
    NOUT.println('Z - U*V test FAILED $NFAIL_Z_XV time(s)');
    NOUT.println('Max error ||Z-U*V||_F was $TMP_ZXW');
    NFAIL_TOTAL = NFAIL_TOTAL + NFAIL_Z_XV;
  }
  if (NFAIL_AU == 0) {
    NOUT.println('>>>> A*U test PASSED. ');
  } else {
    NOUT.println('A*U test FAILED $NFAIL_AU time(s)');
    NOUT.println('Max A*U test adjusted error measure was $TMP_AU');
    NOUT.println('It should be up to O(M*N) times EPS, EPS = $EPS');
    NFAIL_TOTAL = NFAIL_TOTAL + NFAIL_AU;
  }

  if (NFAIL_REZ == 0) {
    NOUT.println('>>>> Rezidual computation test PASSED.');
  } else {
    NOUT.println('Rezidual computation test FAILED $NFAIL_REZ time(s)');
    NOUT.println(
        'Max residual computing test adjusted error measure was $TMP_REZ');
    NOUT.println('It should be up to O(M*N) times EPS, EPS = $EPS');
    NFAIL_TOTAL = NFAIL_TOTAL + NFAIL_REZ;
  }

  if (NFAIL_TOTAL == 0) {
    NOUT.println('>>>> zgedmd :: ALL TESTS PASSED.');
  } else {
    NOUT.println('$NFAIL_TOTAL FAILURES!');
    NOUT.println(
        '>>>>>>>>>>>>>> zgedmd :: TESTS FAILED. CHECK THE IMPLEMENTATION.');
  }

  if (TEST_QRDMD) {
    NOUT.println();
    NOUT.println('>>>>>>>>>>>>>>>>>>>>>>>>>>');
    NOUT.println(' Test summary for zgedmdq :');
    NOUT.println('>>>>>>>>>>>>>>>>>>>>>>>>>>');
    NOUT.println();

    if (NFAIL_SVDIFF == 0) {
      NOUT.println(
          '>>>> zgedmd and zgedmdq computed singular               &values test PASSED.');
    } else {
      NOUT.println(
          'zgedmd and zgedmdq discrepancies in              &the singular values unacceptable $NFAIL_SVDIFF times. Test FAILED.');
      NOUT.println(
          'The maximal discrepancy in the singular values (relative to the norm) was $SVDIFF');
      NOUT.println('It should be up to O(M*N) times EPS, EPS = $EPS');
      NFAILQ_TOTAL = NFAILQ_TOTAL + NFAIL_SVDIFF;
    }

    if (NFAIL_F_QR == 0) {
      NOUT.println('>>>> F - Q*R test PASSED.');
    } else {
      NOUT.println('F - Q*R test FAILED $NFAIL_F_QR time(s)');
      NOUT.println('The largest relative residual was $TMP_FQR');
      NOUT.println('It should be up to O(M*N) times EPS, EPS = $EPS');
      NFAILQ_TOTAL = NFAILQ_TOTAL + NFAIL_F_QR;
    }

    if (NFAIL_REZQ == 0) {
      NOUT.println('>>>> Rezidual computation test PASSED.');
    } else {
      NOUT.println('Rezidual computation test FAILED $NFAIL_REZQ time(s)');
      NOUT.println(
          'Max residual computing test adjusted error measure was $TMP_REZQ');
      NOUT.println('It should be up to O(M*N) times EPS, EPS = $EPS');
      NFAILQ_TOTAL = NFAILQ_TOTAL + NFAIL_REZQ;
    }

    if (NFAILQ_TOTAL == 0) {
      NOUT.println('>>>>>>> zgedmdq :: ALL TESTS PASSED.');
    } else {
      NOUT.println('$NFAILQ_TOTAL FAILURES!');
      NOUT.println(
          '>>>>>>> zgedmdq :: TESTS FAILED. CHECK THE IMPLEMENTATION.');
    }
  }

  NOUT.println();
  NOUT.println('Test completed.');
}
