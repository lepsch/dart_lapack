import 'dart:io';
import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgedmd.dart';
import 'package:lapack/src/dgedmdq.dart';
import 'package:lapack/src/dgeev.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:test/test.dart';

import '../matgen/dlatmr.dart';
import '../test_driver.dart';

Future<void> dchkdmd(
  final Nin NIN,
  final Nout NOUT,
  final TestDriver test,
) async {
  const ONE = 1.0, ZERO = 0.0;
  final AB = Matrix<double>(2, 2), WDUMMY = Array<double>(2);
  final IDUMMY = Array<int>(2), RJOBDATA = Array<int>(8);
  double SVDIFF = 0,
      TMP_AU = 0,
      TMP_FQR = 0,
      TMP_REZ = 0,
      TMP_REZQ = 0,
      TMP_ZXW = 0,
      TMP_EX;
  final INFO = Box(0);

  // The test is always in pairs : ( DGEDMD and DGEDMDQ )
  // because the test includes comparing the results (in pairs).

  const TEST_QRDMD = true; // This code by default performs tests on DGEDMDQ
  // Since the QR factorizations based algorithm is designed for
  // single trajectory data, only single trajectory tests will
  // be performed with xGEDMDQ.
  const WANTQ = 'Q';
  const WANTR = 'R';

  final EPS = dlamch('P'); // machine precision DP

  // Global counters of failures of some particular tests
  var NFAIL_REZ = 0;
  var NFAIL_REZQ = 0;
  var NFAIL_Z_XV = 0;
  var NFAIL_F_QR = 0;
  var NFAIL_AU = 0;
  var KDIFF = 0;
  var NFAIL_SVDIFF = 0;
  var NFAIL_TOTAL = 0;
  var NFAILQ_TOTAL = 0;

  for (final LLOOP in 1.through(4)) {
    NOUT.println('L Loop Index = $LLOOP');

    // Set the dimensions of the problem ...
    NOUT.println('M = ');
    final M = await NIN.readInt();
    NOUT.println(M.toString());
    // ... and the number of snapshots.
    NOUT.println('N = ');
    final N = await NIN.readInt();
    NOUT.println(N.toString());

    // ... Test the dimensions
    if ((min(M, N) == 0) || (M < N)) {
      NOUT.println('Bad dimensions. Required: M >= N > 0.');
      return;
    }

    // The seed inside the LLOOP so that each pass can be reproduced easily.
    final ISEED = Array.fromList([4, 3, 2, 1]);

    final LDA = M;
    final LDF = M;
    final LDX = max(M, N + 1);
    final LDY = max(M, N + 1);
    final LDW = N;
    final LDZ = M;
    final LDAU = max(M, N + 1);
    final LDS = N;

    TMP_ZXW = ZERO;
    TMP_AU = ZERO;
    TMP_REZ = ZERO;
    TMP_REZQ = ZERO;
    SVDIFF = ZERO;
    TMP_EX = ZERO;

    // Test the subroutines on real data snapshots. All
    // computation is done in real arithmetic, even when
    // Koopman eigenvalues and modes are real.

    final A = Matrix<double>(LDA, M),
        AC = Matrix<double>(LDA, M),
        EIGA = Matrix<double>(M, 2),
        LAMBDA = Matrix<double>(N, 2),
        LAMBDAQ = Matrix<double>(N, 2),
        F = Matrix<double>(LDF, N + 1),
        F1 = Matrix<double>(LDF, N + 1),
        F2 = Matrix<double>(LDF, N + 1),
        Z = Matrix<double>(LDZ, N),
        Z1 = Matrix<double>(LDZ, N),
        S = Matrix<double>(N, N),
        AU = Matrix<double>(LDAU, N),
        W = Matrix<double>(LDW, N),
        VA = Matrix<double>(LDA, M),
        X = Matrix<double>(LDX, N),
        X0 = Matrix<double>(LDX, N),
        Y = Matrix<double>(LDY, N + 1),
        Y0 = Matrix<double>(LDY, N + 1),
        Y1 = Matrix<double>(M, N + 1);
    final DA = Array<double>(M),
        DL = Array<double>(M),
        REIG = Array<double>(N),
        REIGA = Array<double>(M),
        REIGQ = Array<double>(N),
        IEIG = Array<double>(N),
        IEIGA = Array<double>(M),
        IEIGQ = Array<double>(N),
        RES = Array<double>(N),
        RES1 = Array<double>(N),
        RESEX = Array<double>(N),
        SINGVX = Array<double>(N),
        SINGVQX = Array<double>(N);

    final TOL = M * EPS;
    // This mimics O(M*N)*EPS bound for accumulated roundoff error.
    // The factor 10 is somewhat arbitrary.
    final TOL2 = 10 * M * N * EPS;

    for (final K_TRAJ in [1, 2]) {
      //  Number of intial conditions in the simulation/trajectories (1 or 2)

      const COND = 1.0e8;
      const DMAX = 1.0e2;
      const RSIGN = 'F';
      const GRADE = 'N';
      const MODEL = 6;
      const CONDL = 1.0e2;
      const MODER = 6;
      const CONDR = 1.0e2;
      const PIVTNG = 'N';

      // Loop over all parameter MODE values for ZLATMR (+1,..,+6)
      for (final MODE in [1, 2, 3, 4, 5, 6]) {
        final IWORK = Array<int>(2 * M);
        final DR = Array<double>(N);
        dlatmr(
            M,
            M,
            'S',
            ISEED,
            'N',
            DA,
            MODE,
            COND,
            DMAX,
            RSIGN,
            GRADE,
            DL,
            MODEL,
            CONDL,
            DR,
            MODER,
            CONDR,
            PIVTNG,
            IWORK,
            M,
            M,
            ZERO,
            -ONE,
            'N',
            A,
            LDA,
            IWORK(M + 1),
            INFO);

        final LWORK = 4 * M + 1;
        final WORK = Array<double>(LWORK);
        // AC = A;
        for (final i in 1.through(M)) {
          for (final j in 1.through(M)) {
            AC[i][j] = A[i][j];
          }
        }

        dgeev(
            'N', 'V', M, AC, M, REIGA, IEIGA, VA, M, VA, M, WORK, LWORK, INFO);
        var TMP = ZERO;
        for (final i in 1.through(M)) {
          EIGA[i][1] = REIGA[i];
          EIGA[i][2] = IEIGA[i];
          TMP = max(TMP, sqrt(pow(REIGA[i], 2) + pow(IEIGA[i], 2)));
        }

        // Scale A to have the desirable spectral radius.
        dlascl('G', 0, 0, TMP, ONE, M, M, A, M, INFO);
        dlascl('G', 0, 0, TMP, ONE, M, 2, EIGA, M, INFO);

        // Compute the norm of A
        final ANORM = dlange('F', N, N, A, M, WDUMMY);

        if (K_TRAJ == 2) {
          // generate data with two inital conditions
          dlarnv(2, ISEED, M, F1(1, 1).asArray());
          // F1(1:M,1) = 1.0E-10*F1(1:M,1);
          for (final i in 1.through(M)) {
            F1[i][1] = 1.0E-10 * F1[i][1];
          }
          for (final i in 1.through(N ~/ 2)) {
            dgemv('N', M, M, ONE, A, M, F1(1, i).asArray(), 1, ZERO,
                F1(1, i + 1).asArray(), 1);
          }
          // X0(1:M,1:N/2) = F1(1:M,1:N/2);
          // Y0(1:M,1:N/2) = F1(1:M,2:N/2+1);
          for (final i in 1.through(M)) {
            for (final j in 1.through(N ~/ 2)) {
              X0[i][j] = F1[i][j];
              Y0[i][j] = F1[i][j + 1];
            }
          }

          dlarnv(2, ISEED, M, F1(1, 1).asArray());
          for (final i in 1.through(N - N ~/ 2)) {
            dgemv('N', M, M, ONE, A, M, F1(1, i).asArray(), 1, ZERO,
                F1(1, i + 1).asArray(), 1);
          }
          // X0(1:M,N/2+1:N) = F1(1:M,1:N-N/2)
          // Y0(1:M,N/2+1:N) = F1(1:M,2:N-N/2+1)
          for (final i in 1.through(M)) {
            for (final j in 1.through(N - N ~/ 2)) {
              X0[i][j + N ~/ 2] = F1[i][j];
              Y0[i][j + N ~/ 2] = F1[i][j + 1];
            }
          }
        } else {
          dlarnv(2, ISEED, M, F(1, 1).asArray());
          for (final i in 1.through(N)) {
            dgemv('N', M, M, ONE, A, M, F(1, i).asArray(), 1, ZERO,
                F(1, i + 1).asArray(), 1);
          }
          // X0(1:M,1:N) = F(1:M,1:N)
          // Y0(1:M,1:N) = F(1:M,2:N+1)
          for (final i in 1.through(M)) {
            for (final j in 1.through(N)) {
              X0[i][j] = F[i][j];
              Y0[i][j] = F[i][j + 1];
            }
          }
        }

        // dlange('F', M, N, X0, LDX, WDUMMY);
        // dlange('F', M, N, Y0, LDX, WDUMMY);

        final ctx = (
          A: A.copy(),
          F: F.copy(),
          F1: F1.copy(),
          X0: X0.copy(),
          Y0: Y0.copy(),
        );

        for (final (JOBZ, RESIDS) in [
          (
            'V', // Ritz vectors will be computed
            'R', // Residuals will be computed
          ),
          ('V', 'N'),
          (
            'F', // Ritz vectors in factored form
            'N'
          ),
          ('N', 'N'),
        ]) {
          final (:A, :F, :F1, :X0, :Y0) = ctx;

          for (final JOBREF in [
            'R', // Data for refined Ritz vectors
            'E', // Exact DMD vectors
            'N',
          ]) {
            for (final SCALE in [
              'S', // X data normalized
              'C', // X normalized, consist. check
              'Y', // Y data normalized
              'N',
            ]) {
              for (final NRNK in [-1, -2]) {
                // Two truncation strategies. The "-2" case for R&D
                // purposes only - it uses possibly low accuracy small
                // singular values, in which case the formulas used in
                // the DMD are highly sensitive.

                for (final WHTSVD in [1, 2, 3, 4]) {
                  // Check all four options to compute the POD basis
                  // via the SVD.

                  for (final LWMINOPT in [1, 2]) {
                    // Workspace query for the minimal (1) and for the optimal
                    // (2) workspace lengths determined by workspace query.

                    test(
                        'DMD: Dynamic Mode Decomposition (LLOOP=$LLOOP TRAJ=$K_TRAJ MODE=$MODE JOBZ=$JOBZ RESIDS=$RESIDS JOBREF=$JOBREF SCALE=$SCALE NRNK=$NRNK WHTSVD=$WHTSVD LWMINOPT=$LWMINOPT)',
                        () {
                      // X(1:M,1:N) = X0(1:M,1:N)
                      // Y(1:M,1:N) = Y0(1:M,1:N)
                      for (final i in 1.through(M)) {
                        for (final j in 1.through(N)) {
                          X[i][j] = X0[i][j];
                          Y[i][j] = Y0[i][j];
                        }
                      }

                      // DGEDMD: Workspace query and workspace allocation
                      final K = Box(0);
                      dgedmd(
                          SCALE,
                          JOBZ,
                          RESIDS,
                          JOBREF,
                          WHTSVD,
                          M,
                          N,
                          X,
                          LDX,
                          Y,
                          LDY,
                          NRNK,
                          TOL,
                          K,
                          REIG,
                          IEIG,
                          Z,
                          LDZ,
                          RES,
                          AU,
                          LDAU,
                          W,
                          LDW,
                          S,
                          LDS,
                          WDUMMY,
                          -1,
                          IDUMMY,
                          -1,
                          INFO);

                      final LIWORK = IDUMMY[1];
                      final IWORK = Array<int>(LIWORK);
                      final LWORK = WDUMMY[LWMINOPT].toInt();
                      final WORK = Array<double>(LWORK);

                      // DGEDMD test: dgedmd
                      dgedmd(
                          SCALE,
                          JOBZ,
                          RESIDS,
                          JOBREF,
                          WHTSVD,
                          M,
                          N,
                          X,
                          LDX,
                          Y,
                          LDY,
                          NRNK,
                          TOL,
                          K,
                          REIG,
                          IEIG,
                          Z,
                          LDZ,
                          RES,
                          AU,
                          LDAU,
                          W,
                          LDW,
                          S,
                          LDS,
                          WORK,
                          LWORK,
                          IWORK,
                          LIWORK,
                          INFO);

                      // SINGVX(1:N) = WORK(1:N);
                      for (final i in 1.through(N)) {
                        SINGVX[i] = WORK[i];
                      }

                      //. DGEDMD check point
                      if (lsame(JOBZ, 'V')) {
                        // Check that Z = X*W, on return from DGEDMD
                        // This checks that the returned aigenvectors in Z are
                        // the product of the SVD'POD basis returned in X
                        // and the eigenvectors of the rayleigh quotient
                        // returned in W
                        dgemm('N', 'N', M, K.value, K.value, ONE, X, LDX, W,
                            LDW, ZERO, Z1, LDZ);
                        var TMP = ZERO;
                        for (final i in 1.through(K.value)) {
                          daxpy(M, -ONE, Z(1, i).asArray(), 1,
                              Z1(1, i).asArray(), 1);
                          TMP = max(TMP, dnrm2(M, Z1(1, i).asArray(), 1));
                        }
                        TMP_ZXW = max(TMP_ZXW, TMP);

                        test.expect(TMP_ZXW, lessThanOrEqualTo(10 * M * EPS));
                        if (TMP_ZXW > 10 * M * EPS) {
                          NFAIL_Z_XV = NFAIL_Z_XV + 1;
                          NOUT.println(
                              ':( .................DGEDMD FAILED!Check the code for implementation errors.');
                          NOUT.println(
                              'The input parameters were $SCALE $JOBZ $RESIDS $JOBREF ${WHTSVD.i12}${M.i12}${N.i12}${LDX.i12}${LDY.i12}${NRNK.i12}${TOL.d12_4}');
                        }
                      }

                      //. DGEDMD check point
                      if (lsame(JOBREF, 'R')) {
                        // The matrix A*U is returned for computing refined Ritz vectors.
                        // Check that A*U is computed correctly using the formula
                        // A*U = Y * V * inv(SIGMA). This depends on the
                        // accuracy in the computed singular values and vectors of X.
                        // See the paper for an error analysis.
                        // Note that the left singular vectors of the input matrix X
                        // are returned in the array X.
                        dgemm('N', 'N', M, K.value, M, ONE, A, LDA, X, LDX,
                            ZERO, Z1, LDZ);
                        var TMP = ZERO;
                        for (final i in 1.through(K.value)) {
                          daxpy(M, -ONE, AU(1, i).asArray(), 1,
                              Z1(1, i).asArray(), 1);
                          TMP = max(
                              TMP,
                              dnrm2(M, Z1(1, i).asArray(), 1) *
                                  SINGVX[K.value] /
                                  (ANORM * SINGVX[1]));
                        }
                        TMP_AU = max(TMP_AU, TMP);

                        test.expect(TMP, lessThanOrEqualTo(TOL2));
                        if (TMP > TOL2) {
                          NFAIL_AU = NFAIL_AU + 1;
                          NOUT.println(
                              ':( .................DGEDMD FAILED!Check the code for implementation errors.');
                          NOUT.println(
                              'The input parameters were $SCALE $JOBZ $RESIDS $JOBREF ${WHTSVD.i12}${M.i12}${N.i12}${LDX.i12}${LDY.i12}${NRNK.i12}${TOL.d12_4}');
                        }
                      } else if (lsame(JOBREF, 'E')) {
                        // The unscaled vectors of the Exact DMD are computed.
                        // This option is included for the sake of completeness,
                        // for users who prefer the Exact DMD vectors. The
                        // returned vectors are in the real form, in the same way
                        // as the Ritz vectors. Here we just save the vectors
                        // and test them separately using a Matlab script.

                        dgemm('N', 'N', M, K.value, M, ONE, A, LDA, AU, LDAU,
                            ZERO, Y1, M);
                        var i = 1;
                        while (i <= K.value) {
                          if (IEIG[i] == ZERO) {
                            // have a real eigenvalue with real eigenvector
                            daxpy(M, -REIG[i], AU(1, i).asArray(), 1,
                                Y1(1, i).asArray(), 1);
                            RESEX[i] = dnrm2(M, Y1(1, i).asArray(), 1) /
                                dnrm2(M, AU(1, i).asArray(), 1);
                            i = i + 1;
                          } else {
                            // Have a complex conjugate pair
                            // REIG(i) +- sqrt(-1)*IMEIG(i).
                            // Since all computation is done in real
                            // arithmetic, the formula for the residual
                            // is recast for real representation of the
                            // complex conjugate eigenpair. See the
                            // description of RES.
                            AB[1][1] = REIG[i];
                            AB[2][1] = -IEIG[i];
                            AB[1][2] = IEIG[i];
                            AB[2][2] = REIG[i];
                            dgemm('N', 'N', M, 2, 2, -ONE, AU(1, i), M, AB, 2,
                                ONE, Y1(1, i), M);
                            RESEX[i] = dlange('F', M, 2, Y1(1, i), M, WORK) /
                                dlange('F', M, 2, AU(1, i), M, WORK);
                            RESEX[i + 1] = RESEX[i];
                            i = i + 2;
                          }
                        }
                      }

                      //. DGEDMD check point
                      if (lsame(RESIDS, 'R')) {
                        // Compare the residuals returned by DGEDMD with the
                        // explicitly computed residuals using the matrix A.
                        // Compute explicitly Y1 = A*Z
                        dgemm('N', 'N', M, K.value, M, ONE, A, LDA, Z, LDZ,
                            ZERO, Y1, M);
                        // ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
                        // of the invariant subspaces that correspond to complex conjugate
                        // pairs of eigencalues. (See the description of Z in DGEDMD,)
                        var i = 1;
                        while (i <= K.value) {
                          if (IEIG[i] == ZERO) {
                            // have a real eigenvalue with real eigenvector
                            daxpy(M, -REIG[i], Z(1, i).asArray(), 1,
                                Y1(1, i).asArray(), 1);
                            RES1[i] = dnrm2(M, Y1(1, i).asArray(), 1);
                            i = i + 1;
                          } else {
                            // Have a complex conjugate pair
                            // REIG(i) +- sqrt(-1)*IMEIG(i).
                            // Since all computation is done in real
                            // arithmetic, the formula for the residual
                            // is recast for real representation of the
                            // complex conjugate eigenpair. See the
                            // description of RES.
                            AB[1][1] = REIG[i];
                            AB[2][1] = -IEIG[i];
                            AB[1][2] = IEIG[i];
                            AB[2][2] = REIG[i];
                            dgemm('N', 'N', M, 2, 2, -ONE, Z(1, i), M, AB, 2,
                                ONE, Y1(1, i), M);
                            RES1[i] = dlange('F', M, 2, Y1(1, i), M, WORK);
                            RES1[i + 1] = RES1[i];
                            i = i + 2;
                          }
                        }
                        var TMP = ZERO;
                        for (final i in 1.through(K.value)) {
                          TMP = max(
                              TMP,
                              (RES[i] - RES1[i]).abs() *
                                  SINGVX[K.value] /
                                  (ANORM * SINGVX[1]));
                        }
                        TMP_REZ = max(TMP_REZ, TMP);

                        test.expect(TMP, lessThanOrEqualTo(TOL2));
                        if (TMP > TOL2) {
                          NFAIL_REZ = NFAIL_REZ + 1;
                          NOUT.println(
                              ':( ..................DGEDMD FAILED!Check the code for implementation errors.');
                          NOUT.println(
                              'The input parameters were $SCALE $JOBZ $RESIDS $JOBREF ${WHTSVD.i12}${M.i12}${N.i12}${LDX.i12}${LDY.i12}${NRNK.i12}${TOL.d12_4}');
                        }

                        if (lsame(JOBREF, 'E')) {
                          var TMP = ZERO;
                          for (final i in 1.through(K.value)) {
                            TMP = max(
                                TMP,
                                (RES1[i] - RESEX[i]).abs() /
                                    (RES1[i] + RESEX[i]));
                          }
                          TMP_EX = max(TMP_EX, TMP);
                        }
                      }

                      // store the results for inspection
                      for (final i in 1.through(K.value)) {
                        LAMBDA[i][1] = REIG[i];
                        LAMBDA[i][2] = IEIG[i];
                      }

                      //======================================================================
                      //     Now test the DGEDMDQ
                      //======================================================================
                      if (TEST_QRDMD && (K_TRAJ == 1)) {
                        RJOBDATA[2] = 1;
                        // F1 = F;
                        for (final i in 1.through(M)) {
                          for (final j in 1.through(N + 1)) {
                            F1[i][j] = F[i][j];
                          }
                        }

                        // DGEDMDQ test: Workspace query and workspace allocation
                        final KQ = Box(0);
                        dgedmdq(
                            SCALE,
                            JOBZ,
                            RESIDS,
                            WANTQ,
                            WANTR,
                            JOBREF,
                            WHTSVD,
                            M,
                            N + 1,
                            F1,
                            LDF,
                            X,
                            LDX,
                            Y,
                            LDY,
                            NRNK,
                            TOL,
                            KQ,
                            REIGQ,
                            IEIGQ,
                            Z,
                            LDZ,
                            RES,
                            AU,
                            LDAU,
                            W,
                            LDW,
                            S,
                            LDS,
                            WDUMMY,
                            -1,
                            IDUMMY,
                            -1,
                            INFO);
                        final LIWORK = IDUMMY[1];
                        final IWORK = Array<int>(LIWORK);
                        final LWORK = WDUMMY[LWMINOPT].toInt();
                        final WORK = Array<double>(LWORK);
                        // DGEDMDQ test: dgedmdq
                        dgedmdq(
                            SCALE,
                            JOBZ,
                            RESIDS,
                            WANTQ,
                            WANTR,
                            JOBREF,
                            WHTSVD,
                            M,
                            N + 1,
                            F1,
                            LDF,
                            X,
                            LDX,
                            Y,
                            LDY,
                            NRNK,
                            TOL,
                            KQ,
                            REIGQ,
                            IEIGQ,
                            Z,
                            LDZ,
                            RES,
                            AU,
                            LDAU,
                            W,
                            LDW,
                            S,
                            LDS,
                            WORK,
                            LWORK,
                            IWORK,
                            LIWORK,
                            INFO);

                        // SINGVQX(1:KQ) = WORK(min(M,N+1)+1: min(M,N+1)+KQ);
                        for (final i in 1.through(KQ.value)) {
                          SINGVQX[i] = WORK[i + min(M, N + 1)];
                        }

                        // DGEDMDQ check point
                        if (KQ.value != K.value) {
                          KDIFF = KDIFF + 1;
                        }

                        var TMP = ZERO;
                        for (final i in 1.through(min(K.value, KQ.value))) {
                          TMP = max(
                              TMP, (SINGVX[i] - SINGVQX[i]).abs() / SINGVX[1]);
                        }
                        SVDIFF = max(SVDIFF, TMP);
                        test.expect(TMP, lessThanOrEqualTo(M * N * EPS));
                        if (TMP > M * N * EPS) {
                          NOUT.println(
                              'FAILED! Something was wrong with the run.');
                          NFAIL_SVDIFF = NFAIL_SVDIFF + 1;
                          for (final j in 1.through(3)) {
                            NOUT.println(
                                '${j.i12}${SINGVX[j].d12_4}${SINGVQX[j].d12_4}');
                            // pause(); // await NIN.read();
                          }
                        }

                        // DGEDMDQ check point
                        if (lsame(WANTQ, 'Q') && lsame(WANTR, 'R')) {
                          // Check that the QR factors are computed and returned
                          // as requested. The residual ||F-Q*R||_F / ||F||_F
                          // is compared to M*N*EPS.

                          // F2 = F;
                          for (final i in 1.through(M)) {
                            for (final j in 1.through(N + 1)) {
                              F2[i][j] = F[i][j];
                            }
                          }

                          dgemm('N', 'N', M, N + 1, min(M, N + 1), -ONE, F1,
                              LDF, Y, LDY, ONE, F2, LDF);
                          TMP_FQR = dlange('F', M, N + 1, F2, LDF, WORK) /
                              dlange('F', M, N + 1, F, LDF, WORK);
                          test.expect(TMP_FQR, lessThanOrEqualTo(TOL2));
                          if (TMP_FQR > TOL2) {
                            NOUT.println(
                                'FAILED! Something was wrong with the run.');
                            NFAIL_F_QR = NFAIL_F_QR + 1;
                          }
                        }

                        // DGEDMDQ check point
                        if (lsame(RESIDS, 'R')) {
                          // Compare the residuals returned by DGEDMDQ with the
                          // explicitly computed residuals using the matrix A.
                          // Compute explicitly Y1 = A*Z
                          dgemm('N', 'N', M, KQ.value, M, ONE, A, M, Z, M, ZERO,
                              Y1, M);
                          // ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
                          // of the invariant subspaces that correspond to complex conjugate
                          // pairs of eigencalues. (See the description of Z in DGEDMDQ)
                          var i = 1;
                          while (i <= KQ.value) {
                            if (IEIGQ[i] == ZERO) {
                              // have a real eigenvalue with real eigenvector
                              daxpy(M, -REIGQ[i], Z(1, i).asArray(), 1,
                                  Y1(1, i).asArray(), 1);
                              // Y(1:M,i) = Y(1:M,i) - REIG(i)*Z(1:M,i)
                              RES1[i] = dnrm2(M, Y1(1, i).asArray(), 1);
                              i = i + 1;
                            } else {
                              // Have a complex conjugate pair
                              // REIG(i) +- sqrt(-1)*IMEIG(i).
                              // Since all computation is done in real
                              // arithmetic, the formula for the residual
                              // is recast for real representation of the
                              // complex conjugate eigenpair. See the
                              // description of RES.
                              AB[1][1] = REIGQ[i];
                              AB[2][1] = -IEIGQ[i];
                              AB[1][2] = IEIGQ[i];
                              AB[2][2] = REIGQ[i];
                              dgemm('N', 'N', M, 2, 2, -ONE, Z(1, i), M, AB, 2,
                                  ONE, Y1(1, i), M);
                              // Y(1:M,i:i+1) = Y(1:M,i:i+1) - Z(1:M,i:i+1) * AB   // INTRINSIC
                              RES1[i] = dlange('F', M, 2, Y1(1, i), M, WORK);
                              RES1[i + 1] = RES1[i];
                              i = i + 2;
                            }
                          }
                          var TMP = ZERO;
                          for (final i in 1.through(KQ.value)) {
                            TMP = max(
                                TMP,
                                (RES[i] - RES1[i]).abs() *
                                    SINGVQX[K.value] /
                                    (ANORM * SINGVQX[1]));
                          }
                          TMP_REZQ = max(TMP_REZQ, TMP);
                          test.expect(TMP, lessThanOrEqualTo(TOL2));
                          if (TMP > TOL2) {
                            NFAIL_REZQ = NFAIL_REZQ + 1;
                            NOUT.println(
                                '................ DGEDMDQ FAILED!Check the code for implementation errors.');
                            return;
                          }
                        }

                        for (final i in 1.through(KQ.value)) {
                          LAMBDAQ[i][1] = REIGQ[i];
                          LAMBDAQ[i][2] = IEIGQ[i];
                        }
                      } // TEST_QRDMD
                    });
                  } // LWMINOPT
                  //NOUT.println( 'LWMINOPT loop completed');
                } // WHTSVD LOOP
                //NOUT.println( 'WHTSVD loop completed');
              } // NRNK LOOP
              //NOUT.println( 'NRNK loop completed');
            } // SCALE LOOP
            //NOUT.println( 'SCALE loop completed');
          } // JOBF LOOP
          //NOUT.println( 'JOBREF loop completed');
        } // JOBZ LOOP
        //NOUT.println( 'JOBZ loop completed');
      } // MODE -6:6
      //NOUT.println( 'MODE loop completed');
    } // 1 or 2 trajectories
    //NOUT.println( 'trajectories  loop completed');

    //     Generate random M-by-M matrix A. Use DLATMR from
  } // LLOOP

  NOUT.println('>>>>>>>>>>>>>>>>>>>>>>>>>>');
  NOUT.println(' Test summary for DGEDMD :');
  NOUT.println('>>>>>>>>>>>>>>>>>>>>>>>>>>');
  NOUT.println();
  if (NFAIL_Z_XV == 0) {
    NOUT.println('>>>> Z - U*V test PASSED.');
  } else {
    NOUT.println('Z - U*V test FAILED ${NFAIL_Z_XV.i12} time(s)');
    NOUT.println('Max error ||Z-U*V||_F was ${TMP_ZXW.d12_4}');
    NFAIL_TOTAL = NFAIL_TOTAL + NFAIL_Z_XV;
  }
  if (NFAIL_AU == 0) {
    NOUT.println('>>>> A*U test PASSED. ');
  } else {
    NOUT.println('A*U test FAILED ${NFAIL_AU.i12} time(s)');
    NOUT.println('Max A*U test adjusted error measure was ${TMP_AU.d12_4}');
    NOUT.println('It should be up to O(M*N) times EPS, EPS = ${EPS.d12_4}');
    NFAIL_TOTAL = NFAIL_TOTAL + NFAIL_AU;
  }

  if (NFAIL_REZ == 0) {
    NOUT.println('>>>> Rezidual computation test PASSED.');
  } else {
    NOUT.println('Rezidual computation test FAILED ${NFAIL_REZ.i12}time(s)');
    NOUT.println(
        'Max residual computing test adjusted error measure was ${TMP_REZ.d12_4}');
    NOUT.println('It should be up to O(M*N) times EPS, EPS = ${EPS.d12_4}');
    NFAIL_TOTAL = NFAIL_TOTAL + NFAIL_REZ;
  }

  if (NFAIL_TOTAL == 0) {
    NOUT.println('>>>> DGEDMD :: ALL TESTS PASSED.');
  } else {
    NOUT.println('${NFAIL_TOTAL.i12} FAILURES!');
    NOUT.println(
        '>>>>>>>>>>>>>> DGEDMD :: TESTS FAILED. CHECK THE IMPLEMENTATION.');
  }

  if (TEST_QRDMD) {
    NOUT.println();
    NOUT.println('>>>>>>>>>>>>>>>>>>>>>>>>>>');
    NOUT.println(' Test summary for DGEDMDQ :');
    NOUT.println('>>>>>>>>>>>>>>>>>>>>>>>>>>');
    NOUT.println();

    if (NFAIL_SVDIFF == 0) {
      NOUT.println(
          '>>>> DGEDMD and DGEDMDQ computed singular               &values test PASSED.');
    } else {
      NOUT.println(
          'DGEDMD and DGEDMDQ discrepancies in               &the singular values unacceptable ${NFAIL_SVDIFF.i12} times. Test FAILED.');
      NOUT.println(
          'The maximal discrepancy in the singular values (relative to the norm) was ${SVDIFF.d12_4}');
      NOUT.println('It should be up to O(M*N) times EPS, EPS = ${EPS.d12_4}');
      NFAILQ_TOTAL = NFAILQ_TOTAL + NFAIL_SVDIFF;
    }

    if (NFAIL_F_QR == 0) {
      NOUT.println('>>>> F - Q*R test PASSED.');
    } else {
      NOUT.println('F - Q*R test FAILED ${NFAIL_F_QR.i12} time(s)');
      NOUT.println('The largest relative residual was ${TMP_FQR.d12_4}');
      NOUT.println('It should be up to O(M*N) times EPS, EPS = ${EPS.d12_4}');
      NFAILQ_TOTAL = NFAILQ_TOTAL + NFAIL_F_QR;
    }

    if (NFAIL_REZQ == 0) {
      NOUT.println('>>>> Rezidual computation test PASSED.');
    } else {
      NOUT.println('Rezidual computation test FAILED ${NFAIL_REZQ.i12}time(s)');
      NOUT.println(
          'Max residual computing test adjusted error measure was ${TMP_REZQ.d12_4}');
      NOUT.println('It should be up to O(M*N) times EPS, EPS = ${EPS.d12_4}');
      NFAILQ_TOTAL = NFAILQ_TOTAL + NFAIL_REZQ;
    }

    if (NFAILQ_TOTAL == 0) {
      NOUT.println('>>>>>>> DGEDMDQ :: ALL TESTS PASSED.');
    } else {
      NOUT.println('${NFAILQ_TOTAL.i12}FAILURES!');
      NOUT.println(
          '>>>>>>> DGEDMDQ :: TESTS FAILED. CHECK THE IMPLEMENTATION.');
    }
  }

  NOUT.println();
  NOUT.println('Test completed.');
}

void main() async {
  await dchkdmd(Nin(stdin), Nout(stdout), syncTestDriver);
  exit(syncTestDriver.errors);
}
