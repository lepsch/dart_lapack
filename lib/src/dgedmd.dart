// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgeev.dart';
import 'package:dart_lapack/src/dgejsv.dart';
import 'package:dart_lapack/src/dgesdd.dart';
import 'package:dart_lapack/src/dgesvd.dart';
import 'package:dart_lapack/src/dgesvdq.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/dlassq.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/range.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgedmd(
  final String JOBS,
  final String JOBZ,
  final String JOBR,
  final String JOBF,
  final int WHTSVD,
  final int M,
  final int N,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> Y_,
  final int LDY,
  final int NRNK,
  final double TOL,
  final Box<int> K,
  final Array<double> REIG_,
  final Array<double> IMEIG_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> RES_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> W_,
  final int LDW,
  final Matrix<double> S_,
  final int LDS,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
//  -- LAPACK driver routine                                           --
//
//  -- LAPACK is a software package provided by University of          --
//  -- Tennessee, University of California Berkeley, University of     --
//  -- Colorado Denver and NAG Ltd..                                   --

  final X = X_.having(ld: LDX), Y = Y_.having(ld: LDY);
  final Z = Z_.having(ld: LDZ),
      B = B_.having(ld: LDB),
      W = W_.having(ld: LDW),
      S = S_.having(ld: LDS);
  final REIG = REIG_.having(), IMEIG = IMEIG_.having(), RES = RES_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  double OFL, ROOTSC, SMALL, XSCL1 = 0, XSCL2 = 0;
  int i,
      IMINWR = 0,
      LWRKEV,
      LWRSDD,
      LWRSVD,
      LWRSVQ,
      MLWORK = 0,
      MWRKEV,
      MWRSDD,
      MWRSVD,
      MWRSVJ,
      MWRSVQ,
      OLWORK = 0;
  bool BADXY, LQUERY, SCCOLX, SCCOLY, WNTEX, WNTREF, WNTRES, WNTVEC;
  String JOBZL = '', T_OR_N = '';
  String JSVOPT = '';
  final AB = Matrix<double>(2, 2),
      RDUMMY = Array<double>(2),
      RDUMMY2 = Array<double>(2);
  final INFO1 = Box(0), INFO2 = Box(0), NUMRNK = Box(0);

  // Test the input arguments

  WNTRES = lsame(JOBR, 'R');
  SCCOLX = lsame(JOBS, 'S') || lsame(JOBS, 'C');
  SCCOLY = lsame(JOBS, 'Y');
  WNTVEC = lsame(JOBZ, 'V');
  WNTREF = lsame(JOBF, 'R');
  WNTEX = lsame(JOBF, 'E');
  INFO.value = 0;
  LQUERY = ((LWORK == -1) || (LIWORK == -1));

  if (!(SCCOLX || SCCOLY || lsame(JOBS, 'N'))) {
    INFO.value = -1;
  } else if (!(WNTVEC || lsame(JOBZ, 'N') || lsame(JOBZ, 'F'))) {
    INFO.value = -2;
  } else if (!(WNTRES || lsame(JOBR, 'N')) || (WNTRES && !WNTVEC)) {
    INFO.value = -3;
  } else if (!(WNTREF || WNTEX || lsame(JOBF, 'N'))) {
    INFO.value = -4;
  } else if (!((WHTSVD == 1) ||
      (WHTSVD == 2) ||
      (WHTSVD == 3) ||
      (WHTSVD == 4))) {
    INFO.value = -5;
  } else if (M < 0) {
    INFO.value = -6;
  } else if ((N < 0) || (N > M)) {
    INFO.value = -7;
  } else if (LDX < M) {
    INFO.value = -9;
  } else if (LDY < M) {
    INFO.value = -11;
  } else if (!((NRNK == -2) || (NRNK == -1) || ((NRNK >= 1) && (NRNK <= N)))) {
    INFO.value = -12;
  } else if ((TOL < ZERO) || (TOL >= ONE)) {
    INFO.value = -13;
  } else if (LDZ < M) {
    INFO.value = -18;
  } else if ((WNTREF || WNTEX) && (LDB < M)) {
    INFO.value = -21;
  } else if (LDW < N) {
    INFO.value = -23;
  } else if (LDS < N) {
    INFO.value = -25;
  }

  if (INFO.value == 0) {
    // Compute the minimal and the optimal workspace
    // requirements. Simulate running the code and
    // determine minimal and optimal sizes of the
    // workspace at any moment of the run.
    if (N == 0) {
      // Quick return. All output except K is void.
      // INFO=1 signals the void input.
      // In case of a workspace query, the default
      // minimal workspace lengths are returned.
      if (LQUERY) {
        IWORK[1] = 1;
        WORK[1] = 2;
        WORK[2] = 2;
      } else {
        K.value = 0;
      }
      INFO.value = 1;
      return;
    }
    MLWORK = max(2, N);
    OLWORK = max(2, N);
    IMINWR = 1;
    switch (WHTSVD) {
      case 1:
        // The following is specified as the minimal
        // length of WORK in the definition of DGESVD:
        // MWRSVD = max(1,3*min(M,N)+max(M,N),5*min(M,N))
        MWRSVD = max(
            1, max(3 * min(M, N).toInt() + max(M, N).toInt(), 5 * min(M, N)));
        MLWORK = max(MLWORK, N + MWRSVD);
        if (LQUERY) {
          dgesvd(
              'O', 'S', M, N, X, LDX, WORK, B, LDB, W, LDW, RDUMMY, -1, INFO1);
          LWRSVD = max(MWRSVD, RDUMMY[1].toInt());
          OLWORK = max(OLWORK, N + LWRSVD);
        }
        break;
      case 2:
        // The following is specified as the minimal
        // length of WORK in the definition of DGESDD:
        // MWRSDD = 3*min(M,N)*min(M,N) +
        // max( max(M,N),5*min(M,N)*min(M,N)+4*min(M,N) )
        // IMINWR = 8*min(M,N)
        MWRSDD = 3 * min(M, N).toInt() * min(M, N).toInt() +
            max(
                max(M, N),
                5 * min(M, N).toInt() * min(M, N).toInt() +
                    4 * min(M, N).toInt());
        MLWORK = max(MLWORK, N + MWRSDD);
        IMINWR = 8 * min(M, N);
        if (LQUERY) {
          dgesdd('O', M, N, X, LDX, WORK, B, LDB, W, LDW, RDUMMY, -1, IWORK,
              INFO1);
          LWRSDD = max(MWRSDD, RDUMMY[1].toInt());
          OLWORK = max(OLWORK, N + LWRSDD);
        }
        break;
      case 3:
        //LWQP3 = 3*N+1
        //LWORQ = max(N, 1)
        //MWRSVD = max(1,3*min(M,N)+max(M,N),5*min(M,N))
        //MWRSVQ = N + max( LWQP3, MWRSVD, LWORQ ) + max(M,2)
        //MLWORK = N +  MWRSVQ
        //IMINWR = M+N-1
        dgesvdq('H', 'P', 'N', 'R', 'R', M, N, X, LDX, WORK, Z, LDZ, W, LDW,
            NUMRNK, IWORK, LIWORK, RDUMMY, -1, RDUMMY2, -1, INFO1);
        IMINWR = IWORK[1];
        MWRSVQ = RDUMMY[2].toInt();
        MLWORK = max(MLWORK, N + MWRSVQ + RDUMMY2[1].toInt());
        if (LQUERY) {
          LWRSVQ = max(MWRSVQ, RDUMMY[1].toInt());
          OLWORK = max(OLWORK, N + LWRSVQ + RDUMMY2[1].toInt());
        }
        break;
      case 4:
        JSVOPT = 'J';
        //MWRSVJ = max( 7, 2*M+N, 6*N+2*N*N ) ! for JSVOPT='V'
        MWRSVJ = max(max(7, 2 * M + N), max(4 * N + N * N, 2 * N + N * N + 6));
        MLWORK = max(MLWORK, N + MWRSVJ);
        IMINWR = max(3, M + 3 * N);
        if (LQUERY) {
          OLWORK = max(OLWORK, N + MWRSVJ);
        }
        break;
    }
    if (WNTVEC || WNTEX || lsame(JOBZ, 'F')) {
      JOBZL = 'V';
    } else {
      JOBZL = 'N';
    }
    // Workspace calculation to the DGEEV call
    if (lsame(JOBZL, 'V')) {
      MWRKEV = max(1, 4 * N);
    } else {
      MWRKEV = max(1, 3 * N);
    }
    MLWORK = max(MLWORK, N + MWRKEV);
    if (LQUERY) {
      dgeev('N', JOBZL, N, S, LDS, REIG, IMEIG, W, LDW, W, LDW, RDUMMY, -1,
          INFO1);
      LWRKEV = max(MWRKEV, RDUMMY[1].toInt());
      OLWORK = max(OLWORK, N + LWRKEV);
    }

    if (LIWORK < IMINWR && !LQUERY) INFO.value = -29;
    if (LWORK < MLWORK && !LQUERY) INFO.value = -27;
  }

  if (INFO.value != 0) {
    xerbla('DGEDMD', -INFO.value);
    return;
  } else if (LQUERY) {
    // Return minimal and optimal workspace sizes
    IWORK[1] = IMINWR;
    WORK[1] = MLWORK.toDouble();
    WORK[2] = OLWORK.toDouble();
    return;
  }

  OFL = dlamch('O');
  SMALL = dlamch('S');
  BADXY = false;

  // <1> Optional scaling of the snapshots (columns of X, Y)

  if (SCCOLX) {
    // The columns of X will be normalized.
    // To prevent overflows, the column norms of X are
    // carefully computed using DLASSQ.
    K.value = 0;
    for (final i in 1.through(N)) {
      //WORK(i) = dnrm2( M, X(1,i), 1 )
      final SSUM = Box(ONE);
      final SCALE = Box(ZERO);
      dlassq(M, X(1, i).asArray(), 1, SCALE, SSUM);
      if (disnan(SCALE.value) || disnan(SSUM.value)) {
        K.value = 0;
        INFO.value = -8;
        xerbla('DGEDMD', -INFO.value);
      }
      if ((SCALE.value != ZERO) && (SSUM.value != ZERO)) {
        ROOTSC = sqrt(SSUM.value);
        if (SCALE.value >= (OFL / ROOTSC)) {
          // Norm of X(:,i) overflows. First, X(:,i)
          // is scaled by
          // ( ONE / ROOTSC ) / SCALE = 1/||X(:,i)||_2.
          // Next, the norm of X(:,i) is stored without
          // overflow as WORK(i) = - SCALE * (ROOTSC/M),
          // the minus sign indicating the 1/M factor.
          // Scaling is performed without overflow, and
          // underflow may occur in the smallest entries
          // of X(:,i). The relative backward and forward
          // errors are small in the ell_2 norm.
          dlascl('G', 0, 0, SCALE.value, ONE / ROOTSC, M, 1, X(1, i), M, INFO2);
          WORK[i] = -SCALE.value * (ROOTSC / M);
        } else {
          // X(:,i) will be scaled to unit 2-norm
          WORK[i] = SCALE.value * ROOTSC;
          dlascl('G', 0, 0, WORK[i], ONE, M, 1, X(1, i), M, INFO2);
          // X(1:M,i) = (ONE/WORK(i)) * X(1:M,i)
        }
      } else {
        WORK[i] = ZERO;
        K.value++;
      }
    }
    if (K.value == N) {
      // All columns of X are zero. Return error code -8.
      // (the 8th input variable had an illegal value)
      K.value = 0;
      INFO.value = -8;
      xerbla('DGEDMD', -INFO.value);
      return;
    }
    for (final i in 1.through(N)) {
      // Now, apply the same scaling to the columns of Y.
      if (WORK[i] > ZERO) {
        dscal(M, ONE / WORK[i], Y(1, i).asArray(), 1);
        // Y(1:M,i) = (ONE/WORK(i)) * Y(1:M,i)
      } else if (WORK[i] < ZERO) {
        dlascl('G', 0, 0, -WORK[i], ONE / M, M, 1, Y(1, i), M, INFO2);
      } else if (Y[idamax(M, Y(1, i).asArray(), 1)][i] != ZERO) {
        // X(:,i) is zero vector. For consistency,
        // Y(:,i) should also be zero. If Y(:,i) is not
        // zero, then the data might be inconsistent or
        // corrupted. If JOBS == 'C', Y(:,i) is set to
        // zero and a warning flag is raised.
        // The computation continues but the
        // situation will be reported in the output.
        BADXY = true;
        if (lsame(JOBS, 'C')) {
          dscal(M, ZERO, Y(1, i).asArray(), 1);
        }
      }
    }
  }

  if (SCCOLY) {
    // The columns of Y will be normalized.
    // To prevent overflows, the column norms of Y are
    // carefully computed using DLASSQ.
    for (final i in 1.through(N)) {
      // WORK(i) = dnrm2( M, Y(1,i), 1 )
      final SSUM = Box(ONE);
      final SCALE = Box(ZERO);
      dlassq(M, Y(1, i).asArray(), 1, SCALE, SSUM);
      if (disnan(SCALE.value) || disnan(SSUM.value)) {
        K.value = 0;
        INFO.value = -10;
        xerbla('DGEDMD', -INFO.value);
      }
      if (SCALE.value != ZERO && (SSUM.value != ZERO)) {
        ROOTSC = sqrt(SSUM.value);
        if (SCALE.value >= (OFL / ROOTSC)) {
          // Norm of Y(:,i) overflows. First, Y(:,i)
          // is scaled by
          // ( ONE / ROOTSC ) / SCALE = 1/||Y(:,i)||_2.
          // Next, the norm of Y(:,i) is stored without
          // overflow as WORK(i) = - SCALE * (ROOTSC/M),
          // the minus sign indicating the 1/M factor.
          // Scaling is performed without overflow, and
          // underflow may occur in the smallest entries
          // of Y(:,i). The relative backward and forward
          // errors are small in the ell_2 norm.
          dlascl('G', 0, 0, SCALE.value, ONE / ROOTSC, M, 1, Y(1, i), M, INFO2);
          WORK[i] = -SCALE.value * (ROOTSC / M);
        } else {
          // X(:,i) will be scaled to unit 2-norm
          WORK[i] = SCALE.value * ROOTSC;
          dlascl('G', 0, 0, WORK[i], ONE, M, 1, Y(1, i), M, INFO2);
          // Y(1:M,i) = (ONE/WORK[i]) * Y(1:M,i)
        }
      } else {
        WORK[i] = ZERO;
      }
    }
    for (final i in 1.through(N)) {
      // Now, apply the same scaling to the columns of X.
      if (WORK[i] > ZERO) {
        dscal(M, ONE / WORK[i], X(1, i).asArray(), 1);
        // X(1:M,i) = (ONE/WORK(i)) * X(1:M,i)
      } else if (WORK[i] < ZERO) {
        dlascl('G', 0, 0, -WORK[i], ONE / M, M, 1, X(1, i), M, INFO2);
      } else if (X[idamax(M, X(1, i).asArray(), 1)][i] != ZERO) {
        // Y(:,i) is zero vector.  If X(:,i) is not
        // zero, then a warning flag is raised.
        // The computation continues but the
        // situation will be reported in the output.
        BADXY = true;
      }
    }
  }

  // <2> SVD of the data snapshot matrix X.
  // =====================================
  // The left singular vectors are stored in the array X.
  // The right singular vectors are in the array W.
  // The array W will later on contain the eigenvectors
  // of a Rayleigh quotient.
  NUMRNK.value = N;
  switch (WHTSVD) {
    case 1:
      dgesvd('O', 'S', M, N, X, LDX, WORK, B, LDB, W, LDW, WORK(N + 1),
          LWORK - N, INFO1);
      T_OR_N = 'T';
      break;
    case 2:
      dgesdd('O', M, N, X, LDX, WORK, B, LDB, W, LDW, WORK(N + 1), LWORK - N,
          IWORK, INFO1);
      T_OR_N = 'T';
      break;
    case 3:
      dgesvdq(
          'H',
          'P',
          'N',
          'R',
          'R',
          M,
          N,
          X,
          LDX,
          WORK,
          Z,
          LDZ,
          W,
          LDW,
          NUMRNK,
          IWORK,
          LIWORK,
          WORK(N + max(2, M).toInt() + 1),
          LWORK - N - max(2, M),
          WORK(N + 1),
          max(2, M),
          INFO1);
      dlacpy('A', M, NUMRNK.value, Z, LDZ, X, LDX);
      T_OR_N = 'T';
      break;
    case 4:
      dgejsv('F', 'U', JSVOPT, 'N', 'N', 'P', M, N, X, LDX, WORK, Z, LDZ, W,
          LDW, WORK(N + 1), LWORK - N, IWORK, INFO1);
      dlacpy('A', M, N, Z, LDZ, X, LDX);
      T_OR_N = 'N';
      XSCL1 = WORK[N + 1];
      XSCL2 = WORK[N + 2];
      if (XSCL1 != XSCL2) {
        /// This is an exceptional situation. If the
        /// data matrices are not scaled and the
        /// largest singular value of X overflows.
        /// In that case DGEJSV can return the SVD
        /// in scaled form. The scaling factor can be used
        /// to rescale the data (X and Y).
        dlascl('G', 0, 0, XSCL1, XSCL2, M, N, Y, LDY, INFO2);
      }
      break;
  }

  if (INFO1.value > 0) {
    // The SVD selected subroutine did not converge.
    // Return with an error code.
    INFO.value = 2;
    return;
  }

  if (WORK[1] == ZERO) {
    // The largest computed singular value of (scaled)
    // X is zero. Return error code -8
    // (the 8th input variable had an illegal value).
    K.value = 0;
    INFO.value = -8;
    xerbla('DGEDMD', -INFO.value);
    return;
  }

  //<3> Determine the numerical rank of the data
  //    snapshots matrix X. This depends on the
  //    parameters NRNK and TOL.

  switch (NRNK) {
    case -1:
      K.value = 1;
      for (final i in 2.through(NUMRNK.value)) {
        if ((WORK[i] <= WORK[1] * TOL) || (WORK[i] <= SMALL)) break;
        K.value++;
      }
    case -2:
      K.value = 1;
      for (final i in 1.through(NUMRNK.value - 1)) {
        if ((WORK[i + 1] <= WORK[i] * TOL) || (WORK[i] <= SMALL)) break;
        K.value++;
      }
    default:
      K.value = 1;
      for (final i in 2.through(NRNK)) {
        if (WORK[i] <= SMALL) break;
        K.value++;
      }
  }
  //   Now, U = X(1:M,1:K) is the SVD/POD basis for the
  //   snapshot data in the input matrix X.

  //<4> Compute the Rayleigh quotient S = U^T * A * U.
  //    Depending on the requested outputs, the computation
  //    is organized to compute additional auxiliary
  //    matrices (for the residuals and refinements).
  //
  //    In all formulas below, we need V_k*Sigma_k^(-1)
  //    where either V_k is in W(1:N,1:K), or V_k^T is in
  //    W(1:K,1:N). Here Sigma_k=diag(WORK(1:K)).
  if (lsame(T_OR_N, 'N')) {
    for (final i in 1.through(K.value)) {
      dscal(N, ONE / WORK[i], W(1, i).asArray(), 1);
      // W(1:N,i) = (ONE/WORK[i]) * W(1:N,i);
    }
  } else {
    // This non-unit stride access is due to the fact
    // that DGESVD, DGESVDQ and DGESDD return the
    // transposed matrix of the right singular vectors.
    //DO i = 1, K
    // dscal( N, ONE/WORK[i], W(i,1), LDW )
    // // W(i,1:N) = (ONE/WORK[i]) * W(i,1:N)
    for (final i in 1.through(K.value)) {
      WORK[N + i] = ONE / WORK[i];
    }
    for (final j in 1.through(N)) {
      for (final i in 1.through(K.value)) {
        W[i][j] = WORK[N + i] * W[i][j];
      }
    }
  }

  if (WNTREF) {
    //
    // Need A*U(:,1:K)=Y*V_k*inv(diag(WORK(1:K)))
    // for computing the refined Ritz vectors
    // (optionally, outside DGEDMD).
    dgemm('N', T_OR_N, M, K.value, N, ONE, Y, LDY, W, LDW, ZERO, Z, LDZ);
    // Z(1:M,1:K)=MATMUL(Y(1:M,1:N),TRANSPOSE(W(1:K,1:N)))  , for T_OR_N=='T'
    // Z(1:M,1:K)=MATMUL(Y(1:M,1:N),W(1:N,1:K))             , for T_OR_N=='N'
    //
    // At this point Z contains
    // A * U(:,1:K) = Y * V_k * Sigma_k^(-1), and
    // this is needed for computing the residuals.
    // This matrix is  returned in the array B and
    // it can be used to compute refined Ritz vectors.
    dlacpy('A', M, K.value, Z, LDZ, B, LDB);
    // B(1:M,1:K) = Z(1:M,1:K)

    dgemm('T', 'N', K.value, K.value, M, ONE, X, LDX, Z, LDZ, ZERO, S, LDS);
    // S(1:K,1:K) = MATMUL(TANSPOSE(X(1:M,1:K)),Z(1:M,1:K))
    // At this point S = U^T * A * U is the Rayleigh quotient.
  } else {
    // A * U(:,1:K) is not explicitly needed and the
    // computation is organized differently. The Rayleigh
    // quotient is computed more efficiently.
    dgemm('T', 'N', K.value, N, M, ONE, X, LDX, Y, LDY, ZERO, Z, LDZ);
    // Z(1:K,1:N) = MATMUL( TRANSPOSE(X(1:M,1:K)), Y(1:M,1:N) )
    // In the two DGEMM calls here, can use K for LDZ.
    dgemm('N', T_OR_N, K.value, K.value, N, ONE, Z, LDZ, W, LDW, ZERO, S, LDS);
    // S(1:K,1:K) = MATMUL(Z(1:K,1:N),TRANSPOSE(W(1:K,1:N))) , for T_OR_N=='T'
    // S(1:K,1:K) = MATMUL(Z(1:K,1:N),(W(1:N,1:K)))          , for T_OR_N=='N'
    // At this point S = U^T * A * U is the Rayleigh quotient.
    // If the residuals are requested, save scaled V_k into Z.
    // Recall that V_k or V_k^T is stored in W.
    if (WNTRES || WNTEX) {
      if (lsame(T_OR_N, 'N')) {
        dlacpy('A', N, K.value, W, LDW, Z, LDZ);
      } else {
        dlacpy('A', K.value, N, W, LDW, Z, LDZ);
      }
    }
  }

  //<5> Compute the Ritz values and (if requested) the
  //   right eigenvectors of the Rayleigh quotient.

  dgeev('N', JOBZL, K.value, S, LDS, REIG, IMEIG, W, LDW, W, LDW, WORK(N + 1),
      LWORK - N, INFO1);

  // W(1:K,1:K) contains the eigenvectors of the Rayleigh
  // quotient. Even in the case of complex spectrum, all
  // computation is done in real arithmetic. REIG and
  // IMEIG are the real and the imaginary parts of the
  // eigenvalues, so that the spectrum is given as
  // REIG(:) + sqrt(-1)*IMEIG(:). Complex conjugate pairs
  // are listed at consecutive positions. For such a
  // complex conjugate pair of the eigenvalues, the
  // corresponding eigenvectors are also a complex
  // conjugate pair with the real and imaginary parts
  // stored column-wise in W at the corresponding
  // consecutive column indices. See the description of Z.
  // Also, see the description of DGEEV.
  if (INFO1.value > 0) {
    // DGEEV failed to compute the eigenvalues and
    // eigenvectors of the Rayleigh quotient.
    INFO.value = 3;
    return;
  }

  // <6> Compute the eigenvectors (if requested) and,
  // the residuals (if requested).

  if (WNTVEC || WNTEX) {
    if (WNTRES) {
      if (WNTREF) {
        // Here, if the refinement is requested, we have
        // A*U(:,1:K) already computed and stored in Z.
        // For the residuals, need Y = A * U(:,1;K) * W.
        dgemm('N', 'N', M, K.value, K.value, ONE, Z, LDZ, W, LDW, ZERO, Y, LDY);
        // Y(1:M,1:K) = Z(1:M,1:K) * W(1:K,1:K)
        // This frees Z; Y contains A * U(:,1:K) * W.
      } else {
        // Compute S = V_k * Sigma_k^(-1) * W, where
        // V_k * Sigma_k^(-1) is stored in Z
        dgemm(T_OR_N, 'N', N, K.value, K.value, ONE, Z, LDZ, W, LDW, ZERO, S,
            LDS);
        // Then, compute Z = Y * S =
        // = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
        // = A * U(:,1:K) * W(1:K,1:K)
        dgemm('N', 'N', M, K.value, N, ONE, Y, LDY, S, LDS, ZERO, Z, LDZ);
        // Save a copy of Z into Y and free Z for holding
        // the Ritz vectors.
        dlacpy('A', M, K.value, Z, LDZ, Y, LDY);
        if (WNTEX) dlacpy('A', M, K.value, Z, LDZ, B, LDB);
      }
    } else if (WNTEX) {
      // Compute S = V_k * Sigma_k^(-1) * W, where
      // V_k * Sigma_k^(-1) is stored in Z
      dgemm(
          T_OR_N, 'N', N, K.value, K.value, ONE, Z, LDZ, W, LDW, ZERO, S, LDS);
      // Then, compute Z = Y * S =
      // = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
      // = A * U(:,1:K) * W(1:K,1:K)
      dgemm('N', 'N', M, K.value, N, ONE, Y, LDY, S, LDS, ZERO, B, LDB);
      // The above replaces the following two calls
      // that were used in the developing-testing phase.
      // dgemm( 'N', 'N', M, K, N, ONE, Y, LDY, S,             //           LDS, ZERO, Z, LDZ)
      // Save a copy of Z into B and free Z for holding
      // the Ritz vectors.
      // dlacpy( 'A', M, K, Z, LDZ, B, LDB )
    }

    // Compute the real form of the Ritz vectors
    if (WNTVEC) {
      dgemm('N', 'N', M, K.value, K.value, ONE, X, LDX, W, LDW, ZERO, Z, LDZ);
      // Z(1:M,1:K) = MATMUL(X(1:M,1:K), W(1:K,1:K))
    }

    if (WNTRES) {
      i = 1;
      while (i <= K.value) {
        if (IMEIG[i] == ZERO) {
          // have a real eigenvalue with real eigenvector
          daxpy(M, -REIG[i], Z(1, i).asArray(), 1, Y(1, i).asArray(), 1);
          // Y(1:M,i) = Y(1:M,i) - REIG[i] * Z(1:M,i)
          RES[i] = dnrm2(M, Y(1, i).asArray(), 1);
          i = i + 1;
        } else {
          // Have a complex conjugate pair
          // REIG[i] +- sqrt(-1)*IMEIG[i].
          // Since all computation is done in real
          // arithmetic, the formula for the residual
          // is recast for real representation of the
          // complex conjugate eigenpair. See the
          // description of RES.
          AB[1][1] = REIG[i];
          AB[2][1] = -IMEIG[i];
          AB[1][2] = IMEIG[i];
          AB[2][2] = REIG[i];
          dgemm(
              'N', 'N', M, 2, 2, -ONE, Z(1, i), LDZ, AB, 2, ONE, Y(1, i), LDY);
          // Y(1:M,i:i+1) = Y(1:M,i:i+1) - Z(1:M,i:i+1) * AB
          RES[i] = dlange('F', M, 2, Y(1, i), LDY, WORK(N + 1));
          RES[i + 1] = RES[i];
          i = i + 2;
        }
      }
    }
  }

  if (WHTSVD == 4) {
    WORK[N + 1] = XSCL1;
    WORK[N + 2] = XSCL2;
  }

  // Successful exit.
  if (!BADXY) {
    INFO.value = 0;
  } else {
    // A warning on possible data inconsistency.
    // This should be a rare event.
    INFO.value = 4;
  }
}
