// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:collection/collection.dart';
import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/range.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgeev.dart';
import 'package:lapack/src/zgejsv.dart';
import 'package:lapack/src/zgesdd.dart';
import 'package:lapack/src/zgesvd.dart';
import 'package:lapack/src/zgesvdq.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlassq.dart';

void zgedmd(
  final String JOBS,
  final String JOBZ,
  final String JOBR,
  final String JOBF,
  final int WHTSVD,
  final int M,
  final int N,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> Y_,
  final int LDY,
  final int NRNK,
  final double TOL,
  final Box<int> K,
  final Array<Complex> EIGS_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<double> RES_,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> W_,
  final int LDW,
  final Matrix<Complex> S_,
  final int LDS,
  final Array<Complex> ZWORK_,
  final int LZWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
//  -- LAPACK driver routine                                           --

//  -- LAPACK is a software package provided by University of          --
//  -- Tennessee, University of California Berkeley, University of     --
//  -- Colorado Denver and NAG Ltd..                                   --
  final X = X_.having(ld: LDX);
  final Y = Y_.having(ld: LDY);
  final EIGS = EIGS_.having();
  final Z = Z_.having(ld: LDZ);
  final RES = RES_.having();
  final B = B_.having(ld: LDB);
  final W = W_.having(ld: LDW);
  final S = S_.having(ld: LDS);
  final ZWORK = ZWORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

  const ONE = 1.0, ZERO = 0.0;
  double OFL, ROOTSC, SMALL, XSCL1 = 0, XSCL2 = 0;
  int IMINWR = 0,
      LWRKEV,
      LWRSDD,
      LWRSVD,
      LWRSVJ,
      LWRSVQ,
      MLWORK = 0,
      MWRKEV,
      MWRSDD,
      MWRSVD,
      MWRSVJ,
      MWRSVQ,
      OLWORK = 0,
      MLRWRK = 0;
  bool BADXY, LQUERY, SCCOLX, SCCOLY, WNTEX, WNTREF, WNTRES, WNTVEC;
  String JOBZL = '', T_OR_N = '';
  String JSVOPT = '';
  final RDUMMY = Array<double>(2);
  final INFO1 = Box(0), INFO2 = Box(0), NUMRNK = Box(0);
  final SCALE = Box(ZERO), SSUM = Box(ZERO);

  // Test the input arguments

  WNTRES = lsame(JOBR, 'R');
  SCCOLX = lsame(JOBS, 'S') || lsame(JOBS, 'C');
  SCCOLY = lsame(JOBS, 'Y');
  WNTVEC = lsame(JOBZ, 'V');
  WNTREF = lsame(JOBF, 'R');
  WNTEX = lsame(JOBF, 'E');
  INFO.value = 0;
  LQUERY = ((LZWORK == -1) || (LIWORK == -1) || (LRWORK == -1));

  if (!(SCCOLX || SCCOLY || lsame(JOBS, 'N'))) {
    INFO.value = -1;
  } else if (!(WNTVEC || lsame(JOBZ, 'N') || lsame(JOBZ, 'F'))) {
    INFO.value = -2;
  } else if (!(WNTRES || lsame(JOBR, 'N')) || (WNTRES && (!WNTVEC))) {
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
    INFO.value = -17;
  } else if ((WNTREF || WNTEX) && (LDB < M)) {
    INFO.value = -20;
  } else if (LDW < N) {
    INFO.value = -22;
  } else if (LDS < N) {
    INFO.value = -24;
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
        RWORK[1] = 1;
        ZWORK[1] = 2.toComplex();
        ZWORK[2] = 2.toComplex();
      } else {
        K.value = 0;
      }
      INFO.value = 1;
      return;
    }

    IMINWR = 1;
    MLRWRK = max(1, N);
    MLWORK = 2;
    OLWORK = 2;
    switch (WHTSVD) {
      case 1:
        // The following is specified as the minimal
        // length of WORK in the definition of zgesvd:
        // MWRSVD = max(1,2*min(M,N)+max(M,N))
        MWRSVD = max(1, 2 * min(M, N) + max(M, N)).toInt();
        MLWORK = max(MLWORK, MWRSVD);
        MLRWRK = max(MLRWRK, N + 5 * min(M, N));
        if (LQUERY) {
          zgesvd('O', 'S', M, N, X, LDX, RWORK, B, LDB, W, LDW, ZWORK, -1,
              RDUMMY, INFO1);
          LWRSVD = ZWORK[1].toInt();
          OLWORK = max(OLWORK, LWRSVD);
        }
      case 2:
        // The following is specified as the minimal
        // length of WORK in the definition of zgesdd:
        // MWRSDD = 2*min(M,N)*min(M,N)+2*min(M,N)+max(M,N).
        // RWORK length: 5*min(M,N)*min(M,N)+7*min(M,N)
        // In LAPACK 3.10.1 RWORK is defined differently.
        // Below we take max over the two versions.
        // IMINWR = 8*min(M,N)
        MWRSDD =
            (2 * min(M, N) * min(M, N) + 2 * min(M, N) + max(M, N)).toInt();
        MLWORK = max(MLWORK, MWRSDD);
        IMINWR = 8 * min(M, N);
        MLRWRK = max(
                MLRWRK,
                N +
                    [
                      5 * min(M, N) * min(M, N) + 7 * min(M, N),
                      5 * min(M, N) * min(M, N) + 5 * min(M, N),
                      2 * max(M, N) * min(M, N) +
                          2 * min(M, N) * min(M, N) +
                          min(M, N)
                    ].max)
            .toInt();
        if (LQUERY) {
          zgesdd('O', M, N, X, LDX, RWORK, B, LDB, W, LDW, ZWORK, -1, RDUMMY,
              IWORK, INFO1);
          LWRSDD = max(MWRSDD, ZWORK[1].toInt());
          // Possible bug in zgesdd optimal workspace size.
          OLWORK = max(OLWORK, LWRSDD);
        }
      case 3:
        zgesvdq('H', 'P', 'N', 'R', 'R', M, N, X, LDX, RWORK, Z, LDZ, W, LDW,
            NUMRNK, IWORK, -1, ZWORK, -1, RDUMMY, -1, INFO1);
        IMINWR = IWORK[1];
        MWRSVQ = ZWORK[2].toInt();
        MLWORK = max(MLWORK, MWRSVQ);
        MLRWRK = max(MLRWRK, N + RDUMMY[1].toInt());
        if (LQUERY) {
          LWRSVQ = ZWORK[1].toInt();
          OLWORK = max(OLWORK, LWRSVQ);
        }
      case 4:
        JSVOPT = 'J';
        zgejsv('F', 'U', JSVOPT, 'R', 'N', 'P', M, N, X, LDX, RWORK, Z, LDZ, W,
            LDW, ZWORK, -1, RDUMMY, -1, IWORK, INFO1);
        IMINWR = IWORK[1];
        MWRSVJ = ZWORK[2].toInt();
        MLWORK = max(MLWORK, MWRSVJ);
        MLRWRK = max(MLRWRK, N + max(7, RDUMMY[1].toInt()));
        if (LQUERY) {
          LWRSVJ = ZWORK[1].toInt();
          OLWORK = max(OLWORK, LWRSVJ);
        }
    }
    if (WNTVEC || WNTEX || lsame(JOBZ, 'F')) {
      JOBZL = 'V';
    } else {
      JOBZL = 'N';
    }
    // Workspace calculation to the zgeev call
    MWRKEV = max(1, 2 * N);
    MLWORK = max(MLWORK, MWRKEV);
    MLRWRK = max(MLRWRK, N + 2 * N);
    if (LQUERY) {
      zgeev(
          'N', JOBZL, N, S, LDS, EIGS, W, LDW, W, LDW, ZWORK, -1, RWORK, INFO1);
      LWRKEV = ZWORK[1].toInt();
      OLWORK = max(OLWORK, LWRKEV);
    }

    if (LIWORK < IMINWR && (!LQUERY)) INFO.value = -30;
    if (LRWORK < MLRWRK && (!LQUERY)) INFO.value = -28;
    if (LZWORK < MLWORK && (!LQUERY)) INFO.value = -26;
  }

  if (INFO.value != 0) {
    xerbla('ZGEDMD', -INFO.value);
    return;
  } else if (LQUERY) {
    // Return minimal and optimal workspace sizes
    IWORK[1] = IMINWR;
    RWORK[1] = MLRWRK.toDouble();
    ZWORK[1] = MLWORK.toComplex();
    ZWORK[2] = OLWORK.toComplex();
    return;
  }

  OFL = dlamch('O');
  SMALL = dlamch('S');
  BADXY = false;

  // <1> Optional scaling of the snapshots (columns of X, Y)
  // ==========================================================
  if (SCCOLX) {
    // The columns of X will be normalized.
    // To prevent overflows, the column norms of X are
    // carefully computed using zlassq.
    K.value = 0;
    for (final i in 1.through(N)) {
      //WORK(i) = dznrm2( M, X(1,i), 1 )
      SSUM.value = ONE;
      SCALE.value = ZERO;
      zlassq(M, X(1, i).asArray(), 1, SCALE, SSUM);
      if (disnan(SCALE.value) || disnan(SSUM.value)) {
        K.value = 0;
        INFO.value = -8;
        xerbla('ZGEDMD', -INFO.value);
      }
      if ((SCALE.value != ZERO) && (SSUM.value != ZERO)) {
        ROOTSC = sqrt(SSUM.value);
        if (SCALE.value >= (OFL / ROOTSC)) {
          // Norm of X(:,i) overflows. First, X(:,i)
          // is scaled by
          // ( ONE / ROOTSC ) / SCALE = 1/||X(:,i)||_2.
          // Next, the norm of X(:,i) is stored without
          // overflow as RWORK(i) = - SCALE * (ROOTSC/M),
          // the minus sign indicating the 1/M factor.
          // Scaling is performed without overflow, and
          // underflow may occur in the smallest entries
          // of X(:,i). The relative backward and forward
          // errors are small in the ell_2 norm.
          zlascl(
              'G', 0, 0, SCALE.value, ONE / ROOTSC, M, 1, X(1, i), LDX, INFO2);
          RWORK[i] = -SCALE.value * (ROOTSC / M);
        } else {
          // X(:,i) will be scaled to unit 2-norm
          RWORK[i] = SCALE.value * ROOTSC;
          zlascl('G', 0, 0, RWORK[i], ONE, M, 1, X(1, i), LDX, INFO2);
          // X(1:M,i) = (ONE/RWORK[i]) * X(1:M,i)   ! INTRINSIC
        }
      } else {
        RWORK[i] = ZERO;
        K.value = K.value + 1;
      }
    }
    if (K.value == N) {
      // All columns of X are zero. Return error code -8.
      // (the 8th input variable had an illegal value)
      K.value = 0;
      INFO.value = -8;
      xerbla('ZGEDMD', -INFO.value);
      return;
    }
    for (final i in 1.through(N)) {
      // Now, apply the same scaling to the columns of Y.
      if (RWORK[i] > ZERO) {
        zdscal(M, ONE / RWORK[i], Y(1, i).asArray(), 1);
        // Y(1:M,i) = (ONE/RWORK[i]) * Y(1:M,i)       ! INTRINSIC
      } else if (RWORK[i] < ZERO) {
        zlascl('G', 0, 0, -RWORK[i], ONE / M, M, 1, Y(1, i), LDY, INFO2);
      } else if (Y[izamax(M, Y(1, i).asArray(), 1)][i].abs() != ZERO) {
        // X(:,i) is zero vector. For consistency,
        // Y(:,i) should also be zero. If Y(:,i) is not
        // zero, then the data might be inconsistent or
        // corrupted. If JOBS == 'C', Y(:,i) is set to
        // zero and a warning flag is raised.
        // The computation continues but the
        // situation will be reported in the output.
        BADXY = true;
        if (lsame(JOBS, 'C')) zdscal(M, ZERO, Y(1, i).asArray(), 1);
      }
    }
  }

  if (SCCOLY) {
    // The columns of Y will be normalized.
    // To prevent overflows, the column norms of Y are
    // carefully computed using zlassq.
    for (final i in 1.through(N)) {
      //RWORK(i) = dznrm2( M, Y(1,i), 1 )
      SSUM.value = ONE;
      SCALE.value = ZERO;
      zlassq(M, Y(1, i).asArray(), 1, SCALE, SSUM);
      if (disnan(SCALE.value) || disnan(SSUM.value)) {
        K.value = 0;
        INFO.value = -10;
        xerbla('ZGEDMD', -INFO.value);
      }
      if (SCALE.value != ZERO && (SSUM.value != ZERO)) {
        ROOTSC = sqrt(SSUM.value);
        if (SCALE.value >= (OFL / ROOTSC)) {
          // Norm of Y(:,i) overflows. First, Y(:,i)
          // is scaled by
          // ( ONE / ROOTSC ) / SCALE = 1/||Y(:,i)||_2.
          // Next, the norm of Y(:,i) is stored without
          // overflow as RWORK(i) = - SCALE * (ROOTSC/M),
          // the minus sign indicating the 1/M factor.
          // Scaling is performed without overflow, and
          // underflow may occur in the smallest entries
          // of Y(:,i). The relative backward and forward
          // errors are small in the ell_2 norm.
          zlascl(
              'G', 0, 0, SCALE.value, ONE / ROOTSC, M, 1, Y(1, i), LDY, INFO2);
          RWORK[i] = -SCALE.value * (ROOTSC / M);
        } else {
          // Y(:,i) will be scaled to unit 2-norm
          RWORK[i] = SCALE.value * ROOTSC;
          zlascl('G', 0, 0, RWORK[i], ONE, M, 1, Y(1, i), LDY, INFO2);
          // Y(1:M,i) = (ONE/RWORK[i]) * Y(1:M,i)          ! INTRINSIC
        }
      } else {
        RWORK[i] = ZERO;
      }
    }
    for (final i in 1.through(N)) {
      // Now, apply the same scaling to the columns of X.
      if (RWORK[i] > ZERO) {
        zdscal(M, ONE / RWORK[i], X(1, i).asArray(), 1);
        // X(1:M,i) = (ONE/RWORK[i]) * X(1:M,i)      ! INTRINSIC
      } else if (RWORK[i] < ZERO) {
        zlascl('G', 0, 0, -RWORK[i], ONE / M, M, 1, X(1, i), LDX, INFO2);
      } else if (X[izamax(M, X(1, i).asArray(), 1)][i].abs() != ZERO) {
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
      zgesvd('O', 'S', M, N, X, LDX, RWORK, B, LDB, W, LDW, ZWORK, LZWORK,
          RWORK(N + 1), INFO1);
      T_OR_N = 'C';
    case 2:
      zgesdd('O', M, N, X, LDX, RWORK, B, LDB, W, LDW, ZWORK, LZWORK,
          RWORK(N + 1), IWORK, INFO1);
      T_OR_N = 'C';
    case 3:
      zgesvdq(
          'H',
          'P',
          'N',
          'R',
          'R',
          M,
          N,
          X,
          LDX,
          RWORK,
          Z,
          LDZ,
          W,
          LDW,
          NUMRNK,
          IWORK,
          LIWORK,
          ZWORK,
          LZWORK,
          RWORK(N + 1),
          LRWORK - N,
          INFO1);
      zlacpy('A', M, NUMRNK.value, Z, LDZ, X, LDX);
      T_OR_N = 'C';
    case 4:
      zgejsv('F', 'U', JSVOPT, 'R', 'N', 'P', M, N, X, LDX, RWORK, Z, LDZ, W,
          LDW, ZWORK, LZWORK, RWORK(N + 1), LRWORK - N, IWORK, INFO1);
      zlacpy('A', M, N, Z, LDZ, X, LDX);
      T_OR_N = 'N';
      XSCL1 = RWORK[N + 1];
      XSCL2 = RWORK[N + 2];
      if (XSCL1 != XSCL2) {
        // This is an exceptional situation. If the
        // data matrices are not scaled and the
        // largest singular value of X overflows.
        // In that case zgejsv can return the SVD
        // in scaled form. The scaling factor can be used
        // to rescale the data (X and Y).
        zlascl('G', 0, 0, XSCL1, XSCL2, M, N, Y, LDY, INFO2);
      }
  }

  if (INFO1.value > 0) {
    // The SVD selected subroutine did not converge.
    // Return with an error code.
    INFO.value = 2;
    return;
  }

  if (RWORK[1] == ZERO) {
    // The largest computed singular value of (scaled)
    // X is zero. Return error code -8
    // (the 8th input variable had an illegal value).
    K.value = 0;
    INFO.value = -8;
    xerbla('ZGEDMD', -INFO.value);
    return;
  }

  //<3> Determine the numerical rank of the data
  //    snapshots matrix X. This depends on the
  //    parameters NRNK and TOL.

  switch (NRNK) {
    case -1:
      K.value = 1;
      for (final i in 2.through(NUMRNK.value)) {
        if ((RWORK[i] <= RWORK[1] * TOL) || (RWORK[i] <= SMALL)) break;
        K.value = K.value + 1;
      }
    case -2:
      K.value = 1;
      for (final i in 1.through(NUMRNK.value - 1)) {
        if ((RWORK[i + 1] <= RWORK[i] * TOL) || (RWORK[i] <= SMALL)) break;
        K.value = K.value + 1;
      }
    default:
      K.value = 1;
      for (final i in 2.through(NRNK)) {
        if (RWORK[i] <= SMALL) break;
        K.value = K.value + 1;
      }
  }
  // Now, U = X(1:M,1:K) is the SVD/POD basis for the
  // snapshot data in the input matrix X.

  //<4> Compute the Rayleigh quotient S = U^H * A * U.
  // Depending on the requested outputs, the computation
  // is organized to compute additional auxiliary
  // matrices (for the residuals and refinements).
  //
  // In all formulas below, we need V_k*Sigma_k^(-1)
  // where either V_k is in W(1:N,1:K), or V_k^H is in
  // W(1:K,1:N). Here Sigma_k=diag(WORK(1:K)).
  if (lsame(T_OR_N, 'N')) {
    for (final i in 1.through(K.value)) {
      zdscal(N, ONE / RWORK[i], W(1, i).asArray(), 1);
      // W(1:N,i) = (ONE/RWORK[i]) * W(1:N,i)      ! INTRINSIC
    }
  } else {
    // This non-unit stride access is due to the fact
    // that zgesvd, zgesvdq and zgesdd return the
    // adjoint matrix of the right singular vectors.
    //DO i = 1, K
    // zdscal( N, ONE/RWORK[i], W(i,1), LDW )    ;
    // ! W(i,1:N) = (ONE/RWORK[i]) * W(i,1:N)      ! INTRINSIC
    //}
    for (final i in 1.through(K.value)) {
      RWORK[N + i] = ONE / RWORK[i];
    }
    for (final j in 1.through(N)) {
      for (final i in 1.through(K.value)) {
        W[i][j] = RWORK[N + i].toComplex() * W[i][j];
      }
    }
  }

  if (WNTREF) {
    //
    // Need A*U(:,1:K)=Y*V_k*inv(diag(WORK(1:K)))
    // for computing the refined Ritz vectors
    // (optionally, outside ZGEDMD).
    zgemm('N', T_OR_N, M, K.value, N, Complex.one, Y, LDY, W, LDW, Complex.zero,
        Z, LDZ);
    // Z(1:M,1:K)=MATMUL(Y(1:M,1:N),TRANSPOSE(CONJG(W(1:K,1:N)))) ! INTRINSIC, for T_OR_N=='C'
    // Z(1:M,1:K)=MATMUL(Y(1:M,1:N),W(1:N,1:K))                   ! INTRINSIC, for T_OR_N=='N'
    //
    // At this point Z contains
    // A * U(:,1:K) = Y * V_k * Sigma_k^(-1), and
    // this is needed for computing the residuals.
    // This matrix is  returned in the array B and
    // it can be used to compute refined Ritz vectors.
    zlacpy('A', M, K.value, Z, LDZ, B, LDB);
    // B(1:M,1:K) = Z(1:M,1:K)                  ! INTRINSIC

    zgemm('C', 'N', K.value, K.value, M, Complex.one, X, LDX, Z, LDZ,
        Complex.zero, S, LDS);
    // S(1:K,1:K) = MATMUL(TRANSPOSE(CONJG(X(1:M,1:K))),Z(1:M,1:K)) ! INTRINSIC
    // At this point S = U^H * A * U is the Rayleigh quotient.
  } else {
    // A * U(:,1:K) is not explicitly needed and the
    // computation is organized differently. The Rayleigh
    // quotient is computed more efficiently.
    zgemm('C', 'N', K.value, N, M, Complex.one, X, LDX, Y, LDY, Complex.zero, Z,
        LDZ);
    // Z(1:K,1:N) = MATMUL( TRANSPOSE(CONJG(X(1:M,1:K))), Y(1:M,1:N) )  ! INTRINSIC
    //
    zgemm('N', T_OR_N, K.value, K.value, N, Complex.one, Z, LDZ, W, LDW,
        Complex.zero, S, LDS);
    // S(1:K,1:K) = MATMUL(Z(1:K,1:N),TRANSPOSE(CONJG(W(1:K,1:N)))) ! INTRINSIC, for T_OR_N=='T'
    // S(1:K,1:K) = MATMUL(Z(1:K,1:N),(W(1:N,1:K)))                 ! INTRINSIC, for T_OR_N=='N'
    // At this point S = U^H * A * U is the Rayleigh quotient.
    // If the residuals are requested, save scaled V_k into Z.
    // Rethat V_k or V_k^H is stored in W.
    if (WNTRES || WNTEX) {
      if (lsame(T_OR_N, 'N')) {
        zlacpy('A', N, K.value, W, LDW, Z, LDZ);
      } else {
        zlacpy('A', K.value, N, W, LDW, Z, LDZ);
      }
    }
  }

  //<5> Compute the Ritz values and (if requested) the
  //   right eigenvectors of the Rayleigh quotient.
  //
  zgeev('N', JOBZL, K.value, S, LDS, EIGS, W, LDW, W, LDW, ZWORK, LZWORK,
      RWORK(N + 1), INFO1);
  //
  // W(1:K,1:K) contains the eigenvectors of the Rayleigh
  // quotient.  See the description of Z.
  // Also, see the description of zgeev.
  if (INFO1.value > 0) {
    // zgeev failed to compute the eigenvalues and
    // eigenvectors of the Rayleigh quotient.
    INFO.value = 3;
    return;
  }

  // <6> Compute the eigenvectors (if requested) and,
  // the residuals (if requested).
  //
  if (WNTVEC || WNTEX) {
    if (WNTRES) {
      if (WNTREF) {
        // Here, if the refinement is requested, we have
        // A*U(:,1:K) already computed and stored in Z.
        // For the residuals, need Y = A * U(:,1;K) * W.
        zgemm('N', 'N', M, K.value, K.value, Complex.one, Z, LDZ, W, LDW,
            Complex.zero, Y, LDY);
        // Y(1:M,1:K) = Z(1:M,1:K) * W(1:K,1:K)        ! INTRINSIC
        // This frees Z; Y contains A * U(:,1:K) * W.
      } else {
        // Compute S = V_k * Sigma_k^(-1) * W, where
        // V_k * Sigma_k^(-1) (or its adjoint) is stored in Z
        zgemm(T_OR_N, 'N', N, K.value, K.value, Complex.one, Z, LDZ, W, LDW,
            Complex.zero, S, LDS);
        // Then, compute Z = Y * S =
        // = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
        // = A * U(:,1:K) * W(1:K,1:K)
        zgemm('N', 'N', M, K.value, N, Complex.one, Y, LDY, S, LDS,
            Complex.zero, Z, LDZ);
        // Save a copy of Z into Y and free Z for holding
        // the Ritz vectors.
        zlacpy('A', M, K.value, Z, LDZ, Y, LDY);
        if (WNTEX) zlacpy('A', M, K.value, Z, LDZ, B, LDB);
      }
    } else if (WNTEX) {
      // Compute S = V_k * Sigma_k^(-1) * W, where
      // V_k * Sigma_k^(-1) is stored in Z
      zgemm(T_OR_N, 'N', N, K.value, K.value, Complex.one, Z, LDZ, W, LDW,
          Complex.zero, S, LDS);
      // Then, compute Z = Y * S =
      // = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
      // = A * U(:,1:K) * W(1:K,1:K)
      zgemm('N', 'N', M, K.value, N, Complex.one, Y, LDY, S, LDS, Complex.zero,
          B, LDB);
      // The above replaces the following two calls
      // that were used in the developing-testing phase.
      // zgemm( 'N', 'N', M, K, N, Complex.one, Y, LDY, S,             !           LDS, Complex.zero, Z, LDZ)
      // Save a copy of Z into B and free Z for holding
      // the Ritz vectors.
      // zlacpy( 'A', M, K, Z, LDZ, B, LDB )
    }

    // Compute the Ritz vectors
    if (WNTVEC) {
      zgemm('N', 'N', M, K.value, K.value, Complex.one, X, LDX, W, LDW,
          Complex.zero, Z, LDZ);
      // Z(1:M,1:K) = MATMUL(X(1:M,1:K), W(1:K,1:K))         ! INTRINSIC
    }

    if (WNTRES) {
      for (final i in 1.through(K.value)) {
        zaxpy(M, -EIGS[i], Z(1, i).asArray(), 1, Y(1, i).asArray(), 1);
        // Y(1:M,i) = Y(1:M,i) - EIGS(i) * Z(1:M,i)            ! INTRINSIC
        RES[i] = dznrm2(M, Y(1, i).asArray(), 1);
      }
    }
  }

  if (WHTSVD == 4) {
    RWORK[N + 1] = XSCL1;
    RWORK[N + 2] = XSCL2;
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
