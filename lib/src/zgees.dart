// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgebak.dart';
import 'package:lapack/src/zgebal.dart';
import 'package:lapack/src/zgehrd.dart';
import 'package:lapack/src/zhseqr.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/ztrsen.dart';
import 'package:lapack/src/zunghr.dart';

void zgees(
  final String JOBVS,
  final String SORT,
  final bool Function(Complex) SELECT,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> SDIM,
  final Array<Complex> W_,
  final Matrix<Complex> VS_,
  final int LDVS,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final W = W_.having();
  final VS = VS_.having(ld: LDVS);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final BWORK = BWORK_.having();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, SCALEA, WANTST, WANTVS;
  int HSWORK, I, IBAL, ITAU, IWRK, MAXWRK = 0, MINWRK;
  double ANRM, BIGNUM, CSCALE = 0, EPS, SMLNUM;
  final DUM = Array<double>(1);
  final IERR = Box(0),
      IEVAL = Box(0),
      IHI = Box(0),
      ILO = Box(0),
      ICOND = Box(0);
  final S = Box(0.0), SEP = Box(0.0);

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  WANTVS = lsame(JOBVS, 'V');
  WANTST = lsame(SORT, 'S');
  if (!WANTVS && !lsame(JOBVS, 'N')) {
    INFO.value = -1;
  } else if (!WANTST && !lsame(SORT, 'N')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDVS < 1 || (WANTVS && LDVS < N)) {
    INFO.value = -10;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   CWorkspace refers to complex workspace, and RWorkspace to real
  //   workspace. NB refers to the optimal block size for the
  //   immediately following subroutine, as returned by ILAENV.
  //   HSWORK refers to the workspace preferred by ZHSEQR, as
  //   calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
  //   the worst case.)

  if (INFO.value == 0) {
    if (N == 0) {
      MINWRK = 1;
      MAXWRK = 1;
    } else {
      MAXWRK = N + N * ilaenv(1, 'ZGEHRD', ' ', N, 1, N, 0);
      MINWRK = 2 * N;

      zhseqr('S', JOBVS, N, 1, N, A, LDA, W, VS, LDVS, WORK, -1, IEVAL);
      HSWORK = WORK[1].toInt();

      if (!WANTVS) {
        MAXWRK = max(MAXWRK, HSWORK);
      } else {
        MAXWRK =
            max(MAXWRK, N + (N - 1) * ilaenv(1, 'ZUNGHR', ' ', N, 1, N, -1));
        MAXWRK = max(MAXWRK, HSWORK);
      }
    }
    WORK[1] = MAXWRK.toComplex();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -12;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZGEES', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    SDIM.value = 0;
    return;
  }

  // Get machine constants

  EPS = dlamch('P');
  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;
  SMLNUM = sqrt(SMLNUM) / EPS;
  BIGNUM = ONE / SMLNUM;

  // Scale A if max element outside range [SMLNUM,BIGNUM]

  ANRM = zlange('M', N, N, A, LDA, DUM);
  SCALEA = false;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    SCALEA = true;
    CSCALE = SMLNUM;
  } else if (ANRM > BIGNUM) {
    SCALEA = true;
    CSCALE = BIGNUM;
  }
  if (SCALEA) zlascl('G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR);

  // Permute the matrix to make it more nearly triangular
  // (CWorkspace: none)
  // (RWorkspace: need N)

  IBAL = 1;
  zgebal('P', N, A, LDA, ILO, IHI, RWORK(IBAL), IERR);

  // Reduce to upper Hessenberg form
  // (CWorkspace: need 2*N, prefer N+N*NB)
  // (RWorkspace: none)

  ITAU = 1;
  IWRK = N + ITAU;
  zgehrd(N, ILO.value, IHI.value, A, LDA, WORK(ITAU), WORK(IWRK),
      LWORK - IWRK + 1, IERR);

  if (WANTVS) {
    // Copy Householder vectors to VS

    zlacpy('L', N, N, A, LDA, VS, LDVS);

    // Generate unitary matrix in VS
    // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
    // (RWorkspace: none)

    zunghr(N, ILO.value, IHI.value, VS, LDVS, WORK(ITAU), WORK(IWRK),
        LWORK - IWRK + 1, IERR);
  }

  SDIM.value = 0;

  // Perform QR iteration, accumulating Schur vectors in VS if desired
  // (CWorkspace: need 1, prefer HSWORK (see comments) )
  // (RWorkspace: none)

  IWRK = ITAU;
  zhseqr('S', JOBVS, N, ILO.value, IHI.value, A, LDA, W, VS, LDVS, WORK(IWRK),
      LWORK - IWRK + 1, IEVAL);
  if (IEVAL.value > 0) INFO.value = IEVAL.value;

  // Sort eigenvalues if desired

  if (WANTST && INFO.value == 0) {
    if (SCALEA) zlascl('G', 0, 0, CSCALE, ANRM, N, 1, W.asMatrix(N), N, IERR);
    for (I = 1; I <= N; I++) {
      BWORK[I] = SELECT(W[I]);
    }

    // Reorder eigenvalues and transform Schur vectors
    // (CWorkspace: none)
    // (RWorkspace: none)

    ztrsen('N', JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM, S, SEP, WORK(IWRK),
        LWORK - IWRK + 1, ICOND);
  }

  if (WANTVS) {
    // Undo balancing
    // (CWorkspace: none)
    // (RWorkspace: need N)

    zgebak('P', 'R', N, ILO.value, IHI.value, RWORK(IBAL), N, VS, LDVS, IERR);
  }

  if (SCALEA) {
    // Undo scaling for the Schur form of A

    zlascl('U', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR);
    zcopy(N, A.asArray(), LDA + 1, W, 1);
  }

  WORK[1] = MAXWRK.toComplex();
}
