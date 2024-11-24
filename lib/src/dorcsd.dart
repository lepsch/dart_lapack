// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dbbcsd.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlapmr.dart';
import 'package:dart_lapack/src/dlapmt.dart';
import 'package:dart_lapack/src/dorbdb.dart';
import 'package:dart_lapack/src/dorglq.dart';
import 'package:dart_lapack/src/dorgqr.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorcsd(
  final String JOBU1,
  final String JOBU2,
  final String JOBV1T,
  final String JOBV2T,
  final String TRANS,
  final String SIGNS,
  final int M,
  final int P,
  final int Q,
  final Matrix<double> X11_,
  final int LDX11,
  final Matrix<double> X12_,
  final int LDX12,
  final Matrix<double> X21_,
  final int LDX21,
  final Matrix<double> X22_,
  final int LDX22,
  final Array<double> THETA_,
  final Matrix<double> U1_,
  final int LDU1,
  final Matrix<double> U2_,
  final int LDU2,
  final Matrix<double> V1T_,
  final int LDV1T,
  final Matrix<double> V2T_,
  final int LDV2T,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X11 = X11_.having(ld: LDX11);
  final X12 = X12_.having(ld: LDX12);
  final X21 = X21_.having(ld: LDX21);
  final X22 = X22_.having(ld: LDX22);
  final THETA = THETA_.having();
  final U1 = U1_.having(ld: LDU1);
  final U2 = U2_.having(ld: LDU2);
  final V1T = V1T_.having(ld: LDV1T);
  final V2T = V2T_.having(ld: LDV2T);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  String TRANST, SIGNST;
  int I,
      IB11D = 0,
      IB11E = 0,
      IB12D = 0,
      IB12E = 0,
      IB21D = 0,
      IB21E = 0,
      IB22D = 0,
      IB22E = 0,
      IBBCSD = 0,
      IORBDB = 0,
      IORGLQ = 0,
      IORGQR = 0,
      IPHI = 0,
      ITAUP1 = 0,
      ITAUP2 = 0,
      ITAUQ1 = 0,
      ITAUQ2 = 0,
      J,
      LBBCSDWORK = 0,
      LBBCSDWORKMIN = 0,
      LBBCSDWORKOPT = 0,
      LORBDBWORK = 0,
      // LORBDBWORKMIN = 0,
      LORBDBWORKOPT = 0,
      LORGLQWORK = 0,
      LORGLQWORKMIN = 0,
      LORGLQWORKOPT = 0,
      LORGQRWORK = 0,
      LORGQRWORKMIN = 0,
      LORGQRWORKOPT = 0,
      LWORKMIN = 0,
      LWORKOPT = 0;
  bool COLMAJOR, DEFAULTSIGNS, LQUERY, WANTU1, WANTU2, WANTV1T, WANTV2T;
  final CHILDINFO = Box(0);

  // Test input arguments

  INFO.value = 0;
  WANTU1 = lsame(JOBU1, 'Y');
  WANTU2 = lsame(JOBU2, 'Y');
  WANTV1T = lsame(JOBV1T, 'Y');
  WANTV2T = lsame(JOBV2T, 'Y');
  COLMAJOR = !lsame(TRANS, 'T');
  DEFAULTSIGNS = !lsame(SIGNS, 'O');
  LQUERY = LWORK == -1;
  if (M < 0) {
    INFO.value = -7;
  } else if (P < 0 || P > M) {
    INFO.value = -8;
  } else if (Q < 0 || Q > M) {
    INFO.value = -9;
  } else if (COLMAJOR && LDX11 < max(1, P)) {
    INFO.value = -11;
  } else if (!COLMAJOR && LDX11 < max(1, Q)) {
    INFO.value = -11;
  } else if (COLMAJOR && LDX12 < max(1, P)) {
    INFO.value = -13;
  } else if (!COLMAJOR && LDX12 < max(1, M - Q)) {
    INFO.value = -13;
  } else if (COLMAJOR && LDX21 < max(1, M - P)) {
    INFO.value = -15;
  } else if (!COLMAJOR && LDX21 < max(1, Q)) {
    INFO.value = -15;
  } else if (COLMAJOR && LDX22 < max(1, M - P)) {
    INFO.value = -17;
  } else if (!COLMAJOR && LDX22 < max(1, M - Q)) {
    INFO.value = -17;
  } else if (WANTU1 && LDU1 < P) {
    INFO.value = -20;
  } else if (WANTU2 && LDU2 < M - P) {
    INFO.value = -22;
  } else if (WANTV1T && LDV1T < Q) {
    INFO.value = -24;
  } else if (WANTV2T && LDV2T < M - Q) {
    INFO.value = -26;
  }

  // Work with transpose if convenient

  if (INFO.value == 0 && min(P, M - P) < min(Q, M - Q)) {
    if (COLMAJOR) {
      TRANST = 'T';
    } else {
      TRANST = 'N';
    }
    if (DEFAULTSIGNS) {
      SIGNST = 'O';
    } else {
      SIGNST = 'D';
    }
    dorcsd(
        JOBV1T,
        JOBV2T,
        JOBU1,
        JOBU2,
        TRANST,
        SIGNST,
        M,
        Q,
        P,
        X11,
        LDX11,
        X21,
        LDX21,
        X12,
        LDX12,
        X22,
        LDX22,
        THETA,
        V1T,
        LDV1T,
        V2T,
        LDV2T,
        U1,
        LDU1,
        U2,
        LDU2,
        WORK,
        LWORK,
        IWORK,
        INFO);
    return;
  }

  // Work with permutation [ 0 I; I 0 ] * X * [ 0 I; I 0 ] if
  // convenient

  if (INFO.value == 0 && M - Q < Q) {
    if (DEFAULTSIGNS) {
      SIGNST = 'O';
    } else {
      SIGNST = 'D';
    }
    dorcsd(
        JOBU2,
        JOBU1,
        JOBV2T,
        JOBV1T,
        TRANS,
        SIGNST,
        M,
        M - P,
        M - Q,
        X22,
        LDX22,
        X21,
        LDX21,
        X12,
        LDX12,
        X11,
        LDX11,
        THETA,
        U2,
        LDU2,
        U1,
        LDU1,
        V2T,
        LDV2T,
        V1T,
        LDV1T,
        WORK,
        LWORK,
        IWORK,
        INFO);
    return;
  }

  // Compute workspace

  if (INFO.value == 0) {
    IPHI = 2;
    ITAUP1 = IPHI + max(1, Q - 1);
    ITAUP2 = ITAUP1 + max(1, P);
    ITAUQ1 = ITAUP2 + max(1, M - P);
    ITAUQ2 = ITAUQ1 + max(1, Q);
    IORGQR = ITAUQ2 + max(1, M - Q);
    dorgqr(M - Q, M - Q, M - Q, U1, max(1, M - Q), U1.asArray(), WORK, -1,
        CHILDINFO);
    LORGQRWORKOPT = WORK[1].toInt();
    LORGQRWORKMIN = max(1, M - Q);
    IORGLQ = ITAUQ2 + max(1, M - Q);
    dorglq(M - Q, M - Q, M - Q, U1, max(1, M - Q), U1.asArray(), WORK, -1,
        CHILDINFO);
    LORGLQWORKOPT = WORK[1].toInt();
    LORGLQWORKMIN = max(1, M - Q);
    IORBDB = ITAUQ2 + max(1, M - Q);
    dorbdb(
        TRANS,
        SIGNS,
        M,
        P,
        Q,
        X11,
        LDX11,
        X12,
        LDX12,
        X21,
        LDX21,
        X22,
        LDX22,
        THETA,
        V1T.asArray(),
        U1.asArray(),
        U2.asArray(),
        V1T.asArray(),
        V2T.asArray(),
        WORK,
        -1,
        CHILDINFO);
    LORBDBWORKOPT = WORK[1].toInt();
    // LORBDBWORKMIN = LORBDBWORKOPT;
    IB11D = ITAUQ2 + max(1, M - Q);
    IB11E = IB11D + max(1, Q);
    IB12D = IB11E + max(1, Q - 1);
    IB12E = IB12D + max(1, Q);
    IB21D = IB12E + max(1, Q - 1);
    IB21E = IB21D + max(1, Q);
    IB22D = IB21E + max(1, Q - 1);
    IB22E = IB22D + max(1, Q);
    IBBCSD = IB22E + max(1, Q - 1);
    dbbcsd(
        JOBU1,
        JOBU2,
        JOBV1T,
        JOBV2T,
        TRANS,
        M,
        P,
        Q,
        THETA,
        THETA,
        U1,
        LDU1,
        U2,
        LDU2,
        V1T,
        LDV1T,
        V2T,
        LDV2T,
        U1.asArray(),
        U1.asArray(),
        U1.asArray(),
        U1.asArray(),
        U1.asArray(),
        U1.asArray(),
        U1.asArray(),
        U1.asArray(),
        WORK,
        -1,
        CHILDINFO);
    LBBCSDWORKOPT = WORK[1].toInt();
    LBBCSDWORKMIN = LBBCSDWORKOPT;
    LWORKOPT = max(
          max(IORGQR + LORGQRWORKOPT, IORGLQ + LORGLQWORKOPT),
          max(IORBDB + LORBDBWORKOPT, IBBCSD + LBBCSDWORKOPT),
        ) -
        1;
    LWORKMIN = max(
          max(IORGQR + LORGQRWORKMIN, IORGLQ + LORGLQWORKMIN),
          max(IORBDB + LORBDBWORKOPT, IBBCSD + LBBCSDWORKMIN),
        ) -
        1;
    WORK[1] = max(LWORKOPT, LWORKMIN).toDouble();

    if (LWORK < LWORKMIN && !LQUERY) {
      INFO.value = -22;
    } else {
      LORGQRWORK = LWORK - IORGQR + 1;
      LORGLQWORK = LWORK - IORGLQ + 1;
      LORBDBWORK = LWORK - IORBDB + 1;
      LBBCSDWORK = LWORK - IBBCSD + 1;
    }
  }

  // Abort if any illegal arguments

  if (INFO.value != 0) {
    xerbla('DORCSD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Transform to bidiagonal block form

  dorbdb(
      TRANS,
      SIGNS,
      M,
      P,
      Q,
      X11,
      LDX11,
      X12,
      LDX12,
      X21,
      LDX21,
      X22,
      LDX22,
      THETA,
      WORK(IPHI),
      WORK(ITAUP1),
      WORK(ITAUP2),
      WORK(ITAUQ1),
      WORK(ITAUQ2),
      WORK(IORBDB),
      LORBDBWORK,
      CHILDINFO);

  // Accumulate Householder reflectors

  if (COLMAJOR) {
    if (WANTU1 && P > 0) {
      dlacpy('L', P, Q, X11, LDX11, U1, LDU1);
      dorgqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQRWORK, INFO);
    }
    if (WANTU2 && M - P > 0) {
      dlacpy('L', M - P, Q, X21, LDX21, U2, LDU2);
      dorgqr(M - P, M - P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQRWORK,
          INFO);
    }
    if (WANTV1T && Q > 0) {
      dlacpy('U', Q - 1, Q - 1, X11(1, 2), LDX11, V1T(2, 2), LDV1T);
      V1T[1][1] = ONE;
      for (J = 2; J <= Q; J++) {
        V1T[1][J] = ZERO;
        V1T[J][1] = ZERO;
      }
      dorglq(Q - 1, Q - 1, Q - 1, V1T(2, 2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ),
          LORGLQWORK, INFO);
    }
    if (WANTV2T && M - Q > 0) {
      dlacpy('U', P, M - Q, X12, LDX12, V2T, LDV2T);
      if (M - P > Q) {
        dlacpy('U', M - P - Q, M - P - Q, X22(Q + 1, P + 1), LDX22,
            V2T(P + 1, P + 1), LDV2T);
      }
      if (M > Q) {
        dorglq(M - Q, M - Q, M - Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGLQ),
            LORGLQWORK, INFO);
      }
    }
  } else {
    if (WANTU1 && P > 0) {
      dlacpy('U', Q, P, X11, LDX11, U1, LDU1);
      dorglq(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGLQ), LORGLQWORK, INFO);
    }
    if (WANTU2 && M - P > 0) {
      dlacpy('U', Q, M - P, X21, LDX21, U2, LDU2);
      dorglq(M - P, M - P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGLQ), LORGLQWORK,
          INFO);
    }
    if (WANTV1T && Q > 0) {
      dlacpy('L', Q - 1, Q - 1, X11(2, 1), LDX11, V1T(2, 2), LDV1T);
      V1T[1][1] = ONE;
      for (J = 2; J <= Q; J++) {
        V1T[1][J] = ZERO;
        V1T[J][1] = ZERO;
      }
      dorgqr(Q - 1, Q - 1, Q - 1, V1T(2, 2), LDV1T, WORK(ITAUQ1), WORK(IORGQR),
          LORGQRWORK, INFO);
    }
    if (WANTV2T && M - Q > 0) {
      dlacpy('L', M - Q, P, X12, LDX12, V2T, LDV2T);
      dlacpy('L', M - P - Q, M - P - Q, X22(P + 1, Q + 1), LDX22,
          V2T(P + 1, P + 1), LDV2T);
      dorgqr(M - Q, M - Q, M - Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGQR),
          LORGQRWORK, INFO);
    }
  }

  // Compute the CSD of the matrix in bidiagonal-block form

  dbbcsd(
      JOBU1,
      JOBU2,
      JOBV1T,
      JOBV2T,
      TRANS,
      M,
      P,
      Q,
      THETA,
      WORK(IPHI),
      U1,
      LDU1,
      U2,
      LDU2,
      V1T,
      LDV1T,
      V2T,
      LDV2T,
      WORK(IB11D),
      WORK(IB11E),
      WORK(IB12D),
      WORK(IB12E),
      WORK(IB21D),
      WORK(IB21E),
      WORK(IB22D),
      WORK(IB22E),
      WORK(IBBCSD),
      LBBCSDWORK,
      INFO);

  // Permute rows and columns to place identity submatrices in top-
  // left corner of (1,1)-block and/or bottom-right corner of (1,2)-
  // block and/or bottom-right corner of (2,1)-block and/or top-left
  // corner of (2,2)-block

  if (Q > 0 && WANTU2) {
    for (I = 1; I <= Q; I++) {
      IWORK[I] = M - P - Q + I;
    }
    for (I = Q + 1; I <= M - P; I++) {
      IWORK[I] = I - Q;
    }
    if (COLMAJOR) {
      dlapmt(false, M - P, M - P, U2, LDU2, IWORK);
    } else {
      dlapmr(false, M - P, M - P, U2, LDU2, IWORK);
    }
  }
  if (M > 0 && WANTV2T) {
    for (I = 1; I <= P; I++) {
      IWORK[I] = M - P - Q + I;
    }
    for (I = P + 1; I <= M - Q; I++) {
      IWORK[I] = I - P;
    }
    if (!COLMAJOR) {
      dlapmt(false, M - Q, M - Q, V2T, LDV2T, IWORK);
    } else {
      dlapmr(false, M - Q, M - Q, V2T, LDV2T, IWORK);
    }
  }
}
