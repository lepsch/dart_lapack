// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zbbcsd.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlapmr.dart';
import 'package:dart_lapack/src/zlapmt.dart';
import 'package:dart_lapack/src/zunbdb.dart';
import 'package:dart_lapack/src/zunglq.dart';
import 'package:dart_lapack/src/zungqr.dart';

void zuncsd(
  final String JOBU1,
  final String JOBU2,
  final String JOBV1T,
  final String JOBV2T,
  final String TRANS,
  final String SIGNS,
  final int M,
  final int P,
  final int Q,
  final Matrix<Complex> X11_,
  final int LDX11,
  final Matrix<Complex> X12_,
  final int LDX12,
  final Matrix<Complex> X21_,
  final int LDX21,
  final Matrix<Complex> X22_,
  final int LDX22,
  final Array<double> THETA_,
  final Matrix<Complex> U1_,
  final int LDU1,
  final Matrix<Complex> U2_,
  final int LDU2,
  final Matrix<Complex> V1T_,
  final int LDV1T,
  final Matrix<Complex> V2T_,
  final int LDV2T,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final X11 = X11_.having(ld: LDX11);
  final X12 = X12_.having(ld: LDX12);
  final X21 = X21_.having(ld: LDX21);
  final X22 = X22_.having(ld: LDX22);
  final U1 = U1_.having(ld: LDU1);
  final U2 = U2_.having(ld: LDU2);
  final V1T = V1T_.having(ld: LDV1T);
  final V2T = V2T_.having(ld: LDV2T);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final THETA = THETA_.having();
  final RWORK = RWORK_.having();
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
      LBBCSDWORKMIN,
      LBBCSDWORKOPT,
      LORBDBWORK = 0,
      LORBDBWORKMIN,
      LORBDBWORKOPT,
      LORGLQWORK = 0,
      LORGLQWORKMIN,
      LORGLQWORKOPT,
      LORGQRWORK = 0,
      LORGQRWORKMIN,
      LORGQRWORKOPT,
      LWORKMIN = 0,
      LWORKOPT = 0,
      P1,
      Q1;
  bool COLMAJOR, DEFAULTSIGNS, LQUERY, WANTU1, WANTU2, WANTV1T, WANTV2T;
  int LRWORKMIN, LRWORKOPT;
  bool LRQUERY;
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
  LRQUERY = LRWORK == -1;
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
    zuncsd(
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
        RWORK,
        LRWORK,
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
    zuncsd(
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
        RWORK,
        LRWORK,
        IWORK,
        INFO);
    return;
  }

  // Compute workspace

  if (INFO.value == 0) {
    // Real workspace

    IPHI = 2;
    IB11D = IPHI + max(1, Q - 1);
    IB11E = IB11D + max(1, Q);
    IB12D = IB11E + max(1, Q - 1);
    IB12E = IB12D + max(1, Q);
    IB21D = IB12E + max(1, Q - 1);
    IB21E = IB21D + max(1, Q);
    IB22D = IB21E + max(1, Q - 1);
    IB22E = IB22D + max(1, Q);
    IBBCSD = IB22E + max(1, Q - 1);
    zbbcsd(
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
        THETA,
        THETA,
        THETA,
        THETA,
        THETA,
        THETA,
        THETA,
        THETA,
        RWORK,
        -1,
        CHILDINFO);
    LBBCSDWORKOPT = RWORK[1].toInt();
    LBBCSDWORKMIN = LBBCSDWORKOPT;
    LRWORKOPT = IBBCSD + LBBCSDWORKOPT - 1;
    LRWORKMIN = IBBCSD + LBBCSDWORKMIN - 1;
    RWORK[1] = LRWORKOPT.toDouble();

    // Complex workspace

    ITAUP1 = 2;
    ITAUP2 = ITAUP1 + max(1, P);
    ITAUQ1 = ITAUP2 + max(1, M - P);
    ITAUQ2 = ITAUQ1 + max(1, Q);
    IORGQR = ITAUQ2 + max(1, M - Q);
    zungqr(M - Q, M - Q, M - Q, U1, max(1, M - Q), U1.asArray(), WORK, -1,
        CHILDINFO);
    LORGQRWORKOPT = WORK[1].toInt();
    LORGQRWORKMIN = max(1, M - Q);
    IORGLQ = ITAUQ2 + max(1, M - Q);
    zunglq(M - Q, M - Q, M - Q, U1, max(1, M - Q), U1.asArray(), WORK, -1,
        CHILDINFO);
    LORGLQWORKOPT = WORK[1].toInt();
    LORGLQWORKMIN = max(1, M - Q);
    IORBDB = ITAUQ2 + max(1, M - Q);
    zunbdb(
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
        THETA,
        U1.asArray(),
        U2.asArray(),
        V1T.asArray(),
        V2T.asArray(),
        WORK,
        -1,
        CHILDINFO);
    LORBDBWORKOPT = WORK[1].toInt();
    LORBDBWORKMIN = LORBDBWORKOPT;
    LWORKOPT = max(IORGQR + LORGQRWORKOPT,
            max(IORGLQ + LORGLQWORKOPT, IORBDB + LORBDBWORKOPT)) -
        1;
    LWORKMIN = max(IORGQR + LORGQRWORKMIN,
            max(IORGLQ + LORGLQWORKMIN, IORBDB + LORBDBWORKMIN)) -
        1;
    WORK[1] = max(LWORKOPT, LWORKMIN).toComplex();

    if (LWORK < LWORKMIN && !(LQUERY || LRQUERY)) {
      INFO.value = -22;
    } else if (LRWORK < LRWORKMIN && !(LQUERY || LRQUERY)) {
      INFO.value = -24;
    } else {
      LORGQRWORK = LWORK - IORGQR + 1;
      LORGLQWORK = LWORK - IORGLQ + 1;
      LORBDBWORK = LWORK - IORBDB + 1;
      LBBCSDWORK = LRWORK - IBBCSD + 1;
    }
  }

  // Abort if any illegal arguments

  if (INFO.value != 0) {
    xerbla('ZUNCSD', -INFO.value);
    return;
  } else if (LQUERY || LRQUERY) {
    return;
  }

  // Transform to bidiagonal block form

  zunbdb(
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
      RWORK(IPHI),
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
      zlacpy('L', P, Q, X11, LDX11, U1, LDU1);
      zungqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQRWORK, INFO);
    }
    if (WANTU2 && M - P > 0) {
      zlacpy('L', M - P, Q, X21, LDX21, U2, LDU2);
      zungqr(M - P, M - P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQRWORK,
          INFO);
    }
    if (WANTV1T && Q > 0) {
      zlacpy('U', Q - 1, Q - 1, X11(1, 2), LDX11, V1T(2, 2), LDV1T);
      V1T[1][1] = Complex.one;
      for (J = 2; J <= Q; J++) {
        V1T[1][J] = Complex.zero;
        V1T[J][1] = Complex.zero;
      }
      zunglq(Q - 1, Q - 1, Q - 1, V1T(2, 2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ),
          LORGLQWORK, INFO);
    }
    if (WANTV2T && M - Q > 0) {
      zlacpy('U', P, M - Q, X12, LDX12, V2T, LDV2T);
      if (M - P > Q) {
        zlacpy('U', M - P - Q, M - P - Q, X22(Q + 1, P + 1), LDX22,
            V2T(P + 1, P + 1), LDV2T);
      }
      if (M > Q) {
        zunglq(M - Q, M - Q, M - Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGLQ),
            LORGLQWORK, INFO);
      }
    }
  } else {
    if (WANTU1 && P > 0) {
      zlacpy('U', Q, P, X11, LDX11, U1, LDU1);
      zunglq(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGLQ), LORGLQWORK, INFO);
    }
    if (WANTU2 && M - P > 0) {
      zlacpy('U', Q, M - P, X21, LDX21, U2, LDU2);
      zunglq(M - P, M - P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGLQ), LORGLQWORK,
          INFO);
    }
    if (WANTV1T && Q > 0) {
      zlacpy('L', Q - 1, Q - 1, X11(2, 1), LDX11, V1T(2, 2), LDV1T);
      V1T[1][1] = Complex.one;
      for (J = 2; J <= Q; J++) {
        V1T[1][J] = Complex.zero;
        V1T[J][1] = Complex.zero;
      }
      zungqr(Q - 1, Q - 1, Q - 1, V1T(2, 2), LDV1T, WORK(ITAUQ1), WORK(IORGQR),
          LORGQRWORK, INFO);
    }
    if (WANTV2T && M - Q > 0) {
      P1 = min(P + 1, M);
      Q1 = min(Q + 1, M);
      zlacpy('L', M - Q, P, X12, LDX12, V2T, LDV2T);
      if (M > P + Q) {
        zlacpy('L', M - P - Q, M - P - Q, X22(P1, Q1), LDX22, V2T(P + 1, P + 1),
            LDV2T);
      }
      zungqr(M - Q, M - Q, M - Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGQR),
          LORGQRWORK, INFO);
    }
  }

  // Compute the CSD of the matrix in bidiagonal-block form

  zbbcsd(
      JOBU1,
      JOBU2,
      JOBV1T,
      JOBV2T,
      TRANS,
      M,
      P,
      Q,
      THETA,
      RWORK(IPHI),
      U1,
      LDU1,
      U2,
      LDU2,
      V1T,
      LDV1T,
      V2T,
      LDV2T,
      RWORK(IB11D),
      RWORK(IB11E),
      RWORK(IB12D),
      RWORK(IB12E),
      RWORK(IB21D),
      RWORK(IB21E),
      RWORK(IB22D),
      RWORK(IB22E),
      RWORK(IBBCSD),
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
      zlapmt(false, M - P, M - P, U2, LDU2, IWORK);
    } else {
      zlapmr(false, M - P, M - P, U2, LDU2, IWORK);
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
      zlapmt(false, M - Q, M - Q, V2T, LDV2T, IWORK);
    } else {
      zlapmr(false, M - Q, M - Q, V2T, LDV2T, IWORK);
    }
  }
}
