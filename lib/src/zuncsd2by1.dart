// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zbbcsd.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlapmr.dart';
import 'package:dart_lapack/src/zlapmt.dart';
import 'package:dart_lapack/src/zunbdb1.dart';
import 'package:dart_lapack/src/zunbdb2.dart';
import 'package:dart_lapack/src/zunbdb3.dart';
import 'package:dart_lapack/src/zunbdb4.dart';
import 'package:dart_lapack/src/zunglq.dart';
import 'package:dart_lapack/src/zungqr.dart';

void zuncsd2by1(
  final String JOBU1,
  final String JOBU2,
  final String JOBV1T,
  final int M,
  final int P,
  final int Q,
  final Matrix<Complex> X11_,
  final int LDX11,
  final Matrix<Complex> X21_,
  final int LDX21,
  final Array<double> THETA_,
  final Matrix<Complex> U1_,
  final int LDU1,
  final Matrix<Complex> U2_,
  final int LDU2,
  final Matrix<Complex> V1T_,
  final int LDV1T,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final X11 = X11_.having(ld: LDX11);
  final X21 = X21_.having(ld: LDX21);
  final U1 = U1_.having(ld: LDU1);
  final U2 = U2_.having(ld: LDU2);
  final V1T = V1T_.having(ld: LDV1T);
  final THETA = THETA_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  int LRWORKMIN, LRWORKOPT;

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
      J,
      LBBCSD = 0,
      LORBDB = 0,
      LORGLQ = 0,
      LORGLQMIN,
      LORGLQOPT,
      LORGQR = 0,
      LORGQRMIN,
      LORGQROPT,
      LWORKMIN,
      LWORKOPT,
      R;
  bool LQUERY, WANTU1, WANTU2, WANTV1T;
  final DUM = Array<double>(1);
  final CDUM = Matrix<Complex>(1, 1);
  final CHILDINFO = Box(0);

  // Test input arguments

  INFO.value = 0;
  WANTU1 = lsame(JOBU1, 'Y');
  WANTU2 = lsame(JOBU2, 'Y');
  WANTV1T = lsame(JOBV1T, 'Y');
  LQUERY = (LWORK == -1) || (LRWORK == -1);

  if (M < 0) {
    INFO.value = -4;
  } else if (P < 0 || P > M) {
    INFO.value = -5;
  } else if (Q < 0 || Q > M) {
    INFO.value = -6;
  } else if (LDX11 < max(1, P)) {
    INFO.value = -8;
  } else if (LDX21 < max(1, M - P)) {
    INFO.value = -10;
  } else if (WANTU1 && LDU1 < max(1, P)) {
    INFO.value = -13;
  } else if (WANTU2 && LDU2 < max(1, M - P)) {
    INFO.value = -15;
  } else if (WANTV1T && LDV1T < max(1, Q)) {
    INFO.value = -17;
  }

  R = min(min(P, M - P), min(Q, M - Q));

  // Compute workspace

  //   WORK layout:
  // |-----------------------------------------|
  // | LWORKOPT (1)                            |
  // |-----------------------------------------|
  // | TAUP1 (max(1,P))                        |
  // | TAUP2 (max(1,M-P))                      |
  // | TAUQ1 (max(1,Q))                        |
  // |-----------------------------------------|
  // | ZUNBDB WORK | ZUNGQR WORK | ZUNGLQ WORK |
  // |             |             |             |
  // |             |             |             |
  // |             |             |             |
  // |             |             |             |
  // |-----------------------------------------|
  //   RWORK layout:
  // |------------------|
  // | LRWORKOPT (1)    |
  // |------------------|
  // | PHI (max(1,R-1)) |
  // |------------------|
  // | B11D (R)         |
  // | B11E (R-1)       |
  // | B12D (R)         |
  // | B12E (R-1)       |
  // | B21D (R)         |
  // | B21E (R-1)       |
  // | B22D (R)         |
  // | B22E (R-1)       |
  // | ZBBCSD RWORK     |
  // |------------------|

  if (INFO.value == 0) {
    IPHI = 2;
    IB11D = IPHI + max(1, R - 1);
    IB11E = IB11D + max(1, R);
    IB12D = IB11E + max(1, R - 1);
    IB12E = IB12D + max(1, R);
    IB21D = IB12E + max(1, R - 1);
    IB21E = IB21D + max(1, R);
    IB22D = IB21E + max(1, R - 1);
    IB22E = IB22D + max(1, R);
    IBBCSD = IB22E + max(1, R - 1);
    ITAUP1 = 2;
    ITAUP2 = ITAUP1 + max(1, P);
    ITAUQ1 = ITAUP2 + max(1, M - P);
    IORBDB = ITAUQ1 + max(1, Q);
    IORGQR = ITAUQ1 + max(1, Q);
    IORGLQ = ITAUQ1 + max(1, Q);
    LORGQRMIN = 1;
    LORGQROPT = 1;
    LORGLQMIN = 1;
    LORGLQOPT = 1;
    if (R == Q) {
      zunbdb1(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM.asArray(),
          CDUM.asArray(), CDUM.asArray(), WORK, -1, CHILDINFO);
      LORBDB = WORK[1].toInt();
      if (WANTU1 && P > 0) {
        zungqr(P, P, Q, U1, LDU1, CDUM.asArray(), WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTU2 && M - P > 0) {
        zungqr(
            M - P, M - P, Q, U2, LDU2, CDUM.asArray(), WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, M - P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTV1T && Q > 0) {
        zunglq(Q - 1, Q - 1, Q - 1, V1T, LDV1T, CDUM.asArray(), WORK(1), -1,
            CHILDINFO);
        LORGLQMIN = max(LORGLQMIN, Q - 1);
        LORGLQOPT = max(LORGLQOPT, WORK[1].toInt());
      }
      zbbcsd(
          JOBU1,
          JOBU2,
          JOBV1T,
          'N',
          'N',
          M,
          P,
          Q,
          THETA,
          DUM,
          U1,
          LDU1,
          U2,
          LDU2,
          V1T,
          LDV1T,
          CDUM,
          1,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          RWORK(1),
          -1,
          CHILDINFO);
      LBBCSD = RWORK[1].toInt();
    } else if (R == P) {
      zunbdb2(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM.asArray(),
          CDUM.asArray(), CDUM.asArray(), WORK(1), -1, CHILDINFO);
      LORBDB = WORK[1].toInt();
      if (WANTU1 && P > 0) {
        zungqr(P - 1, P - 1, P - 1, U1(2, 2), LDU1, CDUM.asArray(), WORK(1), -1,
            CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, P - 1);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTU2 && M - P > 0) {
        zungqr(
            M - P, M - P, Q, U2, LDU2, CDUM.asArray(), WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, M - P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTV1T && Q > 0) {
        zunglq(Q, Q, R, V1T, LDV1T, CDUM.asArray(), WORK(1), -1, CHILDINFO);
        LORGLQMIN = max(LORGLQMIN, Q);
        LORGLQOPT = max(LORGLQOPT, WORK[1].toInt());
      }
      zbbcsd(
          JOBV1T,
          'N',
          JOBU1,
          JOBU2,
          'T',
          M,
          Q,
          P,
          THETA,
          DUM,
          V1T,
          LDV1T,
          CDUM,
          1,
          U1,
          LDU1,
          U2,
          LDU2,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          RWORK(1),
          -1,
          CHILDINFO);
      LBBCSD = RWORK[1].toInt();
    } else if (R == M - P) {
      zunbdb3(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM.asArray(),
          CDUM.asArray(), CDUM.asArray(), WORK(1), -1, CHILDINFO);
      LORBDB = WORK[1].toInt();
      if (WANTU1 && P > 0) {
        zungqr(P, P, Q, U1, LDU1, CDUM.asArray(), WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTU2 && M - P > 0) {
        zungqr(M - P - 1, M - P - 1, M - P - 1, U2(2, 2), LDU2, CDUM.asArray(),
            WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, M - P - 1);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTV1T && Q > 0) {
        zunglq(Q, Q, R, V1T, LDV1T, CDUM.asArray(), WORK(1), -1, CHILDINFO);
        LORGLQMIN = max(LORGLQMIN, Q);
        LORGLQOPT = max(LORGLQOPT, WORK[1].toInt());
      }
      zbbcsd(
          'N',
          JOBV1T,
          JOBU2,
          JOBU1,
          'T',
          M,
          M - Q,
          M - P,
          THETA,
          DUM,
          CDUM,
          1,
          V1T,
          LDV1T,
          U2,
          LDU2,
          U1,
          LDU1,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          RWORK(1),
          -1,
          CHILDINFO);
      LBBCSD = RWORK[1].toInt();
    } else {
      zunbdb4(
          M,
          P,
          Q,
          X11,
          LDX11,
          X21,
          LDX21,
          THETA,
          DUM,
          CDUM.asArray(),
          CDUM.asArray(),
          CDUM.asArray(),
          CDUM.asArray(),
          WORK(1),
          -1,
          CHILDINFO);
      LORBDB = M + WORK[1].toInt();
      if (WANTU1 && P > 0) {
        zungqr(P, P, M - Q, U1, LDU1, CDUM.asArray(), WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTU2 && M - P > 0) {
        zungqr(M - P, M - P, M - Q, U2, LDU2, CDUM.asArray(), WORK(1), -1,
            CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, M - P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTV1T && Q > 0) {
        zunglq(Q, Q, Q, V1T, LDV1T, CDUM.asArray(), WORK(1), -1, CHILDINFO);
        LORGLQMIN = max(LORGLQMIN, Q);
        LORGLQOPT = max(LORGLQOPT, WORK[1].toInt());
      }
      zbbcsd(
          JOBU2,
          JOBU1,
          'N',
          JOBV1T,
          'N',
          M,
          M - P,
          M - Q,
          THETA,
          DUM,
          U2,
          LDU2,
          U1,
          LDU1,
          CDUM,
          1,
          V1T,
          LDV1T,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          DUM,
          RWORK(1),
          -1,
          CHILDINFO);
      LBBCSD = RWORK[1].toInt();
    }
    LRWORKMIN = IBBCSD + LBBCSD - 1;
    LRWORKOPT = LRWORKMIN;
    RWORK[1] = LRWORKOPT.toDouble();
    LWORKMIN = max(IORBDB + LORBDB - 1,
        max(IORGQR + LORGQRMIN - 1, IORGLQ + LORGLQMIN - 1));
    LWORKOPT = max(IORBDB + LORBDB - 1,
        max(IORGQR + LORGQROPT - 1, IORGLQ + LORGLQOPT - 1));
    WORK[1] = LWORKOPT.toComplex();
    if (LWORK < LWORKMIN && !LQUERY) {
      INFO.value = -19;
    }
    if (LRWORK < LRWORKMIN && !LQUERY) {
      INFO.value = -21;
    }
  }
  if (INFO.value != 0) {
    xerbla('ZUNCSD2BY1', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }
  LORGQR = LWORK - IORGQR + 1;
  LORGLQ = LWORK - IORGLQ + 1;

  // Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q,
  // in which R = min(P,M-P,Q,M-Q)

  if (R == Q) {
    // Case 1: R = Q

    // Simultaneously bidiagonalize X11 and X21

    zunbdb1(M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1),
        WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO);

    // Accumulate Householder reflectors

    if (WANTU1 && P > 0) {
      zlacpy('L', P, Q, X11, LDX11, U1, LDU1);
      zungqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO);
    }
    if (WANTU2 && M - P > 0) {
      zlacpy('L', M - P, Q, X21, LDX21, U2, LDU2);
      zungqr(M - P, M - P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR,
          CHILDINFO);
    }
    if (WANTV1T && Q > 0) {
      V1T[1][1] = Complex.one;
      for (J = 2; J <= Q; J++) {
        V1T[1][J] = Complex.zero;
        V1T[J][1] = Complex.zero;
      }
      zlacpy('U', Q - 1, Q - 1, X21(1, 2), LDX21, V1T(2, 2), LDV1T);
      zunglq(Q - 1, Q - 1, Q - 1, V1T(2, 2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ),
          LORGLQ, CHILDINFO);
    }

    // Simultaneously diagonalize X11 and X21.

    zbbcsd(
        JOBU1,
        JOBU2,
        JOBV1T,
        'N',
        'N',
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
        CDUM,
        1,
        RWORK(IB11D),
        RWORK(IB11E),
        RWORK(IB12D),
        RWORK(IB12E),
        RWORK(IB21D),
        RWORK(IB21E),
        RWORK(IB22D),
        RWORK(IB22E),
        RWORK(IBBCSD),
        LRWORK - IBBCSD + 1,
        CHILDINFO);

    // Permute rows and columns to place zero submatrices in
    // preferred positions

    if (Q > 0 && WANTU2) {
      for (I = 1; I <= Q; I++) {
        IWORK[I] = M - P - Q + I;
      }
      for (I = Q + 1; I <= M - P; I++) {
        IWORK[I] = I - Q;
      }
      zlapmt(false, M - P, M - P, U2, LDU2, IWORK);
    }
  } else if (R == P) {
    // Case 2: R = P

    // Simultaneously bidiagonalize X11 and X21

    zunbdb2(M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1),
        WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO);

    // Accumulate Householder reflectors

    if (WANTU1 && P > 0) {
      U1[1][1] = Complex.one;
      for (J = 2; J <= P; J++) {
        U1[1][J] = Complex.zero;
        U1[J][1] = Complex.zero;
      }
      zlacpy('L', P - 1, P - 1, X11(2, 1), LDX11, U1(2, 2), LDU1);
      zungqr(P - 1, P - 1, P - 1, U1(2, 2), LDU1, WORK(ITAUP1), WORK(IORGQR),
          LORGQR, CHILDINFO);
    }
    if (WANTU2 && M - P > 0) {
      zlacpy('L', M - P, Q, X21, LDX21, U2, LDU2);
      zungqr(M - P, M - P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR,
          CHILDINFO);
    }
    if (WANTV1T && Q > 0) {
      zlacpy('U', P, Q, X11, LDX11, V1T, LDV1T);
      zunglq(
          Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO);
    }

    // Simultaneously diagonalize X11 and X21.

    zbbcsd(
        JOBV1T,
        'N',
        JOBU1,
        JOBU2,
        'T',
        M,
        Q,
        P,
        THETA,
        RWORK(IPHI),
        V1T,
        LDV1T,
        CDUM,
        1,
        U1,
        LDU1,
        U2,
        LDU2,
        RWORK(IB11D),
        RWORK(IB11E),
        RWORK(IB12D),
        RWORK(IB12E),
        RWORK(IB21D),
        RWORK(IB21E),
        RWORK(IB22D),
        RWORK(IB22E),
        RWORK(IBBCSD),
        LBBCSD,
        CHILDINFO);

    // Permute rows and columns to place identity submatrices in
    // preferred positions

    if (Q > 0 && WANTU2) {
      for (I = 1; I <= Q; I++) {
        IWORK[I] = M - P - Q + I;
      }
      for (I = Q + 1; I <= M - P; I++) {
        IWORK[I] = I - Q;
      }
      zlapmt(false, M - P, M - P, U2, LDU2, IWORK);
    }
  } else if (R == M - P) {
    // Case 3: R = M-P

    // Simultaneously bidiagonalize X11 and X21

    zunbdb3(M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1),
        WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO);

    // Accumulate Householder reflectors

    if (WANTU1 && P > 0) {
      zlacpy('L', P, Q, X11, LDX11, U1, LDU1);
      zungqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO);
    }
    if (WANTU2 && M - P > 0) {
      U2[1][1] = Complex.one;
      for (J = 2; J <= M - P; J++) {
        U2[1][J] = Complex.zero;
        U2[J][1] = Complex.zero;
      }
      zlacpy('L', M - P - 1, M - P - 1, X21(2, 1), LDX21, U2(2, 2), LDU2);
      zungqr(M - P - 1, M - P - 1, M - P - 1, U2(2, 2), LDU2, WORK(ITAUP2),
          WORK(IORGQR), LORGQR, CHILDINFO);
    }
    if (WANTV1T && Q > 0) {
      zlacpy('U', M - P, Q, X21, LDX21, V1T, LDV1T);
      zunglq(
          Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO);
    }

    // Simultaneously diagonalize X11 and X21.

    zbbcsd(
        'N',
        JOBV1T,
        JOBU2,
        JOBU1,
        'T',
        M,
        M - Q,
        M - P,
        THETA,
        RWORK(IPHI),
        CDUM,
        1,
        V1T,
        LDV1T,
        U2,
        LDU2,
        U1,
        LDU1,
        RWORK(IB11D),
        RWORK(IB11E),
        RWORK(IB12D),
        RWORK(IB12E),
        RWORK(IB21D),
        RWORK(IB21E),
        RWORK(IB22D),
        RWORK(IB22E),
        RWORK(IBBCSD),
        LBBCSD,
        CHILDINFO);

    // Permute rows and columns to place identity submatrices in
    // preferred positions

    if (Q > R) {
      for (I = 1; I <= R; I++) {
        IWORK[I] = Q - R + I;
      }
      for (I = R + 1; I <= Q; I++) {
        IWORK[I] = I - R;
      }
      if (WANTU1) {
        zlapmt(false, P, Q, U1, LDU1, IWORK);
      }
      if (WANTV1T) {
        zlapmr(false, Q, Q, V1T, LDV1T, IWORK);
      }
    }
  } else {
    // Case 4: R = M-Q

    // Simultaneously bidiagonalize X11 and X21

    zunbdb4(
        M,
        P,
        Q,
        X11,
        LDX11,
        X21,
        LDX21,
        THETA,
        RWORK(IPHI),
        WORK(ITAUP1),
        WORK(ITAUP2),
        WORK(ITAUQ1),
        WORK(IORBDB),
        WORK(IORBDB + M),
        LORBDB - M,
        CHILDINFO);

    // Accumulate Householder reflectors

    if (WANTU2 && M - P > 0) {
      zcopy(M - P, WORK(IORBDB + P), 1, U2.asArray(), 1);
    }
    if (WANTU1 && P > 0) {
      zcopy(P, WORK(IORBDB), 1, U1.asArray(), 1);
      for (J = 2; J <= P; J++) {
        U1[1][J] = Complex.zero;
      }
      zlacpy('L', P - 1, M - Q - 1, X11(2, 1), LDX11, U1(2, 2), LDU1);
      zungqr(
          P, P, M - Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO);
    }
    if (WANTU2 && M - P > 0) {
      for (J = 2; J <= M - P; J++) {
        U2[1][J] = Complex.zero;
      }
      zlacpy('L', M - P - 1, M - Q - 1, X21(2, 1), LDX21, U2(2, 2), LDU2);
      zungqr(M - P, M - P, M - Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR,
          CHILDINFO);
    }
    if (WANTV1T && Q > 0) {
      zlacpy('U', M - Q, Q, X21, LDX21, V1T, LDV1T);
      zlacpy('U', P - (M - Q), Q - (M - Q), X11(M - Q + 1, M - Q + 1), LDX11,
          V1T(M - Q + 1, M - Q + 1), LDV1T);
      zlacpy('U', -P + Q, Q - P, X21(M - Q + 1, P + 1), LDX21,
          V1T(P + 1, P + 1), LDV1T);
      zunglq(
          Q, Q, Q, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO);
    }

    // Simultaneously diagonalize X11 and X21.

    zbbcsd(
        JOBU2,
        JOBU1,
        'N',
        JOBV1T,
        'N',
        M,
        M - P,
        M - Q,
        THETA,
        RWORK(IPHI),
        U2,
        LDU2,
        U1,
        LDU1,
        CDUM,
        1,
        V1T,
        LDV1T,
        RWORK(IB11D),
        RWORK(IB11E),
        RWORK(IB12D),
        RWORK(IB12E),
        RWORK(IB21D),
        RWORK(IB21E),
        RWORK(IB22D),
        RWORK(IB22E),
        RWORK(IBBCSD),
        LBBCSD,
        CHILDINFO);

    // Permute rows and columns to place identity submatrices in
    // preferred positions

    if (P > R) {
      for (I = 1; I <= R; I++) {
        IWORK[I] = P - R + I;
      }
      for (I = R + 1; I <= P; I++) {
        IWORK[I] = I - R;
      }
      if (WANTU1) {
        zlapmt(false, P, P, U1, LDU1, IWORK);
      }
      if (WANTV1T) {
        zlapmr(false, P, Q, V1T, LDV1T, IWORK);
      }
    }
  }
}
