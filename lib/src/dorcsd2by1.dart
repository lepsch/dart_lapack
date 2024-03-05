import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dbbcsd.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlapmr.dart';
import 'package:lapack/src/dlapmt.dart';
import 'package:lapack/src/dorbdb1.dart';
import 'package:lapack/src/dorbdb2.dart';
import 'package:lapack/src/dorbdb3.dart';
import 'package:lapack/src/dorbdb4.dart';
import 'package:lapack/src/dorglq.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dorcsd2by1(
  final String JOBU1,
  final String JOBU2,
  final String JOBV1T,
  final int M,
  final int P,
  final int Q,
  final Matrix<double> X11_,
  final int LDX11,
  final Matrix<double> X21_,
  final int LDX21,
  final Array<double> THETA_,
  final Matrix<double> U1_,
  final int LDU1,
  final Matrix<double> U2_,
  final int LDU2,
  final Matrix<double> V1T_,
  final int LDV1T,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine (3.5.0) --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X11 = X11_.having(ld: LDX11);
  final X21 = X21_.having(ld: LDX21);
  final THETA = THETA_.having();
  final U1 = U1_.having(ld: LDU1);
  final U2 = U2_.having(ld: LDU2);
  final V1T = V1T_.having(ld: LDV1T);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
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
      LORGLQ,
      LORGLQMIN,
      LORGLQOPT,
      LORGQR,
      LORGQRMIN,
      LORGQROPT,
      LWORKMIN,
      LWORKOPT,
      R;
  bool LQUERY, WANTU1, WANTU2, WANTV1T;
  final DUM1 = Array<double>(1);
  final DUM2 = Matrix<double>(1, 1);
  final CHILDINFO = Box(0);

  // Test input arguments

  INFO.value = 0;
  WANTU1 = lsame(JOBU1, 'Y');
  WANTU2 = lsame(JOBU2, 'Y');
  WANTV1T = lsame(JOBV1T, 'Y');
  LQUERY = LWORK == -1;

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

  // WORK layout:
  // |-------------------------------------------------------|
  // | LWORKOPT (1)                                          |
  // |-------------------------------------------------------|
  // | PHI (max(1,R-1))                                      |
  // |-------------------------------------------------------|
  // | TAUP1 (max(1,P))                        | B11D (R)    |
  // | TAUP2 (max(1,M-P))                      | B11E (R-1)  |
  // | TAUQ1 (max(1,Q))                        | B12D (R)    |
  // |-----------------------------------------| B12E (R-1)  |
  // | DORBDB WORK | DORGQR WORK | DORGLQ WORK | B21D (R)    |
  // |             |             |             | B21E (R-1)  |
  // |             |             |             | B22D (R)    |
  // |             |             |             | B22E (R-1)  |
  // |             |             |             | DBBCSD WORK |
  // |-------------------------------------------------------|

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
    ITAUP1 = IPHI + max(1, R - 1);
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
      dorbdb1(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM1, DUM1, DUM1, DUM1,
          WORK, -1, CHILDINFO);
      LORBDB = WORK[1].toInt();
      if (WANTU1 && P > 0) {
        dorgqr(P, P, Q, U1, LDU1, DUM1, WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTU2 && M - P > 0) {
        dorgqr(M - P, M - P, Q, U2, LDU2, DUM1, WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, M - P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTV1T && Q > 0) {
        dorglq(Q - 1, Q - 1, Q - 1, V1T, LDV1T, DUM1, WORK(1), -1, CHILDINFO);
        LORGLQMIN = max(LORGLQMIN, Q - 1);
        LORGLQOPT = max(LORGLQOPT, WORK[1].toInt());
      }
      dbbcsd(
          JOBU1,
          JOBU2,
          JOBV1T,
          'N',
          'N',
          M,
          P,
          Q,
          THETA,
          DUM1,
          U1,
          LDU1,
          U2,
          LDU2,
          V1T,
          LDV1T,
          DUM2,
          1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          WORK(1),
          -1,
          CHILDINFO);
      LBBCSD = WORK[1].toInt();
    } else if (R == P) {
      dorbdb2(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM1, DUM1, DUM1, DUM1,
          WORK(1), -1, CHILDINFO);
      LORBDB = WORK[1].toInt();
      if (WANTU1 && P > 0) {
        dorgqr(
            P - 1, P - 1, P - 1, U1(2, 2), LDU1, DUM1, WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, P - 1);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTU2 && M - P > 0) {
        dorgqr(M - P, M - P, Q, U2, LDU2, DUM1, WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, M - P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTV1T && Q > 0) {
        dorglq(Q, Q, R, V1T, LDV1T, DUM1, WORK(1), -1, CHILDINFO);
        LORGLQMIN = max(LORGLQMIN, Q);
        LORGLQOPT = max(LORGLQOPT, WORK[1].toInt());
      }
      dbbcsd(
          JOBV1T,
          'N',
          JOBU1,
          JOBU2,
          'T',
          M,
          Q,
          P,
          THETA,
          DUM1,
          V1T,
          LDV1T,
          DUM2,
          1,
          U1,
          LDU1,
          U2,
          LDU2,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          WORK(1),
          -1,
          CHILDINFO);
      LBBCSD = WORK[1].toInt();
    } else if (R == M - P) {
      dorbdb3(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM1, DUM1, DUM1, DUM1,
          WORK(1), -1, CHILDINFO);
      LORBDB = WORK[1].toInt();
      if (WANTU1 && P > 0) {
        dorgqr(P, P, Q, U1, LDU1, DUM1, WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTU2 && M - P > 0) {
        dorgqr(M - P - 1, M - P - 1, M - P - 1, U2(2, 2), LDU2, DUM1, WORK(1),
            -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, M - P - 1);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTV1T && Q > 0) {
        dorglq(Q, Q, R, V1T, LDV1T, DUM1, WORK(1), -1, CHILDINFO);
        LORGLQMIN = max(LORGLQMIN, Q);
        LORGLQOPT = max(LORGLQOPT, WORK[1].toInt());
      }
      dbbcsd(
          'N',
          JOBV1T,
          JOBU2,
          JOBU1,
          'T',
          M,
          M - Q,
          M - P,
          THETA,
          DUM1,
          DUM2,
          1,
          V1T,
          LDV1T,
          U2,
          LDU2,
          U1,
          LDU1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          WORK(1),
          -1,
          CHILDINFO);
      LBBCSD = WORK[1].toInt();
    } else {
      dorbdb4(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM1, DUM1, DUM1, DUM1,
          DUM1, WORK(1), -1, CHILDINFO);
      LORBDB = M + WORK[1].toInt();
      if (WANTU1 && P > 0) {
        dorgqr(P, P, M - Q, U1, LDU1, DUM1, WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTU2 && M - P > 0) {
        dorgqr(M - P, M - P, M - Q, U2, LDU2, DUM1, WORK(1), -1, CHILDINFO);
        LORGQRMIN = max(LORGQRMIN, M - P);
        LORGQROPT = max(LORGQROPT, WORK[1].toInt());
      }
      if (WANTV1T && Q > 0) {
        dorglq(Q, Q, Q, V1T, LDV1T, DUM1, WORK(1), -1, CHILDINFO);
        LORGLQMIN = max(LORGLQMIN, Q);
        LORGLQOPT = max(LORGLQOPT, WORK[1].toInt());
      }
      dbbcsd(
          JOBU2,
          JOBU1,
          'N',
          JOBV1T,
          'N',
          M,
          M - P,
          M - Q,
          THETA,
          DUM1,
          U2,
          LDU2,
          U1,
          LDU1,
          DUM2,
          1,
          V1T,
          LDV1T,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          DUM1,
          WORK(1),
          -1,
          CHILDINFO);
      LBBCSD = WORK[1].toInt();
    }
    LWORKMIN = max(max(IORBDB + LORBDB - 1, IORGQR + LORGQRMIN - 1),
        max(IORGLQ + LORGLQMIN - 1, IBBCSD + LBBCSD - 1));
    LWORKOPT = max(
      max(IORBDB + LORBDB - 1, IORGQR + LORGQROPT - 1),
      max(IORGLQ + LORGLQOPT - 1, IBBCSD + LBBCSD - 1),
    );
    WORK[1] = LWORKOPT.toDouble();
    if (LWORK < LWORKMIN && !LQUERY) {
      INFO.value = -19;
    }
  }
  if (INFO.value != 0) {
    xerbla('DORCSD2BY1', -INFO.value);
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

    dorbdb1(M, P, Q, X11, LDX11, X21, LDX21, THETA, WORK(IPHI), WORK(ITAUP1),
        WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO);

    // Accumulate Householder reflectors

    if (WANTU1 && P > 0) {
      dlacpy('L', P, Q, X11, LDX11, U1, LDU1);
      dorgqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO);
    }
    if (WANTU2 && M - P > 0) {
      dlacpy('L', M - P, Q, X21, LDX21, U2, LDU2);
      dorgqr(M - P, M - P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR,
          CHILDINFO);
    }
    if (WANTV1T && Q > 0) {
      V1T[1][1] = ONE;
      for (J = 2; J <= Q; J++) {
        V1T[1][J] = ZERO;
        V1T[J][1] = ZERO;
      }
      dlacpy('U', Q - 1, Q - 1, X21(1, 2), LDX21, V1T(2, 2), LDV1T);
      dorglq(Q - 1, Q - 1, Q - 1, V1T(2, 2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ),
          LORGLQ, CHILDINFO);
    }

    // Simultaneously diagonalize X11 and X21.

    dbbcsd(
        JOBU1,
        JOBU2,
        JOBV1T,
        'N',
        'N',
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
        DUM2,
        1,
        WORK(IB11D),
        WORK(IB11E),
        WORK(IB12D),
        WORK(IB12E),
        WORK(IB21D),
        WORK(IB21E),
        WORK(IB22D),
        WORK(IB22E),
        WORK(IBBCSD),
        LBBCSD,
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
      dlapmt(false, M - P, M - P, U2, LDU2, IWORK);
    }
  } else if (R == P) {
    // Case 2: R = P

    // Simultaneously bidiagonalize X11 and X21

    dorbdb2(M, P, Q, X11, LDX11, X21, LDX21, THETA, WORK(IPHI), WORK(ITAUP1),
        WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO);

    // Accumulate Householder reflectors

    if (WANTU1 && P > 0) {
      U1[1][1] = ONE;
      for (J = 2; J <= P; J++) {
        U1[1][J] = ZERO;
        U1[J][1] = ZERO;
      }
      dlacpy('L', P - 1, P - 1, X11(2, 1), LDX11, U1(2, 2), LDU1);
      dorgqr(P - 1, P - 1, P - 1, U1(2, 2), LDU1, WORK(ITAUP1), WORK(IORGQR),
          LORGQR, CHILDINFO);
    }
    if (WANTU2 && M - P > 0) {
      dlacpy('L', M - P, Q, X21, LDX21, U2, LDU2);
      dorgqr(M - P, M - P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR,
          CHILDINFO);
    }
    if (WANTV1T && Q > 0) {
      dlacpy('U', P, Q, X11, LDX11, V1T, LDV1T);
      dorglq(
          Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO);
    }

    // Simultaneously diagonalize X11 and X21.

    dbbcsd(
        JOBV1T,
        'N',
        JOBU1,
        JOBU2,
        'T',
        M,
        Q,
        P,
        THETA,
        WORK(IPHI),
        V1T,
        LDV1T,
        DUM1.asMatrix(1),
        1,
        U1,
        LDU1,
        U2,
        LDU2,
        WORK(IB11D),
        WORK(IB11E),
        WORK(IB12D),
        WORK(IB12E),
        WORK(IB21D),
        WORK(IB21E),
        WORK(IB22D),
        WORK(IB22E),
        WORK(IBBCSD),
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
      dlapmt(false, M - P, M - P, U2, LDU2, IWORK);
    }
  } else if (R == M - P) {
    // Case 3: R = M-P

    // Simultaneously bidiagonalize X11 and X21

    dorbdb3(M, P, Q, X11, LDX11, X21, LDX21, THETA, WORK(IPHI), WORK(ITAUP1),
        WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO);

    // Accumulate Householder reflectors

    if (WANTU1 && P > 0) {
      dlacpy('L', P, Q, X11, LDX11, U1, LDU1);
      dorgqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO);
    }
    if (WANTU2 && M - P > 0) {
      U2[1][1] = ONE;
      for (J = 2; J <= M - P; J++) {
        U2[1][J] = ZERO;
        U2[J][1] = ZERO;
      }
      dlacpy('L', M - P - 1, M - P - 1, X21(2, 1), LDX21, U2(2, 2), LDU2);
      dorgqr(M - P - 1, M - P - 1, M - P - 1, U2(2, 2), LDU2, WORK(ITAUP2),
          WORK(IORGQR), LORGQR, CHILDINFO);
    }
    if (WANTV1T && Q > 0) {
      dlacpy('U', M - P, Q, X21, LDX21, V1T, LDV1T);
      dorglq(
          Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO);
    }

    // Simultaneously diagonalize X11 and X21.

    dbbcsd(
        'N',
        JOBV1T,
        JOBU2,
        JOBU1,
        'T',
        M,
        M - Q,
        M - P,
        THETA,
        WORK(IPHI),
        DUM1.asMatrix(1),
        1,
        V1T,
        LDV1T,
        U2,
        LDU2,
        U1,
        LDU1,
        WORK(IB11D),
        WORK(IB11E),
        WORK(IB12D),
        WORK(IB12E),
        WORK(IB21D),
        WORK(IB21E),
        WORK(IB22D),
        WORK(IB22E),
        WORK(IBBCSD),
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
        dlapmt(false, P, Q, U1, LDU1, IWORK);
      }
      if (WANTV1T) {
        dlapmr(false, Q, Q, V1T, LDV1T, IWORK);
      }
    }
  } else {
    // Case 4: R = M-Q

    // Simultaneously bidiagonalize X11 and X21

    dorbdb4(
        M,
        P,
        Q,
        X11,
        LDX11,
        X21,
        LDX21,
        THETA,
        WORK(IPHI),
        WORK(ITAUP1),
        WORK(ITAUP2),
        WORK(ITAUQ1),
        WORK(IORBDB),
        WORK(IORBDB + M),
        LORBDB - M,
        CHILDINFO);

    // Accumulate Householder reflectors

    if (WANTU2 && M - P > 0) {
      dcopy(M - P, WORK(IORBDB + P), 1, U2.asArray(), 1);
    }
    if (WANTU1 && P > 0) {
      dcopy(P, WORK(IORBDB), 1, U1.asArray(), 1);
      for (J = 2; J <= P; J++) {
        U1[1][J] = ZERO;
      }
      dlacpy('L', P - 1, M - Q - 1, X11(2, 1), LDX11, U1(2, 2), LDU1);
      dorgqr(
          P, P, M - Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO);
    }
    if (WANTU2 && M - P > 0) {
      for (J = 2; J <= M - P; J++) {
        U2[1][J] = ZERO;
      }
      dlacpy('L', M - P - 1, M - Q - 1, X21(2, 1), LDX21, U2(2, 2), LDU2);
      dorgqr(M - P, M - P, M - Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR,
          CHILDINFO);
    }
    if (WANTV1T && Q > 0) {
      dlacpy('U', M - Q, Q, X21, LDX21, V1T, LDV1T);
      dlacpy('U', P - (M - Q), Q - (M - Q), X11(M - Q + 1, M - Q + 1), LDX11,
          V1T(M - Q + 1, M - Q + 1), LDV1T);
      dlacpy('U', -P + Q, Q - P, X21(M - Q + 1, P + 1), LDX21,
          V1T(P + 1, P + 1), LDV1T);
      dorglq(
          Q, Q, Q, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO);
    }

    // Simultaneously diagonalize X11 and X21.

    dbbcsd(
        JOBU2,
        JOBU1,
        'N',
        JOBV1T,
        'N',
        M,
        M - P,
        M - Q,
        THETA,
        WORK(IPHI),
        U2,
        LDU2,
        U1,
        LDU1,
        DUM1.asMatrix(1),
        1,
        V1T,
        LDV1T,
        WORK(IB11D),
        WORK(IB11E),
        WORK(IB12D),
        WORK(IB12E),
        WORK(IB21D),
        WORK(IB21E),
        WORK(IB22D),
        WORK(IB22E),
        WORK(IBBCSD),
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
        dlapmt(false, P, P, U1, LDU1, IWORK);
      }
      if (WANTV1T) {
        dlapmr(false, P, Q, V1T, LDV1T, IWORK);
      }
    }
  }
}
