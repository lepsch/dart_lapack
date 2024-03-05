import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlag2.dart';
import 'package:lapack/src/dlaqz0.dart';
import 'package:lapack/src/dlaqz2.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dtgexc.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlaqz3(
  final bool ILSCHUR,
  final bool ILQ,
  final bool ILZ,
  final int N,
  final int ILO,
  final int IHI,
  final int NW,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> Z_,
  final int LDZ,
  final Box<int> NS,
  final Box<int> ND,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Matrix<double> QC_,
  final int LDQC,
  final Matrix<double> ZC_,
  final int LDZC,
  final Array<double> WORK_,
  final int LWORK,
  final int REC,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final QC = QC_.having(ld: LDQC);
  final ZC = ZC_.having(ld: LDZC);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool BULGE;
  int KWBOT, ISTOPM, ISTARTM, K, K2, LWORKREQ = 0;
  double S, SMLNUM, ULP, SAFMIN
      // SAFMAX
      ;
  final DTGEXC_INFO = Box(0),
      QZ_SMALL_INFO = Box(0),
      IFST = Box(0),
      ILST = Box(0);
  final C1 = Box(0.0), S1 = Box(0.0), TEMP = Box(0.0);

  INFO.value = 0;

  // Set up deflation window
  final JW = min(NW, IHI - ILO + 1);
  final KWTOP = IHI - JW + 1;
  if (KWTOP == ILO) {
    S = ZERO;
  } else {
    S = A[KWTOP][KWTOP - 1];
  }

  // Determine required workspace
  IFST.value = 1;
  ILST.value = JW;
  dtgexc(true, true, JW, A, LDA, B, LDB, QC, LDQC, ZC, LDZC, IFST, ILST, WORK,
      -1, DTGEXC_INFO);
  LWORKREQ = WORK[1].toInt();
  dlaqz0(
      'S',
      'V',
      'V',
      JW,
      1,
      JW,
      A(KWTOP, KWTOP),
      LDA,
      B(KWTOP, KWTOP),
      LDB,
      ALPHAR,
      ALPHAI,
      BETA,
      QC,
      LDQC,
      ZC,
      LDZC,
      WORK,
      -1,
      REC + 1,
      QZ_SMALL_INFO);
  LWORKREQ = max(LWORKREQ, WORK[1].toInt() + 2 * pow(JW, 2).toInt());
  LWORKREQ = max(LWORKREQ, max(N * NW, 2 * pow(NW, 2).toInt() + N));
  if (LWORK == -1) {
    // workspace query, quick return;
    WORK[1] = LWORKREQ.toDouble();
    return;
  } else if (LWORK < LWORKREQ) {
    INFO.value = -26;
  }

  if (INFO.value != 0) {
    xerbla('DLAQZ3', -INFO.value);
    return;
  }

  // Get machine constants
  SAFMIN = dlamch('SAFE MINIMUM');
  // SAFMAX = ONE / SAFMIN;
  ULP = dlamch('PRECISION');
  SMLNUM = SAFMIN * (N.toDouble() / ULP);

  if (IHI == KWTOP) {
    // 1 by 1 deflation window, just try a regular deflation
    ALPHAR[KWTOP] = A[KWTOP][KWTOP];
    ALPHAI[KWTOP] = ZERO;
    BETA[KWTOP] = B[KWTOP][KWTOP];
    NS.value = 1;
    ND.value = 0;
    if ((S).abs() <= max(SMLNUM, ULP * (A[KWTOP][KWTOP]).abs())) {
      NS.value = 0;
      ND.value = 1;
      if (KWTOP > ILO) {
        A[KWTOP][KWTOP - 1] = ZERO;
      }
    }
  }

  // Store window in case of convergence failure
  dlacpy('ALL', JW, JW, A(KWTOP, KWTOP), LDA, WORK.asMatrix(JW), JW);
  dlacpy('ALL', JW, JW, B(KWTOP, KWTOP), LDB,
      WORK(pow(JW, 2).toInt() + 1).asMatrix(JW), JW);

  // Transform window to real schur form
  dlaset('FULL', JW, JW, ZERO, ONE, QC, LDQC);
  dlaset('FULL', JW, JW, ZERO, ONE, ZC, LDZC);
  dlaqz0(
      'S',
      'V',
      'V',
      JW,
      1,
      JW,
      A(KWTOP, KWTOP),
      LDA,
      B(KWTOP, KWTOP),
      LDB,
      ALPHAR,
      ALPHAI,
      BETA,
      QC,
      LDQC,
      ZC,
      LDZC,
      WORK(2 * pow(JW, 2).toInt() + 1),
      LWORK - 2 * pow(JW, 2).toInt(),
      REC + 1,
      QZ_SMALL_INFO);

  if (QZ_SMALL_INFO.value != 0) {
    // Convergence failure, restore the window and exit
    ND.value = 0;
    NS.value = JW - QZ_SMALL_INFO.value;
    dlacpy('ALL', JW, JW, WORK.asMatrix(JW), JW, A(KWTOP, KWTOP), LDA);
    dlacpy('ALL', JW, JW, WORK(pow(JW, 2).toInt() + 1).asMatrix(JW), JW,
        B(KWTOP, KWTOP), LDB);
    return;
  }

  // Deflation detection loop
  if (KWTOP == ILO || S == ZERO) {
    KWBOT = KWTOP - 1;
  } else {
    KWBOT = IHI;
    K = 1;
    K2 = 1;
    while (K <= JW) {
      BULGE = false;
      if (KWBOT - KWTOP + 1 >= 2) {
        BULGE = A[KWBOT][KWBOT - 1] != ZERO;
      }
      if (BULGE) {
        // Try to deflate complex conjugate eigenvalue pair
        TEMP.value = A[KWBOT][KWBOT].abs() +
            sqrt(A[KWBOT][KWBOT - 1].abs()) * sqrt((A[KWBOT - 1][KWBOT]).abs());
        if (TEMP.value == ZERO) {
          TEMP.value = (S).abs();
        }
        if (max((S * QC[1][KWBOT - KWTOP]).abs(),
                (S * QC[1][KWBOT - KWTOP + 1]).abs()) <=
            max(SMLNUM, ULP * TEMP.value)) {
          // Deflatable
          KWBOT = KWBOT - 2;
        } else {
          // Not deflatable, move out of the way
          IFST.value = KWBOT - KWTOP + 1;
          ILST.value = K2;
          dtgexc(true, true, JW, A(KWTOP, KWTOP), LDA, B(KWTOP, KWTOP), LDB, QC,
              LDQC, ZC, LDZC, IFST, ILST, WORK, LWORK, DTGEXC_INFO);
          K2 = K2 + 2;
        }
        K = K + 2;
      } else {
        // Try to deflate real eigenvalue
        TEMP.value = A[KWBOT][KWBOT].abs();
        if (TEMP.value == ZERO) {
          TEMP.value = (S).abs();
        }
        if (((S * QC[1][KWBOT - KWTOP + 1]).abs()) <=
            max(ULP * TEMP.value, SMLNUM)) {
          // Deflatable
          KWBOT = KWBOT - 1;
        } else {
          // Not deflatable, move out of the way
          IFST.value = KWBOT - KWTOP + 1;
          ILST.value = K2;
          dtgexc(true, true, JW, A(KWTOP, KWTOP), LDA, B(KWTOP, KWTOP), LDB, QC,
              LDQC, ZC, LDZC, IFST, ILST, WORK, LWORK, DTGEXC_INFO);
          K2 = K2 + 1;
        }

        K = K + 1;
      }
    }
  }

  // Store eigenvalues
  ND.value = IHI - KWBOT;
  NS.value = JW - ND.value;
  K = KWTOP;
  while (K <= IHI) {
    BULGE = false;
    if (K < IHI) {
      if (A[K + 1][K] != ZERO) {
        BULGE = true;
      }
    }
    if (BULGE) {
      // 2x2 eigenvalue block
      dlag2(A(K, K), LDA, B(K, K), LDB, SAFMIN, BETA.box(K), BETA.box(K + 1),
          ALPHAR.box(K), ALPHAR.box(K + 1), ALPHAI.box(K));
      ALPHAI[K + 1] = -ALPHAI[K];
      K = K + 2;
    } else {
      // 1x1 eigenvalue block
      ALPHAR[K] = A[K][K];
      ALPHAI[K] = ZERO;
      BETA[K] = B[K][K];
      K = K + 1;
    }
  }

  if (KWTOP != ILO && S != ZERO) {
    // Reflect spike back, this will create optimally packed bulges
    final v = A[KWTOP][KWTOP - 1];
    for (var i = KWTOP, j = 1; i <= KWBOT && j <= JW - ND.value; i++, j++) {
      A[i][KWTOP - 1] = v * QC[1][j];
    }
    for (K = KWBOT - 1; K >= KWTOP; K--) {
      dlartg(A[K][KWTOP - 1], A[K + 1][KWTOP - 1], C1, S1, TEMP);
      A[K][KWTOP - 1] = TEMP.value;
      A[K + 1][KWTOP - 1] = ZERO;
      K2 = max(KWTOP, K - 1);
      drot(IHI - K2 + 1, A(K, K2).asArray(), LDA, A(K + 1, K2).asArray(), LDA,
          C1.value, S1.value);
      drot(IHI - (K - 1) + 1, B(K, K - 1).asArray(), LDB,
          B(K + 1, K - 1).asArray(), LDB, C1.value, S1.value);
      drot(JW, QC(1, K - KWTOP + 1).asArray(), 1,
          QC(1, K + 1 - KWTOP + 1).asArray(), 1, C1.value, S1.value);
    }

    // Chase bulges down
    ISTARTM = KWTOP;
    ISTOPM = IHI;
    K = KWBOT - 1;
    while (K >= KWTOP) {
      if ((K >= KWTOP + 1) && A[K + 1][K - 1] != ZERO) {
        // Move double pole block down and remove it
        for (K2 = K - 1; K2 <= KWBOT - 2; K2++) {
          dlaqz2(true, true, K2, KWTOP, KWTOP + JW - 1, KWBOT, A, LDA, B, LDB,
              JW, KWTOP, QC, LDQC, JW, KWTOP, ZC, LDZC);
        }

        K = K - 2;
      } else {
        // k points to single shift
        for (K2 = K; K2 <= KWBOT - 2; K2++) {
          // Move shift down
          dlartg(B[K2 + 1][K2 + 1], B[K2 + 1][K2], C1, S1, TEMP);
          B[K2 + 1][K2 + 1] = TEMP.value;
          B[K2 + 1][K2] = ZERO;
          drot(K2 + 2 - ISTARTM + 1, A(ISTARTM, K2 + 1).asArray(), 1,
              A(ISTARTM, K2).asArray(), 1, C1.value, S1.value);
          drot(K2 - ISTARTM + 1, B(ISTARTM, K2 + 1).asArray(), 1,
              B(ISTARTM, K2).asArray(), 1, C1.value, S1.value);
          drot(JW, ZC(1, K2 + 1 - KWTOP + 1).asArray(), 1,
              ZC(1, K2 - KWTOP + 1).asArray(), 1, C1.value, S1.value);

          dlartg(A[K2 + 1][K2], A[K2 + 2][K2], C1, S1, TEMP);
          A[K2 + 1][K2] = TEMP.value;
          A[K2 + 2][K2] = ZERO;
          drot(ISTOPM - K2, A(K2 + 1, K2 + 1).asArray(), LDA,
              A(K2 + 2, K2 + 1).asArray(), LDA, C1.value, S1.value);
          drot(ISTOPM - K2, B(K2 + 1, K2 + 1).asArray(), LDB,
              B(K2 + 2, K2 + 1).asArray(), LDB, C1.value, S1.value);
          drot(JW, QC(1, K2 + 1 - KWTOP + 1).asArray(), 1,
              QC(1, K2 + 2 - KWTOP + 1).asArray(), 1, C1.value, S1.value);
        }

        // Remove the shift
        dlartg(B[KWBOT][KWBOT], B[KWBOT][KWBOT - 1], C1, S1, TEMP);
        B[KWBOT][KWBOT] = TEMP.value;
        B[KWBOT][KWBOT - 1] = ZERO;
        drot(KWBOT - ISTARTM, B(ISTARTM, KWBOT).asArray(), 1,
            B(ISTARTM, KWBOT - 1).asArray(), 1, C1.value, S1.value);
        drot(KWBOT - ISTARTM + 1, A(ISTARTM, KWBOT).asArray(), 1,
            A(ISTARTM, KWBOT - 1).asArray(), 1, C1.value, S1.value);
        drot(JW, ZC(1, KWBOT - KWTOP + 1).asArray(), 1,
            ZC(1, KWBOT - 1 - KWTOP + 1).asArray(), 1, C1.value, S1.value);

        K = K - 1;
      }
    }
  }

  // Apply Qc and Zc to rest of the matrix
  if (ILSCHUR) {
    ISTARTM = 1;
    ISTOPM = N;
  } else {
    ISTARTM = ILO;
    ISTOPM = IHI;
  }

  if (ISTOPM - IHI > 0) {
    dgemm('T', 'N', JW, ISTOPM - IHI, JW, ONE, QC, LDQC, A(KWTOP, IHI + 1), LDA,
        ZERO, WORK.asMatrix(JW), JW);
    dlacpy(
        'ALL', JW, ISTOPM - IHI, WORK.asMatrix(JW), JW, A(KWTOP, IHI + 1), LDA);
    dgemm('T', 'N', JW, ISTOPM - IHI, JW, ONE, QC, LDQC, B(KWTOP, IHI + 1), LDB,
        ZERO, WORK.asMatrix(JW), JW);
    dlacpy(
        'ALL', JW, ISTOPM - IHI, WORK.asMatrix(JW), JW, B(KWTOP, IHI + 1), LDB);
  }
  if (ILQ) {
    dgemm('N', 'N', N, JW, JW, ONE, Q(1, KWTOP), LDQ, QC, LDQC, ZERO,
        WORK.asMatrix(N), N);
    dlacpy('ALL', N, JW, WORK.asMatrix(N), N, Q(1, KWTOP), LDQ);
  }

  if (KWTOP - 1 - ISTARTM + 1 > 0) {
    dgemm('N', 'N', KWTOP - ISTARTM, JW, JW, ONE, A(ISTARTM, KWTOP), LDA, ZC,
        LDZC, ZERO, WORK.asMatrix(KWTOP), KWTOP - ISTARTM);
    dlacpy('ALL', KWTOP - ISTARTM, JW, WORK.asMatrix(KWTOP), KWTOP - ISTARTM,
        A(ISTARTM, KWTOP), LDA);
    dgemm('N', 'N', KWTOP - ISTARTM, JW, JW, ONE, B(ISTARTM, KWTOP), LDB, ZC,
        LDZC, ZERO, WORK.asMatrix(KWTOP), KWTOP - ISTARTM);
    dlacpy('ALL', KWTOP - ISTARTM, JW, WORK.asMatrix(KWTOP), KWTOP - ISTARTM,
        B(ISTARTM, KWTOP), LDB);
  }
  if (ILZ) {
    dgemm('N', 'N', N, JW, JW, ONE, Z(1, KWTOP), LDZ, ZC, LDZC, ZERO,
        WORK.asMatrix(N), N);
    dlacpy('ALL', N, JW, WORK.asMatrix(N), N, Z(1, KWTOP), LDZ);
  }
}
