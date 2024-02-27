import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaqz0.dart';
import 'package:lapack/src/zlaqz1.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zrot.dart';
import 'package:lapack/src/ztgexc.dart';

void zlaqz2(
  final bool ILSCHUR,
  final bool ILQ,
  final bool ILZ,
  final int N,
  final int ILO,
  final int IHI,
  final int NW,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Box<int> NS,
  final Box<int> ND,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Matrix<Complex> QC_,
  final int LDQC,
  final Matrix<Complex> ZC_,
  final int LDZC,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int REC,
  final Box<int> INFO,
) {
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final Q = Q_.dim(LDQ);
  final Z = Z_.dim(LDZ);
  final ALPHA = ALPHA_.dim();
  final BETA = BETA_.dim();
  final QC = QC_.dim(LDQC);
  final ZC = ZC_.dim(LDZC);
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  const ZERO = 0.0;
  int JW, KWTOP, KWBOT, ISTOPM, ISTARTM, K, K2, IFST, LWORKREQ;
  double SMLNUM, ULP, SAFMIN, TEMPR;
  Complex S;
  final QZ_SMALL_INFO = Box(0), ZTGEXC_INFO = Box(0), ILST = Box(0);
  final C1 = Box(0.0);
  final S1 = Box(Complex.zero), TEMP = Box(Complex.zero);

  INFO.value = 0;

  // Set up deflation window
  JW = min(NW, IHI - ILO + 1);
  KWTOP = IHI - JW + 1;
  if (KWTOP == ILO) {
    S = Complex.zero;
  } else {
    S = A[KWTOP][KWTOP - 1];
  }

  // Determine required workspace
  IFST = 1;
  ILST.value = JW;
  zlaqz0('S', 'V', 'V', JW, 1, JW, A(KWTOP, KWTOP), LDA, B(KWTOP, KWTOP), LDB,
      ALPHA, BETA, QC, LDQC, ZC, LDZC, WORK, -1, RWORK, REC + 1, QZ_SMALL_INFO);
  LWORKREQ = WORK[1].toInt() + 2 * pow(JW, 2).toInt();
  LWORKREQ = max(LWORKREQ, max(N * NW, 2 * pow(NW, 2).toInt() + N));
  if (LWORK == -1) {
    // workspace query, quick return;
    WORK[1] = LWORKREQ.toComplex();
    return;
  } else if (LWORK < LWORKREQ) {
    INFO.value = -26;
  }

  if (INFO.value != 0) {
    xerbla('ZLAQZ2', -INFO.value);
    return;
  }

  // Get machine constants
  SAFMIN = dlamch('SAFE MINIMUM');
  ULP = dlamch('PRECISION');
  SMLNUM = SAFMIN * (N.toDouble() / ULP);

  if (IHI == KWTOP) {
    // 1 by 1 deflation window, just try a regular deflation
    ALPHA[KWTOP] = A[KWTOP][KWTOP];
    BETA[KWTOP] = B[KWTOP][KWTOP];
    NS.value = 1;
    ND.value = 0;
    if (S.abs() <= max(SMLNUM, ULP * (A[KWTOP][KWTOP]).abs())) {
      NS.value = 0;
      ND.value = 1;
      if (KWTOP > ILO) {
        A[KWTOP][KWTOP - 1] = Complex.zero;
      }
    }
  }

  // Store window in case of convergence failure
  zlacpy('ALL', JW, JW, A(KWTOP, KWTOP), LDA, WORK.asMatrix(JW), JW);
  zlacpy('ALL', JW, JW, B(KWTOP, KWTOP), LDB,
      WORK(pow(JW, 2).toInt() + 1).asMatrix(JW), JW);

  // Transform window to real schur form
  zlaset('FULL', JW, JW, Complex.zero, Complex.one, QC, LDQC);
  zlaset('FULL', JW, JW, Complex.zero, Complex.one, ZC, LDZC);
  zlaqz0(
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
      ALPHA,
      BETA,
      QC,
      LDQC,
      ZC,
      LDZC,
      WORK(2 * pow(JW, 2).toInt() + 1),
      LWORK - 2 * pow(JW, 2).toInt(),
      RWORK,
      REC + 1,
      QZ_SMALL_INFO);

  if (QZ_SMALL_INFO.value != 0) {
    // Convergence failure, restore the window and exit
    ND.value = 0;
    NS.value = JW - QZ_SMALL_INFO.value;
    zlacpy('ALL', JW, JW, WORK.asMatrix(JW), JW, A(KWTOP, KWTOP), LDA);
    zlacpy('ALL', JW, JW, WORK(pow(JW, 2).toInt() + 1).asMatrix(JW), JW,
        B(KWTOP, KWTOP), LDB);
    return;
  }

  // Deflation detection loop
  if (KWTOP == ILO || S == Complex.zero) {
    KWBOT = KWTOP - 1;
  } else {
    KWBOT = IHI;
    K = 1;
    K2 = 1;
    while (K <= JW) {
      // Try to deflate eigenvalue
      TEMPR = (A[KWBOT][KWBOT]).abs();
      if (TEMPR == ZERO) {
        TEMPR = (S).abs();
      }
      if (((S * QC[1][KWBOT - KWTOP + 1]).abs()) <= max(ULP * TEMPR, SMLNUM)) {
        // Deflatable
        KWBOT = KWBOT - 1;
      } else {
        // Not deflatable, move out of the way
        IFST = KWBOT - KWTOP + 1;
        ILST.value = K2;
        ztgexc(true, true, JW, A(KWTOP, KWTOP), LDA, B(KWTOP, KWTOP), LDB, QC,
            LDQC, ZC, LDZC, IFST, ILST, ZTGEXC_INFO);
        K2 = K2 + 1;
      }

      K = K + 1;
    }
  }

  // Store eigenvalues
  ND.value = IHI - KWBOT;
  NS.value = JW - ND.value;
  K = KWTOP;
  while (K <= IHI) {
    ALPHA[K] = A[K][K];
    BETA[K] = B[K][K];
    K = K + 1;
  }

  if (KWTOP != ILO && S != Complex.zero) {
    // Reflect spike back, this will create optimally packed bulges
    for (var i = KWTOP; i <= KWBOT; i++) {
      A[i][KWTOP - 1] = A[KWTOP][KWTOP - 1] * QC[1][i - KWTOP + 1].conjugate();
    }
    for (K = KWBOT - 1; K >= KWTOP; K--) {
      zlartg(A[K][KWTOP - 1], A[K + 1][KWTOP - 1], C1, S1, TEMP);
      A[K][KWTOP - 1] = TEMP.value;
      A[K + 1][KWTOP - 1] = Complex.zero;
      K2 = max(KWTOP, K - 1);
      zrot(IHI - K2 + 1, A(K, K2).asArray(), LDA, A(K + 1, K2).asArray(), LDA,
          C1.value, S1.value);
      zrot(IHI - (K - 1) + 1, B(K, K - 1).asArray(), LDB,
          B(K + 1, K - 1).asArray(), LDB, C1.value, S1.value);
      zrot(
          JW,
          QC(1, K - KWTOP + 1).asArray(),
          1,
          QC(1, K + 1 - KWTOP + 1).asArray(),
          1,
          C1.value,
          S1.value.conjugate());
    }

    // Chase bulges down
    ISTARTM = KWTOP;
    ISTOPM = IHI;
    K = KWBOT - 1;
    while (K >= KWTOP) {
      // Move bulge down and remove it
      for (K2 = K; K2 <= KWBOT - 1; K2++) {
        zlaqz1(true, true, K2, KWTOP, KWTOP + JW - 1, KWBOT, A, LDA, B, LDB, JW,
            KWTOP, QC, LDQC, JW, KWTOP, ZC, LDZC);
      }

      K = K - 1;
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
    zgemm('C', 'N', JW, ISTOPM - IHI, JW, Complex.one, QC, LDQC,
        A(KWTOP, IHI + 1), LDA, Complex.zero, WORK.asMatrix(JW), JW);
    zlacpy(
        'ALL', JW, ISTOPM - IHI, WORK.asMatrix(JW), JW, A(KWTOP, IHI + 1), LDA);
    zgemm('C', 'N', JW, ISTOPM - IHI, JW, Complex.one, QC, LDQC,
        B(KWTOP, IHI + 1), LDB, Complex.zero, WORK.asMatrix(JW), JW);
    zlacpy(
        'ALL', JW, ISTOPM - IHI, WORK.asMatrix(JW), JW, B(KWTOP, IHI + 1), LDB);
  }
  if (ILQ) {
    zgemm('N', 'N', N, JW, JW, Complex.one, Q(1, KWTOP), LDQ, QC, LDQC,
        Complex.zero, WORK.asMatrix(N), N);
    zlacpy('ALL', N, JW, WORK.asMatrix(N), N, Q(1, KWTOP), LDQ);
  }

  if (KWTOP - 1 - ISTARTM + 1 > 0) {
    zgemm(
        'N',
        'N',
        KWTOP - ISTARTM,
        JW,
        JW,
        Complex.one,
        A(ISTARTM, KWTOP),
        LDA,
        ZC,
        LDZC,
        Complex.zero,
        WORK.asMatrix(KWTOP - ISTARTM),
        KWTOP - ISTARTM);
    zlacpy('ALL', KWTOP - ISTARTM, JW, WORK.asMatrix(KWTOP - ISTARTM),
        KWTOP - ISTARTM, A(ISTARTM, KWTOP), LDA);
    zgemm(
        'N',
        'N',
        KWTOP - ISTARTM,
        JW,
        JW,
        Complex.one,
        B(ISTARTM, KWTOP),
        LDB,
        ZC,
        LDZC,
        Complex.zero,
        WORK.asMatrix(KWTOP - ISTARTM),
        KWTOP - ISTARTM);
    zlacpy('ALL', KWTOP - ISTARTM, JW, WORK.asMatrix(KWTOP - ISTARTM),
        KWTOP - ISTARTM, B(ISTARTM, KWTOP), LDB);
  }
  if (ILZ) {
    zgemm('N', 'N', N, JW, JW, Complex.one, Z(1, KWTOP), LDZ, ZC, LDZC,
        Complex.zero, WORK.asMatrix(N), N);
    zlacpy('ALL', N, JW, WORK.asMatrix(N), N, Z(1, KWTOP), LDZ);
  }
}
