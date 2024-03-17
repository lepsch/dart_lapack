import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaqz1.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zrot.dart';

void zlaqz3(
  final bool ILSCHUR,
  final bool ILQ,
  final bool ILZ,
  final int N,
  final int ILO,
  final int IHI,
  final int NSHIFTS,
  final int NBLOCK_DESIRED,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Matrix<Complex> QC_,
  final int LDQC,
  final Matrix<Complex> ZC_,
  final int LDZC,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final QC = QC_.having(ld: LDQC);
  final ZC = ZC_.having(ld: LDZC);
  final WORK = WORK_.having();

  // Parameters
  const ONE = 1.0;
  int I,
      J,
      NS,
      ISTARTM,
      ISTOPM,
      SHEIGHT,
      SWIDTH,
      K,
      NP,
      ISTARTB,
      ISTOPB,
      ISHIFT,
      NBLOCK,
      NPOS;
  double SAFMIN, SAFMAX, SCALE;
  Complex TEMP2, TEMP3;
  final C = Box(0.0);
  final S = Box(Complex.zero), TEMP = Box(Complex.zero);

  INFO.value = 0;
  if (NBLOCK_DESIRED < NSHIFTS + 1) {
    INFO.value = -8;
  }
  if (LWORK == -1) {
    // workspace query, quick return;
    WORK[1] = (N * NBLOCK_DESIRED).toComplex();
    return;
  } else if (LWORK < N * NBLOCK_DESIRED) {
    INFO.value = -25;
  }

  if (INFO.value != 0) {
    xerbla('ZLAQZ3', -INFO.value);
    return;
  }

  // Executable statements

  // Get machine constants
  SAFMIN = dlamch('SAFE MINIMUM');
  SAFMAX = ONE / SAFMIN;

  if (ILO >= IHI) {
    return;
  }

  if (ILSCHUR) {
    ISTARTM = 1;
    ISTOPM = N;
  } else {
    ISTARTM = ILO;
    ISTOPM = IHI;
  }

  NS = NSHIFTS;
  NPOS = max(NBLOCK_DESIRED - NS, 1);

  // The following block introduces the shifts and chases
  // them down one by one just enough to make space for
  // the other shifts. The near-the-diagonal block is
  // of size (ns+1) x ns.

  zlaset('FULL', NS + 1, NS + 1, Complex.zero, Complex.one, QC, LDQC);
  zlaset('FULL', NS, NS, Complex.zero, Complex.one, ZC, LDZC);

  for (I = 1; I <= NS; I++) {
    // Introduce the shift
    SCALE = sqrt((ALPHA[I]).abs()) * sqrt((BETA[I]).abs());
    if (SCALE >= SAFMIN && SCALE <= SAFMAX) {
      ALPHA[I] /= SCALE.toComplex();
      BETA[I] /= SCALE.toComplex();
    }

    TEMP2 = BETA[I] * A[ILO][ILO] - ALPHA[I] * B[ILO][ILO];
    TEMP3 = BETA[I] * A[ILO + 1][ILO];
    if ((TEMP2).abs() > SAFMAX || (TEMP3).abs() > SAFMAX) {
      TEMP2 = Complex.one;
      TEMP3 = Complex.zero;
    }

    zlartg(TEMP2, TEMP3, C, S, TEMP);
    zrot(NS, A(ILO, ILO).asArray(), LDA, A(ILO + 1, ILO).asArray(), LDA,
        C.value, S.value);
    zrot(NS, B(ILO, ILO).asArray(), LDB, B(ILO + 1, ILO).asArray(), LDB,
        C.value, S.value);
    zrot(NS + 1, QC(1, 1).asArray(), 1, QC(1, 2).asArray(), 1, C.value,
        S.value.conjugate());

    // Chase the shift down
    for (J = 1; J <= NS - I; J++) {
      zlaqz1(true, true, J, 1, NS, IHI - ILO + 1, A(ILO, ILO), LDA, B(ILO, ILO),
          LDB, NS + 1, 1, QC, LDQC, NS, 1, ZC, LDZC);
    }
  }

  // Update the rest of the pencil

  // Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm)
  // from the left with Qc(1:ns+1,1:ns+1)'
  SHEIGHT = NS + 1;
  SWIDTH = ISTOPM - (ILO + NS) + 1;
  if (SWIDTH > 0) {
    zgemm('C', 'N', SHEIGHT, SWIDTH, SHEIGHT, Complex.one, QC, LDQC,
        A(ILO, ILO + NS), LDA, Complex.zero, WORK.asMatrix(SHEIGHT), SHEIGHT);
    zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        A(ILO, ILO + NS), LDA);
    zgemm('C', 'N', SHEIGHT, SWIDTH, SHEIGHT, Complex.one, QC, LDQC,
        B(ILO, ILO + NS), LDB, Complex.zero, WORK.asMatrix(SHEIGHT), SHEIGHT);
    zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        B(ILO, ILO + NS), LDB);
  }
  if (ILQ) {
    zgemm('N', 'N', N, SHEIGHT, SHEIGHT, Complex.one, Q(1, ILO), LDQ, QC, LDQC,
        Complex.zero, WORK.asMatrix(N), N);
    zlacpy('ALL', N, SHEIGHT, WORK.asMatrix(N), N, Q(1, ILO), LDQ);
  }

  // Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1)
  // from the right with Zc(1:ns,1:ns)
  SHEIGHT = ILO - 1 - ISTARTM + 1;
  SWIDTH = NS;
  if (SHEIGHT > 0) {
    zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, Complex.one, A(ISTARTM, ILO), LDA,
        ZC, LDZC, Complex.zero, WORK.asMatrix(SHEIGHT), SHEIGHT);
    zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        A(ISTARTM, ILO), LDA);
    zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, Complex.one, B(ISTARTM, ILO), LDB,
        ZC, LDZC, Complex.zero, WORK.asMatrix(SHEIGHT), SHEIGHT);
    zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        B(ISTARTM, ILO), LDB);
  }
  if (ILZ) {
    zgemm('N', 'N', N, SWIDTH, SWIDTH, Complex.one, Z(1, ILO), LDZ, ZC, LDZC,
        Complex.zero, WORK.asMatrix(N), N);
    zlacpy('ALL', N, SWIDTH, WORK.asMatrix(N), N, Z(1, ILO), LDZ);
  }

  // The following block chases the shifts down to the bottom
  // right block. If possible, a shift is moved down npos
  // positions at a time

  K = ILO;
  while (K < IHI - NS) {
    NP = min(IHI - NS - K, NPOS);
    // Size of the near-the-diagonal block
    NBLOCK = NS + NP;
    // istartb points to the first row we will be updating
    ISTARTB = K + 1;
    // istopb points to the last column we will be updating
    ISTOPB = K + NBLOCK - 1;

    zlaset('FULL', NS + NP, NS + NP, Complex.zero, Complex.one, QC, LDQC);
    zlaset('FULL', NS + NP, NS + NP, Complex.zero, Complex.one, ZC, LDZC);

    // Near the diagonal shift chase
    for (I = NS - 1; I >= 0; I--) {
      for (J = 0; J <= NP - 1; J++) {
        // Move down the block with index k+i+j, updating
        // the (ns+np x ns+np) block:
        // (k:k+ns+np,k:k+ns+np-1)
        zlaqz1(true, true, K + I + J, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB,
            NBLOCK, K + 1, QC, LDQC, NBLOCK, K, ZC, LDZC);
      }
    }

    // Update rest of the pencil

    // Update A(k+1:k+ns+np, k+ns+np:istopm) and
    // B(k+1:k+ns+np, k+ns+np:istopm)
    // from the left with Qc(1:ns+np,1:ns+np)'
    SHEIGHT = NS + NP;
    SWIDTH = ISTOPM - (K + NS + NP) + 1;
    if (SWIDTH > 0) {
      zgemm(
          'C.value',
          'N',
          SHEIGHT,
          SWIDTH,
          SHEIGHT,
          Complex.one,
          QC,
          LDQC,
          A(K + 1, K + NS + NP),
          LDA,
          Complex.zero,
          WORK.asMatrix(SHEIGHT),
          SHEIGHT);
      zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
          A(K + 1, K + NS + NP), LDA);
      zgemm(
          'C.value',
          'N',
          SHEIGHT,
          SWIDTH,
          SHEIGHT,
          Complex.one,
          QC,
          LDQC,
          B(K + 1, K + NS + NP),
          LDB,
          Complex.zero,
          WORK.asMatrix(SHEIGHT),
          SHEIGHT);
      zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
          B(K + 1, K + NS + NP), LDB);
    }
    if (ILQ) {
      zgemm('N', 'N', N, NBLOCK, NBLOCK, Complex.one, Q(1, K + 1), LDQ, QC,
          LDQC, Complex.zero, WORK.asMatrix(N), N);
      zlacpy('ALL', N, NBLOCK, WORK.asMatrix(N), N, Q(1, K + 1), LDQ);
    }

    // Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1)
    // from the right with Zc(1:ns+np,1:ns+np)
    SHEIGHT = K - ISTARTM + 1;
    SWIDTH = NBLOCK;
    if (SHEIGHT > 0) {
      zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, Complex.one, A(ISTARTM, K), LDA,
          ZC, LDZC, Complex.zero, WORK.asMatrix(SHEIGHT), SHEIGHT);
      zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
          A(ISTARTM, K), LDA);
      zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, Complex.one, B(ISTARTM, K), LDB,
          ZC, LDZC, Complex.zero, WORK.asMatrix(SHEIGHT), SHEIGHT);
      zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
          B(ISTARTM, K), LDB);
    }
    if (ILZ) {
      zgemm('N', 'N', N, NBLOCK, NBLOCK, Complex.one, Z(1, K), LDZ, ZC, LDZC,
          Complex.zero, WORK.asMatrix(N), N);
      zlacpy('ALL', N, NBLOCK, WORK.asMatrix(N), N, Z(1, K), LDZ);
    }

    K += NP;
  }

  // The following block removes the shifts from the bottom right corner
  // one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi).

  zlaset('FULL', NS, NS, Complex.zero, Complex.one, QC, LDQC);
  zlaset('FULL', NS + 1, NS + 1, Complex.zero, Complex.one, ZC, LDZC);

  // istartb points to the first row we will be updating
  ISTARTB = IHI - NS + 1;
  // istopb points to the last column we will be updating
  ISTOPB = IHI;

  for (I = 1; I <= NS; I++) {
    // Chase the shift down to the bottom right corner
    for (ISHIFT = IHI - I; ISHIFT <= IHI - 1; ISHIFT++) {
      zlaqz1(true, true, ISHIFT, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NS,
          IHI - NS + 1, QC, LDQC, NS + 1, IHI - NS, ZC, LDZC);
    }
  }

  // Update rest of the pencil

  // Update A(ihi-ns+1:ihi, ihi+1:istopm)
  // from the left with Qc(1:ns,1:ns)'
  SHEIGHT = NS;
  SWIDTH = ISTOPM - (IHI + 1) + 1;
  if (SWIDTH > 0) {
    zgemm(
        'C.value',
        'N',
        SHEIGHT,
        SWIDTH,
        SHEIGHT,
        Complex.one,
        QC,
        LDQC,
        A(IHI - NS + 1, IHI + 1),
        LDA,
        Complex.zero,
        WORK.asMatrix(SHEIGHT),
        SHEIGHT);
    zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        A(IHI - NS + 1, IHI + 1), LDA);
    zgemm(
        'C.value',
        'N',
        SHEIGHT,
        SWIDTH,
        SHEIGHT,
        Complex.one,
        QC,
        LDQC,
        B(IHI - NS + 1, IHI + 1),
        LDB,
        Complex.zero,
        WORK.asMatrix(SHEIGHT),
        SHEIGHT);
    zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        B(IHI - NS + 1, IHI + 1), LDB);
  }
  if (ILQ) {
    zgemm('N', 'N', N, NS, NS, Complex.one, Q(1, IHI - NS + 1), LDQ, QC, LDQC,
        Complex.zero, WORK.asMatrix(N), N);
    zlacpy('ALL', N, NS, WORK.asMatrix(N), N, Q(1, IHI - NS + 1), LDQ);
  }

  // Update A(istartm:ihi-ns,ihi-ns:ihi)
  // from the right with Zc(1:ns+1,1:ns+1)
  SHEIGHT = IHI - NS - ISTARTM + 1;
  SWIDTH = NS + 1;
  if (SHEIGHT > 0) {
    zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, Complex.one, A(ISTARTM, IHI - NS),
        LDA, ZC, LDZC, Complex.zero, WORK.asMatrix(SHEIGHT), SHEIGHT);
    zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        A(ISTARTM, IHI - NS), LDA);
    zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, Complex.one, B(ISTARTM, IHI - NS),
        LDB, ZC, LDZC, Complex.zero, WORK.asMatrix(SHEIGHT), SHEIGHT);
    zlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        B(ISTARTM, IHI - NS), LDB);
  }
  if (ILZ) {
    zgemm('N', 'N', N, NS + 1, NS + 1, Complex.one, Z(1, IHI - NS), LDZ, ZC,
        LDZC, Complex.zero, WORK.asMatrix(N), N);
    zlacpy('ALL', N, NS + 1, WORK.asMatrix(N), N, Z(1, IHI - NS), LDZ);
  }
}
