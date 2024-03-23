import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaqz1.dart';
import 'package:lapack/src/dlaqz2.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlaqz4(
  final bool ILSCHUR,
  final bool ILQ,
  final bool ILZ,
  final int N,
  final int ILO,
  final int IHI,
  final int NSHIFTS,
  final int NBLOCK_DESIRED,
  final Array<double> SR_,
  final Array<double> SI_,
  final Array<double> SS_,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> Z_,
  final int LDZ,
  final Matrix<double> QC_,
  final int LDQC,
  final Matrix<double> ZC_,
  final int LDZC,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final SR = SR_.having();
  final SI = SI_.having();
  final SS = SS_.having();
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final QC = QC_.having(ld: LDQC);
  final ZC = ZC_.having(ld: LDZC);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Local scalars
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
  double SWAP;
  final C1 = Box(0.0),
      S1 = Box(0.0),
      C2 = Box(0.0),
      S2 = Box(0.0),
      TEMP = Box(0.0);
  final V = Array<double>(3);

  // External functions
  // EXTERNAL :: XERBLA, DGEMM, DLAQZ1, DLAQZ2, DLASET, DLARTG, DROT, DLACPY

  INFO.value = 0;
  if (NBLOCK_DESIRED < NSHIFTS + 1) {
    INFO.value = -8;
  }
  if (LWORK == -1) {
    // workspace query, quick return;
    WORK[1] = N * NBLOCK_DESIRED.toDouble();
    return;
  } else if (LWORK < N * NBLOCK_DESIRED) {
    INFO.value = -25;
  }

  if (INFO.value != 0) {
    xerbla('DLAQZ4', -INFO.value);
    return;
  }

  // Executable statements

  if (NSHIFTS < 2) {
    return;
  }

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

  // Shuffle shifts into pairs of real shifts and pairs
  // of complex conjugate shifts assuming complex
  // conjugate shifts are already adjacent to one
  // another

  for (I = 1; I <= NSHIFTS - 2; I += 2) {
    if (SI[I] != -SI[I + 1]) {
      SWAP = SR[I];
      SR[I] = SR[I + 1];
      SR[I + 1] = SR[I + 2];
      SR[I + 2] = SWAP;

      SWAP = SI[I];
      SI[I] = SI[I + 1];
      SI[I + 1] = SI[I + 2];
      SI[I + 2] = SWAP;

      SWAP = SS[I];
      SS[I] = SS[I + 1];
      SS[I + 1] = SS[I + 2];
      SS[I + 2] = SWAP;
    }
  }

  // NSHFTS is supposed to be even, but if it is odd,
  // then simply reduce it by one.  The shuffle above
  // ensures that the dropped shift is real and that
  // the remaining shifts are paired.

  NS = NSHIFTS - (NSHIFTS % 2);
  NPOS = max(NBLOCK_DESIRED - NS, 1);

  // The following block introduces the shifts and chases
  // them down one by one just enough to make space for
  // the other shifts. The near-the-diagonal block is
  // of size (ns+1) x ns.

  dlaset('FULL', NS + 1, NS + 1, ZERO, ONE, QC, LDQC);
  dlaset('FULL', NS, NS, ZERO, ONE, ZC, LDZC);

  for (I = 1; I <= NS; I += 2) {
    // Introduce the shift
    dlaqz1(A(ILO, ILO), LDA, B(ILO, ILO), LDB, SR[I], SR[I + 1], SI[I], SS[I],
        SS[I + 1], V);

    TEMP.value = V[2];
    dlartg(TEMP.value, V[3], C1, S1, V.box(2));
    dlartg(V[1], V[2], C2, S2, TEMP);
    drot(NS, A(ILO + 1, ILO).asArray(), LDA, A(ILO + 2, ILO).asArray(), LDA,
        C1.value, S1.value);
    drot(NS, A(ILO, ILO).asArray(), LDA, A(ILO + 1, ILO).asArray(), LDA,
        C2.value, S2.value);
    drot(NS, B(ILO + 1, ILO).asArray(), LDB, B(ILO + 2, ILO).asArray(), LDB,
        C1.value, S1.value);
    drot(NS, B(ILO, ILO).asArray(), LDB, B(ILO + 1, ILO).asArray(), LDB,
        C2.value, S2.value);
    drot(NS + 1, QC(1, 2).asArray(), 1, QC(1, 3).asArray(), 1, C1.value,
        S1.value);
    drot(NS + 1, QC(1, 1).asArray(), 1, QC(1, 2).asArray(), 1, C2.value,
        S2.value);

    // Chase the shift down
    for (J = 1; J <= NS - 1 - I; J++) {
      dlaqz2(true, true, J, 1, NS, IHI - ILO + 1, A(ILO, ILO), LDA, B(ILO, ILO),
          LDB, NS + 1, 1, QC, LDQC, NS, 1, ZC, LDZC);
    }
  }

  // Update the rest of the pencil

  // Update A[ilo:ilo+ns][ilo+ns:istopm] and B[ilo:ilo+ns][ilo+ns:istopm]
  // from the left with Qc[1:ns+1][1:ns+1]'
  SHEIGHT = NS + 1;
  SWIDTH = ISTOPM - (ILO + NS) + 1;
  if (SWIDTH > 0) {
    dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A(ILO, ILO + NS),
        LDA, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
    dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        A(ILO, ILO + NS), LDA);
    dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B(ILO, ILO + NS),
        LDB, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
    dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        B(ILO, ILO + NS), LDB);
  }
  if (ILQ) {
    dgemm('N', 'N', N, SHEIGHT, SHEIGHT, ONE, Q(1, ILO), LDQ, QC, LDQC, ZERO,
        WORK.asMatrix(N), N);
    dlacpy('ALL', N, SHEIGHT, WORK.asMatrix(N), N, Q(1, ILO), LDQ);
  }

  // Update A[istartm:ilo-1][ilo:ilo+ns-1] and B[istartm:ilo-1][ilo:ilo+ns-1]
  // from the right with Zc[1:ns][1:ns]
  SHEIGHT = ILO - 1 - ISTARTM + 1;
  SWIDTH = NS;
  if (SHEIGHT > 0) {
    dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A(ISTARTM, ILO), LDA, ZC,
        LDZC, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
    dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        A(ISTARTM, ILO), LDA);
    dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B(ISTARTM, ILO), LDB, ZC,
        LDZC, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
    dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        B(ISTARTM, ILO), LDB);
  }
  if (ILZ) {
    dgemm('N', 'N', N, SWIDTH, SWIDTH, ONE, Z(1, ILO), LDZ, ZC, LDZC, ZERO,
        WORK.asMatrix(N), N);
    dlacpy('ALL', N, SWIDTH, WORK.asMatrix(N), N, Z(1, ILO), LDZ);
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

    dlaset('FULL', NS + NP, NS + NP, ZERO, ONE, QC, LDQC);
    dlaset('FULL', NS + NP, NS + NP, ZERO, ONE, ZC, LDZC);

    // Near the diagonal shift chase
    for (I = NS - 1; I >= 0; I -= 2) {
      for (J = 0; J <= NP - 1; J++) {
        // Move down the block with index k+i+j-1, updating
        // the (ns+np x ns+np) block:
        // (k:k+ns+np,k:k+ns+np-1)
        dlaqz2(true, true, K + I + J - 1, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB,
            NBLOCK, K + 1, QC, LDQC, NBLOCK, K, ZC, LDZC);
      }
    }

    // Update rest of the pencil

    // Update A[k+1:k+ns+np][ k+ns+np:istopm] and
    // B[k+1:k+ns+np][ k+ns+np:istopm]
    // from the left with Qc[1:ns+np][1:ns+np]'
    SHEIGHT = NS + NP;
    SWIDTH = ISTOPM - (K + NS + NP) + 1;
    if (SWIDTH > 0) {
      dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC,
          A(K + 1, K + NS + NP), LDA, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
      dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
          A(K + 1, K + NS + NP), LDA);
      dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC,
          B(K + 1, K + NS + NP), LDB, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
      dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
          B(K + 1, K + NS + NP), LDB);
    }
    if (ILQ) {
      dgemm('N', 'N', N, NBLOCK, NBLOCK, ONE, Q(1, K + 1), LDQ, QC, LDQC, ZERO,
          WORK.asMatrix(N), N);
      dlacpy('ALL', N, NBLOCK, WORK.asMatrix(N), N, Q(1, K + 1), LDQ);
    }

    // Update A[istartm:k,k:k+ns+npos-1] and B[istartm:k,k:k+ns+npos-1]
    // from the right with Zc[1:ns+np,1:ns+np]
    SHEIGHT = K - ISTARTM + 1;
    SWIDTH = NBLOCK;
    if (SHEIGHT > 0) {
      dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A(ISTARTM, K), LDA, ZC,
          LDZC, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
      dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
          A(ISTARTM, K), LDA);
      dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B(ISTARTM, K), LDB, ZC,
          LDZC, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
      dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
          B(ISTARTM, K), LDB);
    }
    if (ILZ) {
      dgemm('N', 'N', N, NBLOCK, NBLOCK, ONE, Z(1, K), LDZ, ZC, LDZC, ZERO,
          WORK.asMatrix(N), N);
      dlacpy('ALL', N, NBLOCK, WORK.asMatrix(N), N, Z(1, K), LDZ);
    }

    K += NP;
  }

  // The following block removes the shifts from the bottom right corner
  // one by one. Updates are initially applied to A[ihi-ns+1:ihi,ihi-ns:ihi].

  dlaset('FULL', NS, NS, ZERO, ONE, QC, LDQC);
  dlaset('FULL', NS + 1, NS + 1, ZERO, ONE, ZC, LDZC);

  // istartb points to the first row we will be updating
  ISTARTB = IHI - NS + 1;
  // istopb points to the last column we will be updating
  ISTOPB = IHI;

  for (I = 1; I <= NS; I += 2) {
    // Chase the shift down to the bottom right corner
    for (ISHIFT = IHI - I - 1; ISHIFT <= IHI - 2; ISHIFT++) {
      dlaqz2(true, true, ISHIFT, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NS,
          IHI - NS + 1, QC, LDQC, NS + 1, IHI - NS, ZC, LDZC);
    }
  }

  // Update rest of the pencil

  // Update A[ihi-ns+1:ihi, ihi+1:istopm]
  // from the left with Qc[1:ns,1:ns]'
  SHEIGHT = NS;
  SWIDTH = ISTOPM - (IHI + 1) + 1;
  if (SWIDTH > 0) {
    dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC,
        A(IHI - NS + 1, IHI + 1), LDA, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
    dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        A(IHI - NS + 1, IHI + 1), LDA);
    dgemm('T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC,
        B(IHI - NS + 1, IHI + 1), LDB, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
    dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        B(IHI - NS + 1, IHI + 1), LDB);
  }
  if (ILQ) {
    dgemm('N', 'N', N, NS, NS, ONE, Q(1, IHI - NS + 1), LDQ, QC, LDQC, ZERO,
        WORK.asMatrix(N), N);
    dlacpy('ALL', N, NS, WORK.asMatrix(N), N, Q(1, IHI - NS + 1), LDQ);
  }

  // Update A[istartm:ihi-ns,ihi-ns:ihi]
  // from the right with Zc[1:ns+1,1:ns+1]
  SHEIGHT = IHI - NS - ISTARTM + 1;
  SWIDTH = NS + 1;
  if (SHEIGHT > 0) {
    dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A(ISTARTM, IHI - NS), LDA, ZC,
        LDZC, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
    dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        A(ISTARTM, IHI - NS), LDA);
    dgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B(ISTARTM, IHI - NS), LDB, ZC,
        LDZC, ZERO, WORK.asMatrix(SHEIGHT), SHEIGHT);
    dlacpy('ALL', SHEIGHT, SWIDTH, WORK.asMatrix(SHEIGHT), SHEIGHT,
        B(ISTARTM, IHI - NS), LDB);
  }
  if (ILZ) {
    dgemm('N', 'N', N, NS + 1, NS + 1, ONE, Z(1, IHI - NS), LDZ, ZC, LDZC, ZERO,
        WORK.asMatrix(N), N);
    dlacpy('ALL', N, NS + 1, WORK.asMatrix(N), N, Z(1, IHI - NS), LDZ);
  }
}
