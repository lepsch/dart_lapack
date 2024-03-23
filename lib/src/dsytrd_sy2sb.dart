import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsymm.dart';
import 'package:lapack/src/blas/dsyr2k.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgelqf.dart';
import 'package:lapack/src/dlarft.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/dgeqrf.dart';
import 'package:lapack/src/xerbla.dart';

void dsytrd_sy2sb(
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AB = AB_.having(ld: LDAB);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const RONE = 1.0, ZERO = 0.0, ONE = 1.0, HALF = 0.5;
  bool LQUERY, UPPER;
  int I,
      J,
      LWMIN,
      PN,
      PK,
      LK,
      LDT,
      LDW,
      LDS2,
      LDS1,
      LS2,
      LS1,
      LW,
      LT,
      TPOS,
      WPOS,
      S2POS,
      S1POS;
  final IINFO = Box(0);

  // Determine the minimal workspace size required
  // and test the input parameters

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);
  if (N <= KD + 1) {
    LWMIN = 1;
  } else {
    LWMIN = ilaenv2stage(4, 'DSYTRD_SY2SB', ' ', N, KD, -1, -1);
  }

  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KD < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDAB < max(1, KD + 1)) {
    INFO.value = -7;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -10;
  }

  if (INFO.value != 0) {
    xerbla('DSYTRD_SY2SB', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = LWMIN.toDouble();
    return;
  }

  // Quick return if possible
  // Copy the upper/lower portion of A into AB

  if (N <= KD + 1) {
    if (UPPER) {
      for (I = 1; I <= N; I++) {
        LK = min(KD + 1, I);
        dcopy(LK, A(I - LK + 1, I).asArray(), 1,
            AB(KD + 1 - LK + 1, I).asArray(), 1);
      }
    } else {
      for (I = 1; I <= N; I++) {
        LK = min(KD + 1, N - I + 1);
        dcopy(LK, A(I, I).asArray(), 1, AB(1, I).asArray(), 1);
      }
    }
    WORK[1] = 1;
    return;
  }

  // Determine the pointer position for the workspace

  LDT = KD;
  LDS1 = KD;
  LT = LDT * KD;
  LW = N * KD;
  LS1 = LDS1 * KD;
  LS2 = LWMIN - LT - LW - LS1;
  // LS2 = N*max(KD,FACTOPTNB)
  TPOS = 1;
  WPOS = TPOS + LT;
  S1POS = WPOS + LW;
  S2POS = S1POS + LS1;
  if (UPPER) {
    LDW = KD;
    LDS2 = KD;
  } else {
    LDW = N;
    LDS2 = N;
  }

  // Set the workspace of the triangular matrix T to zero once such a
  // way every time T is generated the upper/lower portion will be always zero

  dlaset('A', LDT, KD, ZERO, ZERO, WORK(TPOS).asMatrix(LDT), LDT);

  if (UPPER) {
    for (I = 1; I <= N - KD; I += KD) {
      PN = N - I - KD + 1;
      PK = min(N - I - KD + 1, KD);

      // Compute the LQ factorization of the current block

      dgelqf(KD, PN, A(I, I + KD), LDA, TAU(I), WORK(S2POS), LS2, IINFO);

      // Copy the upper portion of A into AB

      for (J = I; J <= I + PK - 1; J++) {
        LK = min(KD, N - J) + 1;
        dcopy(LK, A(J, J).asArray(), LDA, AB(KD + 1, J).asArray(), LDAB - 1);
      }

      dlaset('Lower', PK, PK, ZERO, ONE, A(I, I + KD), LDA);

      // Form the matrix T

      dlarft('Forward', 'Rowwise', PN, PK, A(I, I + KD), LDA, TAU(I),
          WORK(TPOS).asMatrix(LDT), LDT);

      // Compute W:

      dgemm(
          'Conjugate',
          'No transpose',
          PK,
          PN,
          PK,
          ONE,
          WORK(TPOS).asMatrix(LDT),
          LDT,
          A(I, I + KD),
          LDA,
          ZERO,
          WORK(S2POS).asMatrix(LDS2),
          LDS2);

      dsymm(
          'Right',
          UPLO,
          PK,
          PN,
          ONE,
          A(I + KD, I + KD),
          LDA,
          WORK(S2POS).asMatrix(LDS2),
          LDS2,
          ZERO,
          WORK(WPOS).asMatrix(LDW),
          LDW);

      dgemm(
          'No transpose',
          'Conjugate',
          PK,
          PK,
          PN,
          ONE,
          WORK(WPOS).asMatrix(LDW),
          LDW,
          WORK(S2POS).asMatrix(LDS2),
          LDS2,
          ZERO,
          WORK(S1POS).asMatrix(LDS1),
          LDS1);

      dgemm(
          'No transpose',
          'No transpose',
          PK,
          PN,
          PK,
          -HALF,
          WORK(S1POS).asMatrix(LDS1),
          LDS1,
          A(I, I + KD),
          LDA,
          ONE,
          WORK(WPOS).asMatrix(LDW),
          LDW);

      // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
      // an update of the form:  A := A - V'*W - W'*V

      dsyr2k(UPLO, 'Conjugate', PN, PK, -ONE, A(I, I + KD), LDA,
          WORK(WPOS).asMatrix(LDW), LDW, RONE, A(I + KD, I + KD), LDA);
    }

    // Copy the upper band to AB which is the band storage matrix

    for (J = N - KD + 1; J <= N; J++) {
      LK = min(KD, N - J) + 1;
      dcopy(LK, A(J, J).asArray(), LDA, AB(KD + 1, J).asArray(), LDAB - 1);
    }
  } else {
    // Reduce the lower triangle of A to lower band matrix

    for (I = 1; I <= N - KD; I += KD) {
      PN = N - I - KD + 1;
      PK = min(N - I - KD + 1, KD);

      // Compute the QR factorization of the current block

      dgeqrf(PN, KD, A(I + KD, I), LDA, TAU(I), WORK(S2POS), LS2, IINFO);

      // Copy the upper portion of A into AB

      for (J = I; J <= I + PK - 1; J++) {
        LK = min(KD, N - J) + 1;
        dcopy(LK, A(J, J).asArray(), 1, AB(1, J).asArray(), 1);
      }

      dlaset('Upper', PK, PK, ZERO, ONE, A(I + KD, I), LDA);

      // Form the matrix T

      dlarft('Forward', 'Columnwise', PN, PK, A(I + KD, I), LDA, TAU(I),
          WORK(TPOS).asMatrix(LDT), LDT);

      // Compute W:

      dgemm(
          'No transpose',
          'No transpose',
          PN,
          PK,
          PK,
          ONE,
          A(I + KD, I),
          LDA,
          WORK(TPOS).asMatrix(LDT),
          LDT,
          ZERO,
          WORK(S2POS).asMatrix(LDS2),
          LDS2);

      dsymm(
          'Left',
          UPLO,
          PN,
          PK,
          ONE,
          A(I + KD, I + KD),
          LDA,
          WORK(S2POS).asMatrix(LDS2),
          LDS2,
          ZERO,
          WORK(WPOS).asMatrix(LDW),
          LDW);

      dgemm(
          'Conjugate',
          'No transpose',
          PK,
          PK,
          PN,
          ONE,
          WORK(S2POS).asMatrix(LDS2),
          LDS2,
          WORK(WPOS).asMatrix(LDW),
          LDW,
          ZERO,
          WORK(S1POS).asMatrix(LDS1),
          LDS1);

      dgemm(
          'No transpose',
          'No transpose',
          PN,
          PK,
          PK,
          -HALF,
          A(I + KD, I),
          LDA,
          WORK(S1POS).asMatrix(LDS1),
          LDS1,
          ONE,
          WORK(WPOS).asMatrix(LDW),
          LDW);

      // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
      // an update of the form:  A := A - V*W' - W*V'

      dsyr2k(UPLO, 'No transpose', PN, PK, -ONE, A(I + KD, I), LDA,
          WORK(WPOS).asMatrix(LDW), LDW, RONE, A(I + KD, I + KD), LDA);
      // ==================================================================
      // RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED
      //  DO 45 J = I, I+PK-1
      //     LK = min( KD, N-J ) + 1
      //     CALL DCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 )
// 45        CONTINUE
      // ==================================================================
    }

    // Copy the lower band to AB which is the band storage matrix

    for (J = N - KD + 1; J <= N; J++) {
      LK = min(KD, N - J) + 1;
      dcopy(LK, A(J, J).asArray(), 1, AB(1, J).asArray(), 1);
    }
  }

  WORK[1] = LWMIN.toDouble();
}
