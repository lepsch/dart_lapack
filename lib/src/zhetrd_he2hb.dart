import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zhemm.dart';
import 'package:lapack/src/blas/zher2k.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/qr/ll/zgeqrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgelqf.dart';
import 'package:lapack/src/zlarft.dart';
import 'package:lapack/src/zlaset.dart';

void zhetrd_he2hb(
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final AB = AB_.having(ld: LDAB);
  final TAU = TAU_.having();
  final WORK = WORK_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const RONE = 1.0, HALF = Complex(0.5, 0.0);
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
    LWMIN = ilaenv2stage(4, 'ZHETRD_HE2HB', '', N, KD, -1, -1);
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
    xerbla('ZHETRD_HE2HB', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = LWMIN.toComplex();
    return;
  }

  // Quick return if possible
  // Copy the upper/lower portion of A into AB

  if (N <= KD + 1) {
    if (UPPER) {
      for (I = 1; I <= N; I++) {
        // 100
        LK = min(KD + 1, I);
        zcopy(LK, A(I - LK + 1, I).asArray(), 1,
            AB(KD + 1 - LK + 1, I).asArray(), 1);
      } // 100
    } else {
      for (I = 1; I <= N; I++) {
        // 110
        LK = min(KD + 1, N - I + 1);
        zcopy(LK, A(I, I).asArray(), 1, AB(1, I).asArray(), 1);
      } // 110
    }
    WORK[1] = Complex.one;
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

  zlaset(
      'A', LDT, KD, Complex.zero, Complex.zero, WORK(TPOS).asMatrix(LDT), LDT);

  if (UPPER) {
    for (I = 1; KD < 0 ? I >= N - KD : I <= N - KD; I += KD) {
      // 10
      PN = N - I - KD + 1;
      PK = min(N - I - KD + 1, KD);

      // Compute the LQ factorization of the current block

      zgelqf(KD, PN, A(I, I + KD), LDA, TAU(I), WORK(S2POS), LS2, IINFO);

      // Copy the upper portion of A into AB

      for (J = I; J <= I + PK - 1; J++) {
        // 20
        LK = min(KD, N - J) + 1;
        zcopy(LK, A(J, J).asArray(), LDA, AB(KD + 1, J).asArray(), LDAB - 1);
      } // 20

      zlaset('Lower', PK, PK, Complex.zero, Complex.one, A(I, I + KD), LDA);

      // Form the matrix T

      zlarft('Forward', 'Rowwise', PN, PK, A(I, I + KD), LDA, TAU(I),
          WORK(TPOS).asMatrix(LDT), LDT);

      // Compute W:

      zgemm(
          'Conjugate',
          'No transpose',
          PK,
          PN,
          PK,
          Complex.one,
          WORK(TPOS).asMatrix(LDT),
          LDT,
          A(I, I + KD),
          LDA,
          Complex.zero,
          WORK(S2POS).asMatrix(LDS2),
          LDS2);

      zhemm(
          'Right',
          UPLO,
          PK,
          PN,
          Complex.one,
          A(I + KD, I + KD),
          LDA,
          WORK(S2POS).asMatrix(LDS2),
          LDS2,
          Complex.zero,
          WORK(WPOS).asMatrix(LDW),
          LDW);

      zgemm(
          'No transpose',
          'Conjugate',
          PK,
          PK,
          PN,
          Complex.one,
          WORK(WPOS).asMatrix(LDW),
          LDW,
          WORK(S2POS).asMatrix(LDS2),
          LDS2,
          Complex.zero,
          WORK(S1POS).asMatrix(LDS1),
          LDS1);

      zgemm(
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
          Complex.one,
          WORK(WPOS).asMatrix(LDW),
          LDW);

      // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
      // an update of the form:  A := A - V'*W - W'*V

      zher2k(UPLO, 'Conjugate', PN, PK, -Complex.one, A(I, I + KD), LDA,
          WORK(WPOS).asMatrix(LDW), LDW, RONE, A(I + KD, I + KD), LDA);
    } // 10

    // Copy the upper band to AB which is the band storage matrix

    for (J = N - KD + 1; J <= N; J++) {
      // 30
      LK = min(KD, N - J) + 1;
      zcopy(LK, A(J, J).asArray(), LDA, AB(KD + 1, J).asArray(), LDAB - 1);
    } // 30
  } else {
    // Reduce the lower triangle of A to lower band matrix

    for (I = 1; KD < 0 ? I >= N - KD : I <= N - KD; I += KD) {
      // 40
      PN = N - I - KD + 1;
      PK = min(N - I - KD + 1, KD);

      // Compute the QR factorization of the current block

      zgeqrf(PN, KD, A(I + KD, I), LDA, TAU(I), WORK(S2POS), LS2, IINFO);

      // Copy the upper portion of A into AB

      for (J = I; J <= I + PK - 1; J++) {
        // 50
        LK = min(KD, N - J) + 1;
        zcopy(LK, A(J, J).asArray(), 1, AB(1, J).asArray(), 1);
      } // 50

      zlaset('Upper', PK, PK, Complex.zero, Complex.one, A(I + KD, I), LDA);

      // Form the matrix T

      zlarft('Forward', 'Columnwise', PN, PK, A(I + KD, I), LDA, TAU(I),
          WORK(TPOS).asMatrix(LDT), LDT);

      // Compute W:

      zgemm(
          'No transpose',
          'No transpose',
          PN,
          PK,
          PK,
          Complex.one,
          A(I + KD, I),
          LDA,
          WORK(TPOS).asMatrix(LDT),
          LDT,
          Complex.zero,
          WORK(S2POS).asMatrix(LDS2),
          LDS2);

      zhemm(
          'Left',
          UPLO,
          PN,
          PK,
          Complex.one,
          A(I + KD, I + KD),
          LDA,
          WORK(S2POS).asMatrix(LDS2),
          LDS2,
          Complex.zero,
          WORK(WPOS).asMatrix(LDW),
          LDW);

      zgemm(
          'Conjugate',
          'No transpose',
          PK,
          PK,
          PN,
          Complex.one,
          WORK(S2POS).asMatrix(LDS2),
          LDS2,
          WORK(WPOS).asMatrix(LDW),
          LDW,
          Complex.zero,
          WORK(S1POS).asMatrix(LDS1),
          LDS1);

      zgemm(
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
          Complex.one,
          WORK(WPOS).asMatrix(LDW),
          LDW);

      // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
      // an update of the form:  A := A - V*W' - W*V'

      zher2k(UPLO, 'No transpose', PN, PK, -Complex.one, A(I + KD, I), LDA,
          WORK(WPOS).asMatrix(LDW), LDW, RONE, A(I + KD, I + KD), LDA);
      // ==================================================================
      // RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED
      //  DO 45 J = I, I+PK-1
      //     LK = min( KD, N-J ) + 1
      //     CALL ZCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 )
      // 45        CONTINUE
      // ==================================================================
    } // 40

    // Copy the lower band to AB which is the band storage matrix

    for (J = N - KD + 1; J <= N; J++) {
      // 60
      LK = min(KD, N - J) + 1;
      zcopy(LK, A(J, J).asArray(), 1, AB(1, J).asArray(), 1);
    } // 60
  }

  WORK[1] = LWMIN.toComplex();
}
