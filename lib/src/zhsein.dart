import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlaein.dart';
import 'package:lapack/src/zlanhs.dart';

void zhsein(
  final String SIDE,
  final String EIGSRC,
  final String INITV,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<Complex> H_,
  final int LDH,
  final Array<Complex> W_,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final int MM,
  final Box<int> M,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<int> IFAILL_,
  final Array<int> IFAILR_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SELECT = SELECT_.having();
  final W = W_.having();
  final H = H_.having(ld: LDH);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IFAILL = IFAILL_.having();
  final IFAILR = IFAILR_.having();
  const RZERO = 0.0;
  bool BOTHV, FROMQR, LEFTV, NOINIT, RIGHTV;
  int I, K, KL, KLN, KR, KS, LDWORK;
  double EPS3 = 0, HNORM, SMLNUM, ULP, UNFL;
  Complex WK = Complex.zero;
  final IINFO = Box(0);

  double CABS1(Complex CDUM) => CDUM.toDouble().abs() + CDUM.imaginary.abs();

  // Decode and test the input parameters.

  BOTHV = lsame(SIDE, 'B');
  RIGHTV = lsame(SIDE, 'R') || BOTHV;
  LEFTV = lsame(SIDE, 'L') || BOTHV;

  FROMQR = lsame(EIGSRC, 'Q');

  NOINIT = lsame(INITV, 'N');

  // Set M.value to the number of columns required to store the selected
  // eigenvectors.

  M.value = 0;
  for (K = 1; K <= N; K++) {
    // 10
    if (SELECT[K]) M.value = M.value + 1;
  } // 10

  INFO.value = 0;
  if (!RIGHTV && !LEFTV) {
    INFO.value = -1;
  } else if (!FROMQR && !lsame(EIGSRC, 'N')) {
    INFO.value = -2;
  } else if (!NOINIT && !lsame(INITV, 'U')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDH < max(1, N)) {
    INFO.value = -7;
  } else if (LDVL < 1 || (LEFTV && LDVL < N)) {
    INFO.value = -10;
  } else if (LDVR < 1 || (RIGHTV && LDVR < N)) {
    INFO.value = -12;
  } else if (MM < M.value) {
    INFO.value = -13;
  }
  if (INFO.value != 0) {
    xerbla('ZHSEIN', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  // Set machine-dependent constants.

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');
  SMLNUM = UNFL * (N / ULP);

  LDWORK = N;

  KL = 1;
  KLN = 0;
  if (FROMQR) {
    KR = 0;
  } else {
    KR = N;
  }
  KS = 1;

  for (K = 1; K <= N; K++) {
    // 100
    if (SELECT[K]) {
      // Compute eigenvector(s) corresponding to W(K).

      if (FROMQR) {
        // If affiliation of eigenvalues is known, check whether
        // the matrix splits.

        // Determine KL and KR such that 1 <= KL <= K <= KR <= N
        // and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or
        // KR = N).

        // Then inverse iteration can be performed with the
        // submatrix H(KL:N,KL:N) for a left eigenvector, and with
        // the submatrix H(1:KR,1:KR) for a right eigenvector.

        for (I = K; I >= KL + 1; I--) {
          // 20
          if (H[I][I - 1] == Complex.zero) break;
        } // 20
        KL = I;
        if (K > KR) {
          for (I = K; I <= N - 1; I++) {
            // 40
            if (H[I + 1][I] == Complex.zero) break;
          } // 40
          KR = I;
        }
      }

      if (KL != KLN) {
        KLN = KL;

        // Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it
        // has not ben computed before.

        HNORM = zlanhs('I', KR - KL + 1, H(KL, KL), LDH, RWORK);
        if (disnan(HNORM)) {
          INFO.value = -6;
          return;
        } else if (HNORM > RZERO) {
          EPS3 = HNORM * ULP;
        } else {
          EPS3 = SMLNUM;
        }
      }

      // Perturb eigenvalue if it is close to any previous
      // selected eigenvalues affiliated to the submatrix
      // H(KL:KR,KL:KR). Close roots are modified by EPS3.

      WK = W[K];
      repeat:
      while (true) {
        for (I = K - 1; I >= KL; I--) {
          // 70
          if (SELECT[I] && CABS1(W[I] - WK) < EPS3) {
            WK = WK + EPS3.toComplex();
            continue repeat;
          }
        } // 70
        break;
      }
      W[K] = WK;

      if (LEFTV) {
        // Compute left eigenvector.

        zlaein(
            false,
            NOINIT,
            N - KL + 1,
            H(KL, KL),
            LDH,
            WK,
            VL(KL, KS).asArray(),
            WORK.asMatrix(),
            LDWORK,
            RWORK,
            EPS3,
            SMLNUM,
            IINFO);
        if (IINFO.value > 0) {
          INFO.value = INFO.value + 1;
          IFAILL[KS] = K;
        } else {
          IFAILL[KS] = 0;
        }
        for (I = 1; I <= KL - 1; I++) {
          // 80
          VL[I][KS] = Complex.zero;
        } // 80
      }
      if (RIGHTV) {
        // Compute right eigenvector.

        zlaein(true, NOINIT, KR, H, LDH, WK, VR(1, KS).asArray(),
            WORK.asMatrix(), LDWORK, RWORK, EPS3, SMLNUM, IINFO);
        if (IINFO.value > 0) {
          INFO.value = INFO.value + 1;
          IFAILR[KS] = K;
        } else {
          IFAILR[KS] = 0;
        }
        for (I = KR + 1; I <= N; I++) {
          // 90
          VR[I][KS] = Complex.zero;
        } // 90
      }
      KS = KS + 1;
    }
  } // 100
}
