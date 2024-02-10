import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlaein.dart';
import 'package:lapack/src/dlanhs.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dhsein(
  final String SIDE,
  final String EIGSRC,
  final String INITV,
  final Array<bool> SELECT,
  final int N,
  final Matrix<double> H,
  final int LDH,
  final Array<double> WR,
  final Array<double> WI,
  final Matrix<double> VL,
  final int LDVL,
  final Matrix<double> VR,
  final int LDVR,
  final int MM,
  final Box<int> M,
  final Array<double> WORK,
  final Array<int> IFAILL,
  final Array<int> IFAILR,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool BOTHV, FROMQR, LEFTV, NOINIT, PAIR, RIGHTV;
  int I, K, KL, KLN, KR, KSI, KSR, LDWORK;
  double BIGNUM, EPS3 = 0, HNORM, SMLNUM, ULP, UNFL, WKI, WKR;
  final IINFO = Box(0);

  // Decode and test the input parameters.

  BOTHV = lsame(SIDE, 'B');
  RIGHTV = lsame(SIDE, 'R') || BOTHV;
  LEFTV = lsame(SIDE, 'L') || BOTHV;

  FROMQR = lsame(EIGSRC, 'Q');

  NOINIT = lsame(INITV, 'N');

  // Set M to the number of columns required to store the selected
  // eigenvectors, and standardize the array SELECT.

  M.value = 0;
  PAIR = false;
  for (K = 1; K <= N; K++) {
    if (PAIR) {
      PAIR = false;
      SELECT[K] = false;
    } else {
      if (WI[K] == ZERO) {
        if (SELECT[K]) M.value = M.value + 1;
      } else {
        PAIR = true;
        if (SELECT[K] || SELECT[K + 1]) {
          SELECT[K] = true;
          M.value = M.value + 2;
        }
      }
    }
  }

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
    INFO.value = -11;
  } else if (LDVR < 1 || (RIGHTV && LDVR < N)) {
    INFO.value = -13;
  } else if (MM < M.value) {
    INFO.value = -14;
  }
  if (INFO.value != 0) {
    xerbla('DHSEIN', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  // Set machine-dependent constants.

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');
  SMLNUM = UNFL * (N / ULP);
  BIGNUM = (ONE - ULP) / SMLNUM;

  LDWORK = N + 1;

  KL = 1;
  KLN = 0;
  if (FROMQR) {
    KR = 0;
  } else {
    KR = N;
  }
  KSR = 1;

  for (K = 1; K <= N; K++) {
    if (SELECT[K]) {
      // Compute eigenvector(s) corresponding to W(K).

      if (FROMQR) {
        // If affiliation of eigenvalues is known, check whether
        // the matrix splits.

        // Determine KL and KR such that 1 <= KL <= K <= KR <= N
        // and H[KL][KL-1] and H[KR+1][KR] are zero (or KL = 1 or
        // KR = N).

        // Then inverse iteration can be performed with the
        // submatrix H[KL:N][KL:N] for a left eigenvector, and with
        // the submatrix H[1:KR][1:KR] for a right eigenvector.

        for (I = K; I >= KL + 1; I--) {
          if (H[I][I - 1] == ZERO) break;
        }
        //  }
        KL = I;
        if (K > KR) {
          for (I = K; I <= N - 1; I++) {
            if (H[I + 1][I] == ZERO) break;
          }
          // }
          KR = I;
        }
      }

      if (KL != KLN) {
        KLN = KL;

        // Compute infinity-norm of submatrix H[KL:KR][KL:KR] if it
        // has not ben computed before.

        HNORM = dlanhs('I', KR - KL + 1, H(KL, KL), LDH, WORK);
        if (disnan(HNORM)) {
          INFO.value = -6;
          return;
        } else if (HNORM > ZERO) {
          EPS3 = HNORM * ULP;
        } else {
          EPS3 = SMLNUM;
        }
      }

      // Perturb eigenvalue if it is close to any previous
      // selected eigenvalues affiliated to the submatrix
      // H[KL:KR][KL:KR]. Close roots are modified by EPS3.

      WKR = WR[K];
      WKI = WI[K];
      // }
      restart:
      while (true) {
        for (I = K - 1; I >= KL; I--) {
          if (SELECT[I] && (WR[I] - WKR).abs() + (WI[I] - WKI).abs() < EPS3) {
            WKR = WKR + EPS3;
            continue restart;
          }
        }
        break;
      }
      WR[K] = WKR;

      PAIR = WKI != ZERO;
      if (PAIR) {
        KSI = KSR + 1;
      } else {
        KSI = KSR;
      }
      if (LEFTV) {
        // Compute left eigenvector.

        dlaein(
            false,
            NOINIT,
            N - KL + 1,
            H(KL, KL),
            LDH,
            WKR,
            WKI,
            VL(KL, KSR).asArray(),
            VL(KL, KSI).asArray(),
            WORK.asMatrix(LDWORK),
            LDWORK,
            WORK(N * N + N + 1),
            EPS3,
            SMLNUM,
            BIGNUM,
            IINFO);
        if (IINFO.value > 0) {
          if (PAIR) {
            INFO.value = INFO.value + 2;
          } else {
            INFO.value = INFO.value + 1;
          }
          IFAILL[KSR] = K;
          IFAILL[KSI] = K;
        } else {
          IFAILL[KSR] = 0;
          IFAILL[KSI] = 0;
        }
        for (I = 1; I <= KL - 1; I++) {
          VL[I][KSR] = ZERO;
        }
        if (PAIR) {
          for (I = 1; I <= KL - 1; I++) {
            VL[I][KSI] = ZERO;
          }
        }
      }
      if (RIGHTV) {
        // Compute right eigenvector.

        dlaein(
            true,
            NOINIT,
            KR,
            H,
            LDH,
            WKR,
            WKI,
            VR(1, KSR).asArray(),
            VR(1, KSI).asArray(),
            WORK.asMatrix(LDWORK),
            LDWORK,
            WORK(N * N + N + 1),
            EPS3,
            SMLNUM,
            BIGNUM,
            IINFO);
        if (IINFO.value > 0) {
          if (PAIR) {
            INFO.value = INFO.value + 2;
          } else {
            INFO.value = INFO.value + 1;
          }
          IFAILR[KSR] = K;
          IFAILR[KSI] = K;
        } else {
          IFAILR[KSR] = 0;
          IFAILR[KSI] = 0;
        }
        for (I = KR + 1; I <= N; I++) {
          VR[I][KSR] = ZERO;
        }
        if (PAIR) {
          for (I = KR + 1; I <= N; I++) {
            VR[I][KSI] = ZERO;
          }
        }
      }

      if (PAIR) {
        KSR = KSR + 2;
      } else {
        KSR = KSR + 1;
      }
    }
  }
}
