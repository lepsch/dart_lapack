import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zgebak(
  final String JOB,
  final String SIDE,
  final int N,
  final int ILO,
  final int IHI,
  final Array<double> SCALE_,
  final int M,
  final Matrix<Complex> V_,
  final int LDV,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.dim(LDV);
  final SCALE = SCALE_.dim();
  const ONE = 1.0;
  bool LEFTV, RIGHTV;
  int I, II, K;
  double S;

  // Decode and Test the input parameters

  RIGHTV = lsame(SIDE, 'R');
  LEFTV = lsame(SIDE, 'L');

  INFO.value = 0;
  if (!lsame(JOB, 'N') &&
      !lsame(JOB, 'P') &&
      !lsame(JOB, 'S') &&
      !lsame(JOB, 'B')) {
    INFO.value = -1;
  } else if (!RIGHTV && !LEFTV) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (ILO < 1 || ILO > max(1, N)) {
    INFO.value = -4;
  } else if (IHI < min(ILO, N) || IHI > N) {
    INFO.value = -5;
  } else if (M < 0) {
    INFO.value = -7;
  } else if (LDV < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('ZGEBAK', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;
  if (M == 0) return;
  if (lsame(JOB, 'N')) return;

  if (ILO != IHI) {
    // Backward balance

    if (lsame(JOB, 'S') || lsame(JOB, 'B')) {
      if (RIGHTV) {
        for (I = ILO; I <= IHI; I++) {
          // 10
          S = SCALE[I];
          zdscal(M, S, V(I, 1).asArray(), LDV);
        } // 10
      }

      if (LEFTV) {
        for (I = ILO; I <= IHI; I++) {
          // 20
          S = ONE / SCALE[I];
          zdscal(M, S, V(I, 1).asArray(), LDV);
        } // 20
      }
    }
  }

  // Backward permutation

  // For  I = ILO-1 step -1 until 1,
  //          IHI+1 step 1 until N do --

  if (lsame(JOB, 'P') || lsame(JOB, 'B')) {
    if (RIGHTV) {
      for (II = 1; II <= N; II++) {
        // 40
        I = II;
        if (I >= ILO && I <= IHI) continue;
        if (I < ILO) I = ILO - II;
        K = SCALE[I].toInt();
        if (K == I) continue;
        zswap(M, V(I, 1).asArray(), LDV, V(K, 1).asArray(), LDV);
      } // 40
    }

    if (LEFTV) {
      for (II = 1; II <= N; II++) {
        // 50
        I = II;
        if (I >= ILO && I <= IHI) continue;
        if (I < ILO) I = ILO - II;
        K = SCALE[I].toInt();
        if (K == I) continue;
        zswap(M, V(I, 1).asArray(), LDV, V(K, 1).asArray(), LDV);
      } // 50
    }
  }
}
