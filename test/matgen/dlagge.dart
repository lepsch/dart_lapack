import 'dart:math';

import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlagge(
  final int M,
  final int N,
  final int KL,
  final int KU,
  final Array<double> D_,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> ISEED_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final A = A_.having(ld: LDA);
  final ISEED = ISEED_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double TAU, WA, WB, WN;

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KL < 0 || KL > M - 1) {
    INFO.value = -3;
  } else if (KU < 0 || KU > N - 1) {
    INFO.value = -4;
  } else if (LDA < max(1, M)) {
    INFO.value = -7;
  }
  if (INFO.value < 0) {
    xerbla('DLAGGE', -INFO.value);
    return;
  }

  // initialize A to diagonal matrix

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      A[I][J] = ZERO;
    }
  }
  for (I = 1; I <= min(M, N); I++) {
    A[I][I] = D[I];
  }

  // Quick exit if the user wants a diagonal matrix

  if ((KL == 0) && (KU == 0)) return;

  // pre- and post-multiply A by random orthogonal matrices

  for (I = min(M, N); I >= 1; I--) {
    if (I < M) {
      // generate random reflection

      dlarnv(3, ISEED, M - I + 1, WORK);
      WN = dnrm2(M - I + 1, WORK, 1);
      WA = sign(WN, WORK[1]).toDouble();
      if (WN == ZERO) {
        TAU = ZERO;
      } else {
        WB = WORK[1] + WA;
        dscal(M - I, ONE / WB, WORK(2), 1);
        WORK[1] = ONE;
        TAU = WB / WA;
      }

      // multiply A(i:m,i:n) by random reflection from the left

      dgemv('Transpose', M - I + 1, N - I + 1, ONE, A(I, I), LDA, WORK, 1, ZERO,
          WORK(M + 1), 1);
      dger(M - I + 1, N - I + 1, -TAU, WORK, 1, WORK(M + 1), 1, A(I, I), LDA);
    }
    if (I < N) {
      // generate random reflection

      dlarnv(3, ISEED, N - I + 1, WORK);
      WN = dnrm2(N - I + 1, WORK, 1);
      WA = sign(WN, WORK[1]).toDouble();
      if (WN == ZERO) {
        TAU = ZERO;
      } else {
        WB = WORK[1] + WA;
        dscal(N - I, ONE / WB, WORK(2), 1);
        WORK[1] = ONE;
        TAU = WB / WA;
      }

      // multiply A(i:m,i:n) by random reflection from the right

      dgemv('No transpose', M - I + 1, N - I + 1, ONE, A(I, I), LDA, WORK, 1,
          ZERO, WORK(N + 1), 1);
      dger(M - I + 1, N - I + 1, -TAU, WORK(N + 1), 1, WORK, 1, A(I, I), LDA);
    }
  }

  // Reduce number of subdiagonals to KL and number of superdiagonals
  // to KU

  for (I = 1; I <= max(M - 1 - KL, N - 1 - KU); I++) {
    if (KL <= KU) {
      // annihilate subdiagonal elements first (necessary if KL = 0)

      if (I <= min(M - 1 - KL, N)) {
        // generate reflection to annihilate A(kl+i+1:m,i)

        WN = dnrm2(M - KL - I + 1, A(KL + I, I).asArray(), 1);
        WA = sign(WN, A[KL + I][I]).toDouble();
        if (WN == ZERO) {
          TAU = ZERO;
        } else {
          WB = A[KL + I][I] + WA;
          dscal(M - KL - I, ONE / WB, A(KL + I + 1, I).asArray(), 1);
          A[KL + I][I] = ONE;
          TAU = WB / WA;
        }

        // apply reflection to A(kl+i:m,i+1:n) from the left

        dgemv('Transpose', M - KL - I + 1, N - I, ONE, A(KL + I, I + 1), LDA,
            A(KL + I, I).asArray(), 1, ZERO, WORK, 1);
        dger(M - KL - I + 1, N - I, -TAU, A(KL + I, I).asArray(), 1, WORK, 1,
            A(KL + I, I + 1), LDA);
        A[KL + I][I] = -WA;
      }

      if (I <= min(N - 1 - KU, M)) {
        // generate reflection to annihilate A(i,ku+i+1:n)

        WN = dnrm2(N - KU - I + 1, A(I, KU + I).asArray(), LDA);
        WA = sign(WN, A[I][KU + I]).toDouble();
        if (WN == ZERO) {
          TAU = ZERO;
        } else {
          WB = A[I][KU + I] + WA;
          dscal(N - KU - I, ONE / WB, A(I, KU + I + 1).asArray(), LDA);
          A[I][KU + I] = ONE;
          TAU = WB / WA;
        }

        // apply reflection to A(i+1:m,ku+i:n) from the right

        dgemv('No transpose', M - I, N - KU - I + 1, ONE, A(I + 1, KU + I), LDA,
            A(I, KU + I).asArray(), LDA, ZERO, WORK, 1);
        dger(M - I, N - KU - I + 1, -TAU, WORK, 1, A(I, KU + I).asArray(), LDA,
            A(I + 1, KU + I), LDA);
        A[I][KU + I] = -WA;
      }
    } else {
      // annihilate superdiagonal elements first (necessary if
      // KU = 0)

      if (I <= min(N - 1 - KU, M)) {
        // generate reflection to annihilate A(i,ku+i+1:n)

        WN = dnrm2(N - KU - I + 1, A(I, KU + I).asArray(), LDA);
        WA = sign(WN, A[I][KU + I]).toDouble();
        if (WN == ZERO) {
          TAU = ZERO;
        } else {
          WB = A[I][KU + I] + WA;
          dscal(N - KU - I, ONE / WB, A(I, KU + I + 1).asArray(), LDA);
          A[I][KU + I] = ONE;
          TAU = WB / WA;
        }

        // apply reflection to A(i+1:m,ku+i:n) from the right

        dgemv('No transpose', M - I, N - KU - I + 1, ONE, A(I + 1, KU + I), LDA,
            A(I, KU + I).asArray(), LDA, ZERO, WORK, 1);
        dger(M - I, N - KU - I + 1, -TAU, WORK, 1, A(I, KU + I).asArray(), LDA,
            A(I + 1, KU + I), LDA);
        A[I][KU + I] = -WA;
      }

      if (I <= min(M - 1 - KL, N)) {
        // generate reflection to annihilate A(kl+i+1:m,i)

        WN = dnrm2(M - KL - I + 1, A(KL + I, I).asArray(), 1);
        WA = sign(WN, A[KL + I][I]).toDouble();
        if (WN == ZERO) {
          TAU = ZERO;
        } else {
          WB = A[KL + I][I] + WA;
          dscal(M - KL - I, ONE / WB, A(KL + I + 1, I).asArray(), 1);
          A[KL + I][I] = ONE;
          TAU = WB / WA;
        }

        // apply reflection to A(kl+i:m,i+1:n) from the left

        dgemv('Transpose', M - KL - I + 1, N - I, ONE, A(KL + I, I + 1), LDA,
            A(KL + I, I).asArray(), 1, ZERO, WORK, 1);
        dger(M - KL - I + 1, N - I, -TAU, A(KL + I, I).asArray(), 1, WORK, 1,
            A(KL + I, I + 1), LDA);
        A[KL + I][I] = -WA;
      }
    }

    if (I <= N) {
      for (J = KL + I + 1; J <= M; J++) {
        A[J][I] = ZERO;
      }
    }

    if (I <= M) {
      for (J = KU + I + 1; J <= N; J++) {
        A[I][J] = ZERO;
      }
    }
  }
}
