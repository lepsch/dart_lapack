import 'dart:math';

import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/blas/ztrmv.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlarfg.dart';

void zlahr2(
  final int N,
  final int K,
  final int NB,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> Y_,
  final int LDY,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final TAU = TAU_.having();
  final Y = Y_.having(ld: LDY);
  int I;
  Complex EI = Complex.zero;

  // Quick return if possible

  if (N <= 1) return;

  for (I = 1; I <= NB; I++) {
    // 10
    if (I > 1) {
      // Update A(K+1:N,I)

      // Update I-th column of A - Y * V**H

      zlacgv(I - 1, A(K + I - 1, 1).asArray(), LDA);
      zgemv(
          'NO TRANSPOSE',
          N - K,
          I - 1,
          -Complex.one,
          Y(K + 1, 1),
          LDY,
          A(K + I - 1, 1).asArray(),
          LDA,
          Complex.one,
          A(K + 1, I).asArray(),
          1);
      zlacgv(I - 1, A(K + I - 1, 1).asArray(), LDA);

      // Apply I - V * T**H * V**H to this column (call it b) from the
      // left, using the last column of T as workspace

      // Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
      //          ( V2 )             ( b2 )

      // where V1 is unit lower triangular

      // w := V1**H * b1

      zcopy(I - 1, A(K + 1, I).asArray(), 1, T(1, NB).asArray(), 1);
      ztrmv('Lower', 'Conjugate transpose', 'UNIT', I - 1, A(K + 1, 1), LDA,
          T(1, NB).asArray(), 1);

      // w := w + V2**H * b2

      zgemv(
          'Conjugate transpose',
          N - K - I + 1,
          I - 1,
          Complex.one,
          A(K + I, 1),
          LDA,
          A(K + I, I).asArray(),
          1,
          Complex.one,
          T(1, NB).asArray(),
          1);

      // w := T**H * w

      ztrmv('Upper', 'Conjugate transpose', 'NON-UNIT', I - 1, T, LDT,
          T(1, NB).asArray(), 1);

      // b2 := b2 - V2*w

      zgemv('NO TRANSPOSE', N - K - I + 1, I - 1, -Complex.one, A(K + I, 1),
          LDA, T(1, NB).asArray(), 1, Complex.one, A(K + I, I).asArray(), 1);

      // b1 := b1 - V1*w

      ztrmv('Lower', 'NO TRANSPOSE', 'UNIT', I - 1, A(K + 1, 1), LDA,
          T(1, NB).asArray(), 1);
      zaxpy(
          I - 1, -Complex.one, T(1, NB).asArray(), 1, A(K + 1, I).asArray(), 1);

      A[K + I - 1][I - 1] = EI;
    }

    // Generate the elementary reflector H(I) to annihilate
    // A(K+I+1:N,I)

    zlarfg(N - K - I + 1, A(K + I, I), A(min(K + I + 1, N), I).asArray(), 1,
        TAU(I));
    EI = A[K + I][I];
    A[K + I][I] = Complex.one;

    // Compute  Y(K+1:N,I)

    zgemv('NO TRANSPOSE', N - K, N - K - I + 1, Complex.one, A(K + 1, I + 1),
        LDA, A(K + I, I).asArray(), 1, Complex.zero, Y(K + 1, I).asArray(), 1);
    zgemv('Conjugate transpose', N - K - I + 1, I - 1, Complex.one, A(K + I, 1),
        LDA, A(K + I, I).asArray(), 1, Complex.zero, T(1, I).asArray(), 1);
    zgemv('NO TRANSPOSE', N - K, I - 1, -Complex.one, Y(K + 1, 1), LDY,
        T(1, I).asArray(), 1, Complex.one, Y(K + 1, I).asArray(), 1);
    zscal(N - K, TAU[I], Y(K + 1, I).asArray(), 1);

    // Compute T(1:I,I)

    zscal(I - 1, -TAU[I], T(1, I).asArray(), 1);
    ztrmv('Upper', 'No Transpose', 'NON-UNIT', I - 1, T, LDT, T(1, I).asArray(),
        1);
    T[I][I] = TAU[I];
  } // 10
  A[K + NB][NB] = EI;

  // Compute Y(1:K,1:NB)

  zlacpy('ALL', K, NB, A(1, 2), LDA, Y, LDY);
  ztrmm('RIGHT', 'Lower', 'NO TRANSPOSE', 'UNIT', K, NB, Complex.one,
      A(K + 1, 1), LDA, Y, LDY);
  if (N > K + NB) {
    zgemm('NO TRANSPOSE', 'NO TRANSPOSE', K, NB, N - K - NB, Complex.one,
        A(1, 2 + NB), LDA, A(K + 1 + NB, 1), LDA, Complex.one, Y, LDY);
  }
  ztrmm('RIGHT', 'Upper', 'NO TRANSPOSE', 'NON-UNIT', K, NB, Complex.one, T,
      LDT, Y, LDY);
}
