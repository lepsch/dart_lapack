import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dtrmm.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';

void dlahr2(
  final int N,
  final int K,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> Y_,
  final int LDY,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final T = T_.having(ld: LDT);
  final Y = Y_.having(ld: LDY);
  const ZERO = 0.0, ONE = 1.0;
  int I;
  double EI = 0;

  // Quick return if possible

  if (N <= 1) return;

  for (I = 1; I <= NB; I++) {
    if (I > 1) {
      // Update A(K+1:N,I)

      // Update I-th column of A - Y * V**T

      dgemv('NO TRANSPOSE', N - K, I - 1, -ONE, Y(K + 1, 1), LDY,
          A(K + I - 1, 1).asArray(), LDA, ONE, A(K + 1, I).asArray(), 1);

      // Apply I - V * T**T * V**T to this column (call it b) from the
      // left, using the last column of T as workspace
      //
      // Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
      //          ( V2 )             ( b2 )
      //
      // where V1 is unit lower triangular
      //
      // w := V1**T * b1

      dcopy(I - 1, A(K + 1, I).asArray(), 1, T(1, NB).asArray(), 1);
      dtrmv('Lower', 'Transpose', 'UNIT', I - 1, A(K + 1, 1), LDA,
          T(1, NB).asArray(), 1);

      // w := w + V2**T * b2

      dgemv('Transpose', N - K - I + 1, I - 1, ONE, A(K + I, 1), LDA,
          A(K + I, I).asArray(), 1, ONE, T(1, NB).asArray(), 1);

      // w := T**T * w

      dtrmv('Upper', 'Transpose', 'NON-UNIT', I - 1, T, LDT, T(1, NB).asArray(),
          1);

      // b2 := b2 - V2*w

      dgemv('NO TRANSPOSE', N - K - I + 1, I - 1, -ONE, A(K + I, 1), LDA,
          T(1, NB).asArray(), 1, ONE, A(K + I, I).asArray(), 1);

      // b1 := b1 - V1*w

      dtrmv('Lower', 'NO TRANSPOSE', 'UNIT', I - 1, A(K + 1, 1), LDA,
          T(1, NB).asArray(), 1);
      daxpy(I - 1, -ONE, T(1, NB).asArray(), 1, A(K + 1, I).asArray(), 1);

      A[K + I - 1][I - 1] = EI;
    }

    // Generate the elementary reflector H(I) to annihilate
    // A(K+I+1:N,I)

    dlarfg(N - K - I + 1, A.box(K + I, I), A(min(K + I + 1, N), I).asArray(), 1,
        TAU.box(I));
    EI = A[K + I][I];
    A[K + I][I] = ONE;

    // Compute  Y(K+1:N,I)

    dgemv('NO TRANSPOSE', N - K, N - K - I + 1, ONE, A(K + 1, I + 1), LDA,
        A(K + I, I).asArray(), 1, ZERO, Y(K + 1, I).asArray(), 1);
    dgemv('Transpose', N - K - I + 1, I - 1, ONE, A(K + I, 1), LDA,
        A(K + I, I).asArray(), 1, ZERO, T(1, I).asArray(), 1);
    dgemv('NO TRANSPOSE', N - K, I - 1, -ONE, Y(K + 1, 1), LDY,
        T(1, I).asArray(), 1, ONE, Y(K + 1, I).asArray(), 1);
    dscal(N - K, TAU[I], Y(K + 1, I).asArray(), 1);

    // Compute T(1:I,I)

    dscal(I - 1, -TAU[I], T(1, I).asArray(), 1);
    dtrmv('Upper', 'No Transpose', 'NON-UNIT', I - 1, T, LDT, T(1, I).asArray(),
        1);
    T[I][I] = TAU[I];
  }
  A[K + NB][NB] = EI;

  // Compute Y(1:K,1:NB)

  dlacpy('ALL', K, NB, A(1, 2), LDA, Y, LDY);
  dtrmm('RIGHT', 'Lower', 'NO TRANSPOSE', 'UNIT', K, NB, ONE, A(K + 1, 1), LDA,
      Y, LDY);
  if (N > K + NB) {
    dgemm('NO TRANSPOSE', 'NO TRANSPOSE', K, NB, N - K - NB, ONE, A(1, 2 + NB),
        LDA, A(K + 1 + NB, 1), LDA, ONE, Y, LDY);
  }
  dtrmm(
      'RIGHT', 'Upper', 'NO TRANSPOSE', 'NON-UNIT', K, NB, ONE, T, LDT, Y, LDY);
}
