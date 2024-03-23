import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlauu2.dart';

void zlauum(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  const ONE = 1.0;
  bool UPPER;
  int I, IB, NB;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZLAUUM', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Determine the block size for this environment.

  NB = ilaenv(1, 'ZLAUUM', UPLO, N, -1, -1, -1);

  if (NB <= 1 || NB >= N) {
    // Use unblocked code

    zlauu2(UPLO, N, A, LDA, INFO);
  } else {
    // Use blocked code

    if (UPPER) {
      // Compute the product U * U**H.

      for (I = 1; I <= N; I += NB) {
        IB = min(NB, N - I + 1);
        ztrmm('Right', 'Upper', 'Conjugate transpose', 'Non-unit', I - 1, IB,
            Complex.one, A(I, I), LDA, A(1, I), LDA);
        zlauu2('Upper', IB, A(I, I), LDA, INFO);
        if (I + IB <= N) {
          zgemm(
              'No transpose',
              'Conjugate transpose',
              I - 1,
              IB,
              N - I - IB + 1,
              Complex.one,
              A(1, I + IB),
              LDA,
              A(I, I + IB),
              LDA,
              Complex.one,
              A(1, I),
              LDA);
          zherk('Upper', 'No transpose', IB, N - I - IB + 1, ONE, A(I, I + IB),
              LDA, ONE, A(I, I), LDA);
        }
      }
    } else {
      // Compute the product L**H * L.

      for (I = 1; I <= N; I += NB) {
        IB = min(NB, N - I + 1);
        ztrmm('Left', 'Lower', 'Conjugate transpose', 'Non-unit', IB, I - 1,
            Complex.one, A(I, I), LDA, A(I, 1), LDA);
        zlauu2('Lower', IB, A(I, I), LDA, INFO);
        if (I + IB <= N) {
          zgemm(
              'Conjugate transpose',
              'No transpose',
              IB,
              I - 1,
              N - I - IB + 1,
              Complex.one,
              A(I + IB, I),
              LDA,
              A(I + IB, 1),
              LDA,
              Complex.one,
              A(I, 1),
              LDA);
          zherk('Lower', 'Conjugate transpose', IB, N - I - IB + 1, ONE,
              A(I + IB, I), LDA, ONE, A(I, I), LDA);
        }
      }
    }
  }
}
