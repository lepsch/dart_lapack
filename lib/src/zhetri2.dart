import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetri.dart';
import 'package:lapack/src/zhetri2x.dart';

void zhetri2(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final IPIV = IPIV_.dim();
  final WORK = WORK_.dim();
  bool UPPER, LQUERY;
  int MINSIZE, NBMAX;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);

  // Get blocksize

  NBMAX = ilaenv(1, 'ZHETRF', UPLO, N, -1, -1, -1);
  if (N == 0) {
    MINSIZE = 1;
  } else if (NBMAX >= N) {
    MINSIZE = N;
  } else {
    MINSIZE = (N + NBMAX + 1) * (NBMAX + 3);
  }

  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LWORK < MINSIZE && !LQUERY) {
    INFO.value = -7;
  }

  if (INFO.value != 0) {
    xerbla('ZHETRI2', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = MINSIZE.toComplex();
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (NBMAX >= N) {
    zhetri(UPLO, N, A, LDA, IPIV, WORK, INFO);
  } else {
    zhetri2x(UPLO, N, A, LDA, IPIV, WORK.asMatrix(), NBMAX, INFO);
  }
}
