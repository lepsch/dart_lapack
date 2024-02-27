import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zpttrf.dart';
import 'package:lapack/src/zpttrs.dart';

void zptsv(
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<Complex> E_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final B = B_.dim(LDB);
  final E = E_.dim();
  final D = D_.dim();

  // Test the input parameters.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (NRHS < 0) {
    INFO.value = -2;
  } else if (LDB < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZPTSV ', -INFO.value);
    return;
  }

  // Compute the L*D*L**H (or U**H*D*U) factorization of A.

  zpttrf(N, D, E, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    zpttrs('Lower', N, NRHS, D, E, B, LDB, INFO);
  }
}
