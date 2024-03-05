import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dsymv.dart';
import 'package:lapack/src/blas/dsyr2.dart';
import 'package:lapack/src/matrix.dart';

void dlarfy(
  final String UPLO,
  final int N,
  final Array<double> V_,
  final int INCV,
  final double TAU,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.having();
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0, HALF = 0.5;
  double ALPHA;

  if (TAU == ZERO) return;

  // Form  w:= C * v

  dsymv(UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1);

  ALPHA = -HALF * TAU * ddot(N, WORK, 1, V, INCV);
  daxpy(N, ALPHA, V, INCV, WORK, 1);

  // C := C - v * w' - w * v'

  dsyr2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC);
}
