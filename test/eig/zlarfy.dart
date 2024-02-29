import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zhemv.dart';
import 'package:lapack/src/blas/zher2.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlarfy(
  final String UPLO,
  final int N,
  final Array<Complex> V_,
  final int INCV,
  final Complex TAU,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final C = C_.dim(LDC);
  final V = V_.dim();
  final WORK = WORK_.dim();

  const HALF = Complex(0.5, 0.0);
  Complex ALPHA;

  if (TAU == Complex.zero) return;

  // Form  w:= C * v

  zhemv(UPLO, N, Complex.one, C, LDC, V, INCV, Complex.zero, WORK, 1);

  ALPHA = -HALF * TAU * zdotc(N, WORK, 1, V, INCV);
  zaxpy(N, ALPHA, V, INCV, WORK, 1);

  // C := C - v * w' - w * v'

  zher2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC);
}
