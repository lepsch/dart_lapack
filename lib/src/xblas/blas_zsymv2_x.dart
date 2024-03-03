import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void Function(
  int uplo,
  int n,
  Complex alpha,
  Matrix<Complex> a,
  int lda,
  Array<Complex> x_head,
  Array<Complex> x_tail,
  int incx,
  Complex beta,
  Array<Complex> y,
  int incy,
  int prec,
) blas_zsymv2_x = (
  uplo,
  n,
  alpha,
  a,
  lda,
  x_head,
  x_tail,
  incx,
  beta,
  y,
  incy,
  prec,
) {};
