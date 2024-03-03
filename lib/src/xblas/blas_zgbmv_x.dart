import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void Function(
  int trans,
  int m,
  int n,
  int kl,
  int ku,
  Complex alpha,
  Matrix<Complex> A,
  int lda,
  Array<Complex> x,
  int incx,
  Complex beta,
  Array<Complex> y,
  int incy,
  int prec,
) blas_zgbmv_x = (
  trans,
  m,
  n,
  kl,
  ku,
  alpha,
  A,
  lda,
  x,
  incx,
  beta,
  y,
  incy,
  prec,
) {};
