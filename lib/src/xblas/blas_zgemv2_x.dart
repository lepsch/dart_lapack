import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void Function(
  int trans,
  int m,
  int n,
  Complex alpha,
  Matrix<Complex> A,
  int lda,
  Array<Complex> head_x,
  Array<Complex> tail_x,
  int incx,
  Complex beta,
  Array<Complex> y,
  int incy,
  int prec,
) blas_zgemv2_x = (
  trans,
  m,
  n,
  alpha,
  A,
  lda,
  head_x,
  tail_x,
  incx,
  beta,
  y,
  incy,
  prec,
) {};