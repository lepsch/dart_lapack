import 'package:lapack/src/matrix.dart';

void Function(
  int trans,
  int m,
  int n,
  double alpha,
  Matrix<double> A,
  int lda,
  Array<double> x,
  int incx,
  double beta,
  Array<double> y,
  int incy,
  int prec,
) blas_dgemv_x = (
  trans,
  m,
  n,
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
