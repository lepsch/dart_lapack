import 'package:lapack/src/matrix.dart';

void Function(
  int uplo,
  int n,
  double alpha,
  Matrix<double> a,
  int lda,
  Array<double> x_head,
  Array<double> x_tail,
  int incx,
  double beta,
  Array<double> y,
  int incy,
  int prec,
) blas_dsymv2_x = (
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