import 'package:lapack/src/matrix.dart';

void Function(
  int uplo,
  int n,
  double alpha,
  Matrix<double> a,
  int lda,
  Array<double> x,
  int incx,
  double beta,
  Array<double> y,
  int incy,
  int prec,
) blas_dsymv_x = (
  uplo,
  n,
  alpha,
  a,
  lda,
  x,
  incx,
  beta,
  y,
  incy,
  prec,
) {};
