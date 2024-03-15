import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void Function(
  int uplo,
  int n,
  Complex alpha,
  Matrix<Complex> a,
  int lda,
  Array<Complex> x,
  int incx,
  Complex beta,
  Array<Complex> y,
  int incy,
  int prec,
) blas_zsymv_x = (
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