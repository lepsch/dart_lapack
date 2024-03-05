import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlassq.dart';
import 'package:lapack/src/zunbdb6.dart';

void zunbdb5(
  final int M1,
  final int M2,
  final int N,
  final Array<Complex> X1_,
  final int INCX1,
  final Array<Complex> X2_,
  final int INCX2,
  final Matrix<Complex> Q1_,
  final int LDQ1,
  final Matrix<Complex> Q2_,
  final int LDQ2,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Q1 = Q1_.having(ld: LDQ1);
  final Q2 = Q2_.having(ld: LDQ2);
  final WORK = WORK_.having();
  final X1 = X1_.having();
  final X2 = X2_.having();
  const REALZERO = 0.0;
  int I, J;
  double EPS, NORM;
  final CHILDINFO = Box(0);
  final SCL = Box(0.0), SSQ = Box(0.0);

  // Test input arguments

  INFO.value = 0;
  if (M1 < 0) {
    INFO.value = -1;
  } else if (M2 < 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (INCX1 < 1) {
    INFO.value = -5;
  } else if (INCX2 < 1) {
    INFO.value = -7;
  } else if (LDQ1 < max(1, M1)) {
    INFO.value = -9;
  } else if (LDQ2 < max(1, M2)) {
    INFO.value = -11;
  } else if (LWORK < N) {
    INFO.value = -13;
  }

  if (INFO.value != 0) {
    xerbla('ZUNBDB5', -INFO.value);
    return;
  }

  EPS = dlamch('Precision');

  // Project X onto the orthogonal complement of Q if X is nonzero

  SCL.value = REALZERO;
  SSQ.value = REALZERO;
  zlassq(M1, X1, INCX1, SCL, SSQ);
  zlassq(M2, X2, INCX2, SCL, SSQ);
  NORM = SCL.value * sqrt(SSQ.value);

  if (NORM > N * EPS) {
    // Scale vector to unit norm to avoid problems in the caller code.
    // Computing the reciprocal is undesirable but
    //  * xLASCL cannot be used because of the vector increments and
    //  * the round-off error has a negligible impact on
    //    orthogonalization.
    zscal(M1, Complex.one / NORM.toComplex(), X1, INCX1);
    zscal(M2, Complex.one / NORM.toComplex(), X2, INCX2);
    zunbdb6(M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK,
        CHILDINFO);

    // If the projection is nonzero, then return;

    if (dznrm2(M1, X1, INCX1) != REALZERO ||
        dznrm2(M2, X2, INCX2) != REALZERO) {
      return;
    }
  }

  // Project each standard basis vector e_1,...,e_M1 in turn, stopping
  // when a nonzero projection is found

  for (I = 1; I <= M1; I++) {
    for (J = 1; J <= M1; J++) {
      X1[J] = Complex.zero;
    }
    X1[I] = Complex.one;
    for (J = 1; J <= M2; J++) {
      X2[J] = Complex.zero;
    }
    zunbdb6(M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK,
        CHILDINFO);
    if (dznrm2(M1, X1, INCX1) != REALZERO ||
        dznrm2(M2, X2, INCX2) != REALZERO) {
      return;
    }
  }

  // Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn,
  // stopping when a nonzero projection is found

  for (I = 1; I <= M2; I++) {
    for (J = 1; J <= M1; J++) {
      X1[J] = Complex.zero;
    }
    for (J = 1; J <= M2; J++) {
      X2[J] = Complex.zero;
    }
    X2[I] = Complex.one;
    zunbdb6(M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK,
        CHILDINFO);
    if (dznrm2(M1, X1, INCX1) != REALZERO ||
        dznrm2(M2, X2, INCX2) != REALZERO) {
      return;
    }
  }
}
