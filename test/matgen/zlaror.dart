import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zgerc.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlaset.dart';

import 'zlarnd.dart';

void zlaror(
  final String SIDE,
  final String INIT,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> ISEED_,
  final Array<Complex> X_,
  final Box<int> INFO,
) {
  final A = A_.dim(LDA);
  final ISEED = ISEED_.dim(4);
  final X = X_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TOOSML = 1.0e-20;
  int IROW, ITYPE, IXFRM, J, JCOL, KBEG, NXFRM;
  double FACTOR, XABS, XNORM;
  Complex CSIGN, XNORMS;

  INFO.value = 0;
  if (N == 0 || M == 0) return;

  ITYPE = 0;
  if (lsame(SIDE, 'L')) {
    ITYPE = 1;
  } else if (lsame(SIDE, 'R')) {
    ITYPE = 2;
  } else if (lsame(SIDE, 'C')) {
    ITYPE = 3;
  } else if (lsame(SIDE, 'T')) {
    ITYPE = 4;
  }

  // Check for argument errors.

  if (ITYPE == 0) {
    INFO.value = -1;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0 || (ITYPE == 3 && N != M)) {
    INFO.value = -4;
  } else if (LDA < M) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZLAROR', -INFO.value);
    return;
  }

  if (ITYPE == 1) {
    NXFRM = M;
  } else {
    NXFRM = N;
  }

  // Initialize A to the identity matrix if desired

  if (lsame(INIT, 'I')) zlaset('Full', M, N, Complex.zero, Complex.one, A, LDA);

  // If no rotation possible, still multiply by
  // a random complex number from the circle |x| = 1

  // 2)      Compute Rotation by computing Householder
  //         Transformations H(2), H(3), ..., H(n).  Note that the
  //         order in which they are computed is irrelevant.

  for (J = 1; J <= NXFRM; J++) {
    // 10
    X[J] = Complex.zero;
  } // 10

  for (IXFRM = 2; IXFRM <= NXFRM; IXFRM++) {
    // 30
    KBEG = NXFRM - IXFRM + 1;

    // Generate independent normal( 0, 1 ) random numbers

    for (J = KBEG; J <= NXFRM; J++) {
      // 20
      X[J] = zlarnd(3, ISEED);
    } // 20

    // Generate a Householder transformation from the random vector X

    XNORM = dznrm2(IXFRM, X(KBEG), 1);
    XABS = (X[KBEG]).abs();
    if (XABS != ZERO) {
      CSIGN = X[KBEG] / XABS.toComplex();
    } else {
      CSIGN = Complex.one;
    }
    XNORMS = CSIGN * XNORM.toComplex();
    X[NXFRM + KBEG] = -CSIGN;
    FACTOR = XNORM * (XNORM + XABS);
    if ((FACTOR).abs() < TOOSML) {
      INFO.value = 1;
      xerbla('ZLAROR', -INFO.value);
      return;
    } else {
      FACTOR = ONE / FACTOR;
    }
    X[KBEG] = X[KBEG] + XNORMS;

    // Apply Householder transformation to A

    if (ITYPE == 1 || ITYPE == 3 || ITYPE == 4) {
      // Apply H(k) on the left of A

      zgemv('C', IXFRM, N, Complex.one, A(KBEG, 1), LDA, X(KBEG), 1,
          Complex.zero, X(2 * NXFRM + 1), 1);
      zgerc(IXFRM, N, -FACTOR.toComplex(), X(KBEG), 1, X(2 * NXFRM + 1), 1,
          A(KBEG, 1), LDA);
    }

    if (ITYPE >= 2 && ITYPE <= 4) {
      // Apply H(k)* (or H(k)') on the right of A

      if (ITYPE == 4) {
        zlacgv(IXFRM, X(KBEG), 1);
      }

      zgemv('N', M, IXFRM, Complex.one, A(1, KBEG), LDA, X(KBEG), 1,
          Complex.zero, X(2 * NXFRM + 1), 1);
      zgerc(M, IXFRM, -FACTOR.toComplex(), X(2 * NXFRM + 1), 1, X(KBEG), 1,
          A(1, KBEG), LDA);
    }
  } // 30

  X[1] = zlarnd(3, ISEED);
  XABS = (X[1]).abs();
  if (XABS != ZERO) {
    CSIGN = X[1] / XABS.toComplex();
  } else {
    CSIGN = Complex.one;
  }
  X[2 * NXFRM] = CSIGN;

  // Scale the matrix A by D.

  if (ITYPE == 1 || ITYPE == 3 || ITYPE == 4) {
    for (IROW = 1; IROW <= M; IROW++) {
      // 40
      zscal(N, X[NXFRM + IROW].conjugate(), A(IROW, 1).asArray(), LDA);
    } // 40
  }

  if (ITYPE == 2 || ITYPE == 3) {
    for (JCOL = 1; JCOL <= N; JCOL++) {
      // 50
      zscal(M, X[NXFRM + JCOL], A(1, JCOL).asArray(), 1);
    } // 50
  }

  if (ITYPE == 4) {
    for (JCOL = 1; JCOL <= N; JCOL++) {
      // 60
      zscal(M, X[NXFRM + JCOL].conjugate(), A(1, JCOL).asArray(), 1);
    } // 60
  }
}
