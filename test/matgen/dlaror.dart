import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

import 'dlarnd.dart';

void dlaror(
  final String SIDE,
  final String INIT,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> ISEED_,
  final Array<double> X_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final ISEED = ISEED_.dim();
  final X = X_.dim();
  const ZERO = 0.0, ONE = 1.0, TOOSML = 1.0e-20;
  int IROW, ITYPE, IXFRM, J, JCOL, KBEG, NXFRM;
  double FACTOR, XNORM, XNORMS;

  INFO.value = 0;
  if (N == 0 || M == 0) return;

  ITYPE = 0;
  if (lsame(SIDE, 'L')) {
    ITYPE = 1;
  } else if (lsame(SIDE, 'R')) {
    ITYPE = 2;
  } else if (lsame(SIDE, 'C') || lsame(SIDE, 'T')) {
    ITYPE = 3;
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
    xerbla('DLAROR', -INFO.value);
    return;
  }

  if (ITYPE == 1) {
    NXFRM = M;
  } else {
    NXFRM = N;
  }

  // Initialize A to the identity matrix if desired

  if (lsame(INIT, 'I')) dlaset('Full', M, N, ZERO, ONE, A, LDA);

  // If no rotation possible, multiply by random +/-1

  // Compute rotation by computing Householder transformations
  // H(2), H(3), ..., H(nhouse)

  for (J = 1; J <= NXFRM; J++) {
    X[J] = ZERO;
  }

  for (IXFRM = 2; IXFRM <= NXFRM; IXFRM++) {
    KBEG = NXFRM - IXFRM + 1;

    // Generate independent normal( 0, 1 ) random numbers

    for (J = KBEG; J <= NXFRM; J++) {
      X[J] = dlarnd(3, ISEED);
    }

    // Generate a Householder transformation from the random vector X

    XNORM = dnrm2(IXFRM, X(KBEG), 1);
    XNORMS = sign(XNORM, X[KBEG]).toDouble();
    X[KBEG + NXFRM] = sign(ONE, -X[KBEG]).toDouble();
    FACTOR = XNORMS * (XNORMS + X[KBEG]);
    if ((FACTOR).abs() < TOOSML) {
      INFO.value = 1;
      xerbla('DLAROR', INFO.value);
      return;
    }
    FACTOR = ONE / FACTOR;
    X[KBEG] = X[KBEG] + XNORMS;

    // Apply Householder transformation to A

    if (ITYPE == 1 || ITYPE == 3) {
      // Apply H(k) from the left.

      dgemv('T', IXFRM, N, ONE, A(KBEG, 1), LDA, X(KBEG), 1, ZERO,
          X(2 * NXFRM + 1), 1);
      dger(IXFRM, N, -FACTOR, X(KBEG), 1, X(2 * NXFRM + 1), 1, A(KBEG, 1), LDA);
    }

    if (ITYPE == 2 || ITYPE == 3) {
      // Apply H(k) from the right.

      dgemv('N', M, IXFRM, ONE, A(1, KBEG), LDA, X(KBEG), 1, ZERO,
          X(2 * NXFRM + 1), 1);
      dger(M, IXFRM, -FACTOR, X(2 * NXFRM + 1), 1, X(KBEG), 1, A(1, KBEG), LDA);
    }
  }

  X[2 * NXFRM] = sign(ONE, dlarnd(3, ISEED)).toDouble();

  // Scale the matrix A by D.

  if (ITYPE == 1 || ITYPE == 3) {
    for (IROW = 1; IROW <= M; IROW++) {
      dscal(N, X[NXFRM + IROW], A(IROW, 1).asArray(), LDA);
    }
  }

  if (ITYPE == 2 || ITYPE == 3) {
    for (JCOL = 1; JCOL <= N; JCOL++) {
      dscal(M, X[NXFRM + JCOL], A(1, JCOL).asArray(), 1);
    }
  }
}
