import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dormql.dart';
import 'package:lapack/src/dormqr.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dormtr(
  final String SIDE,
  final String UPLO,
  final String TRANS,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  bool LEFT, LQUERY, UPPER;
  int I1, I2, LWKOPT = 0, MI, NB, NI, NQ, NW;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);

  // NQ is the order of Q and NW is the minimum dimension of WORK

  if (LEFT) {
    NQ = M;
    NW = max(1, N);
  } else {
    NQ = N;
    NW = max(1, M);
  }
  if (!LEFT && !lsame(SIDE, 'R')) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T')) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, NQ)) {
    INFO.value = -7;
  } else if (LDC < max(1, M)) {
    INFO.value = -10;
  } else if (LWORK < NW && !LQUERY) {
    INFO.value = -12;
  }

  if (INFO.value == 0) {
    if (UPPER) {
      if (LEFT) {
        NB = ilaenv(1, 'DORMQL', SIDE + TRANS, M - 1, N, M - 1, -1);
      } else {
        NB = ilaenv(1, 'DORMQL', SIDE + TRANS, M, N - 1, N - 1, -1);
      }
    } else {
      if (LEFT) {
        NB = ilaenv(1, 'DORMQR', SIDE + TRANS, M - 1, N, M - 1, -1);
      } else {
        NB = ilaenv(1, 'DORMQR', SIDE + TRANS, M, N - 1, N - 1, -1);
      }
    }
    LWKOPT = NW * NB;
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DORMTR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0 || NQ == 1) {
    WORK[1] = 1;
    return;
  }

  if (LEFT) {
    MI = M - 1;
    NI = N;
  } else {
    MI = M;
    NI = N - 1;
  }

  if (UPPER) {
    // Q was determined by a call to DSYTRD with UPLO = 'U'

    dormql(SIDE, TRANS, MI, NI, NQ - 1, A(1, 2), LDA, TAU, C, LDC, WORK, LWORK,
        IINFO);
  } else {
    // Q was determined by a call to DSYTRD with UPLO = 'L'

    if (LEFT) {
      I1 = 2;
      I2 = 1;
    } else {
      I1 = 1;
      I2 = 2;
    }
    dormqr(SIDE, TRANS, MI, NI, NQ - 1, A(2, 1), LDA, TAU, C(I1, I2), LDC, WORK,
        LWORK, IINFO);
  }
  WORK[1] = LWKOPT.toDouble();
}
