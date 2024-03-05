import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zunmqr.dart';

void zunmhr(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final C = C_.having(ld: LDC);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  bool LEFT, LQUERY;
  int I1, I2, LWKOPT = 0, MI, NB, NH, NI, NQ, NW;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  NH = IHI - ILO;
  LEFT = lsame(SIDE, 'L');
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
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (ILO < 1 || ILO > max(1, NQ)) {
    INFO.value = -5;
  } else if (IHI < min(ILO, NQ) || IHI > NQ) {
    INFO.value = -6;
  } else if (LDA < max(1, NQ)) {
    INFO.value = -8;
  } else if (LDC < max(1, M)) {
    INFO.value = -11;
  } else if (LWORK < NW && !LQUERY) {
    INFO.value = -13;
  }

  if (INFO.value == 0) {
    if (LEFT) {
      NB = ilaenv(1, 'ZUNMQR', SIDE + TRANS, NH, N, NH, -1);
    } else {
      NB = ilaenv(1, 'ZUNMQR', SIDE + TRANS, M, NH, NH, -1);
    }
    LWKOPT = NW * NB;
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZUNMHR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0 || NH == 0) {
    WORK[1] = Complex.one;
    return;
  }

  if (LEFT) {
    MI = NH;
    NI = N;
    I1 = ILO + 1;
    I2 = 1;
  } else {
    MI = M;
    NI = NH;
    I1 = 1;
    I2 = ILO + 1;
  }

  zunmqr(SIDE, TRANS, MI, NI, NH, A(ILO + 1, ILO), LDA, TAU(ILO), C(I1, I2),
      LDC, WORK, LWORK, IINFO);

  WORK[1] = LWKOPT.toComplex();
}
