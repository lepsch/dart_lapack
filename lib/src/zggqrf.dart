import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/qr/ll/zgeqrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgerqf.dart';
import 'package:lapack/src/zunmqr.dart';

void zggqrf(
  final int N,
  final int M,
  final int P,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAUA_,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> TAUB_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final WORK = WORK_.dim();
  final TAUA = TAUA_.dim();
  final TAUB = TAUB_.dim();
  bool LQUERY;
  int LOPT, LWKOPT, NB, NB1, NB2, NB3;

  // Test the input parameters

  INFO.value = 0;
  NB1 = ilaenv(1, 'ZGEQRF', ' ', N, M, -1, -1);
  NB2 = ilaenv(1, 'ZGERQF', ' ', N, P, -1, -1);
  NB3 = ilaenv(1, 'ZUNMQR', ' ', N, M, P, -1);
  NB = max(NB1, max(NB2, NB3));
  LWKOPT = max(1, max(N, max(M, P)) * NB);
  WORK[1] = LWKOPT.toComplex();
  LQUERY = (LWORK == -1);
  if (N < 0) {
    INFO.value = -1;
  } else if (M < 0) {
    INFO.value = -2;
  } else if (P < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  } else if (LWORK < max(max(1, N), max(M, P)) && !LQUERY) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('ZGGQRF', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // QR factorization of N-by-M matrix A: A = Q*R

  zgeqrf(N, M, A, LDA, TAUA, WORK, LWORK, INFO);
  LOPT = WORK[1].toInt();

  // Update B := Q**H*B.

  zunmqr('Left', 'Conjugate Transpose', N, P, min(N, M), A, LDA, TAUA, B, LDB,
      WORK, LWORK, INFO);
  LOPT = max(LOPT, WORK[1].toInt());

  // RQ factorization of N-by-P matrix B: B = T*Z.

  zgerqf(N, P, B, LDB, TAUB, WORK, LWORK, INFO);
  WORK[1] = max(LOPT, WORK[1].toInt()).toComplex();
}
