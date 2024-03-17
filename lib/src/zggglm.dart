import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zggqrf.dart';
import 'package:lapack/src/ztrtrs.dart';
import 'package:lapack/src/zunmqr.dart';
import 'package:lapack/src/zunmrq.dart';

void zggglm(
  final int N,
  final int M,
  final int P,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> D_,
  final Array<Complex> X_,
  final Array<Complex> Y_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  final D = D_.having();
  final X = X_.having();
  final Y = Y_.having();
  bool LQUERY;
  int I, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3, NB4, NP;

  // Test the input parameters

  INFO.value = 0;
  NP = min(N, P);
  LQUERY = (LWORK == -1);
  if (N < 0) {
    INFO.value = -1;
  } else if (M < 0 || M > N) {
    INFO.value = -2;
  } else if (P < 0 || P < N - M) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }

  // Calculate workspace

  if (INFO.value == 0) {
    if (N == 0) {
      LWKMIN = 1;
      LWKOPT = 1;
    } else {
      NB1 = ilaenv(1, 'ZGEQRF', ' ', N, M, -1, -1);
      NB2 = ilaenv(1, 'ZGERQF', ' ', N, M, -1, -1);
      NB3 = ilaenv(1, 'ZUNMQR', ' ', N, M, P, -1);
      NB4 = ilaenv(1, 'ZUNMRQ', ' ', N, M, P, -1);
      NB = max(max(NB1, NB2), max(NB3, NB4));
      LWKMIN = M + N + P;
      LWKOPT = M + NP + max(N, P) * NB;
    }
    WORK[1] = LWKOPT.toComplex();

    if (LWORK < LWKMIN && !LQUERY) {
      INFO.value = -12;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZGGGLM', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    for (I = 1; I <= M; I++) {
      X[I] = Complex.zero;
    }
    for (I = 1; I <= P; I++) {
      Y[I] = Complex.zero;
    }
    return;
  }

  // Compute the GQR factorization of matrices A and B:

  // Q**H*A = ( R11 ) M,    Q**H*B*Z**H = ( T11   T12 ) M
  //          (  0  ) N-M                 (  0    T22 ) N-M
  //             M                         M+P-N  N-M

  // where R11 and T22 are upper triangular, and Q and Z are
  // unitary.

  zggqrf(N, M, P, A, LDA, WORK, B, LDB, WORK(M + 1), WORK(M + NP + 1),
      LWORK - M - NP, INFO);
  LOPT = WORK[M + NP + 1].toInt();

  // Update left-hand-side vector d = Q**H*d = ( d1 ) M
  //                                           ( d2 ) N-M

  zunmqr('Left', 'Conjugate transpose', N, 1, M, A, LDA, WORK,
      D.asMatrix(max(1, N)), max(1, N), WORK(M + NP + 1), LWORK - M - NP, INFO);
  LOPT = max(LOPT, WORK[M + NP + 1].toInt());

  // Solve T22*y2 = d2 for y2

  if (N > M) {
    ztrtrs('Upper', 'No transpose', 'Non unit', N - M, 1,
        B(M + 1, M + P - N + 1), LDB, D(M + 1).asMatrix(N - M), N - M, INFO);

    if (INFO.value > 0) {
      INFO.value = 1;
      return;
    }

    zcopy(N - M, D(M + 1), 1, Y(M + P - N + 1), 1);
  }

  // Set y1 = 0

  for (I = 1; I <= M + P - N; I++) {
    // 10
    Y[I] = Complex.zero;
  } // 10

  // Update d1 -= T12*y2

  zgemv('No transpose', M, N - M, -Complex.one, B(1, M + P - N + 1), LDB,
      Y(M + P - N + 1), 1, Complex.one, D, 1);

  // Solve triangular system: R11*x = d1

  if (M > 0) {
    ztrtrs('Upper', 'No Transpose', 'Non unit', M, 1, A, LDA, D.asMatrix(M), M,
        INFO);

    if (INFO.value > 0) {
      INFO.value = 2;
      return;
    }

    // Copy D to X

    zcopy(M, D, 1, X, 1);
  }

  // Backward transformation y = Z**H *y

  zunmrq(
      'Left',
      'Conjugate transpose',
      P,
      1,
      NP,
      B(max(1, N - P + 1), 1),
      LDB,
      WORK(M + 1),
      Y.asMatrix(max(1, P)),
      max(1, P),
      WORK(M + NP + 1),
      LWORK - M - NP,
      INFO);
  WORK[1] = (M + NP + max(LOPT, WORK[M + NP + 1].toInt())).toComplex();
}
