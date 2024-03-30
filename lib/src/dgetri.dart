import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dtrtri.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgetri(
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY;
  int I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB, NBMIN, NN;

  // Test the input parameters.

  INFO.value = 0;
  NB = ilaenv(1, 'DGETRI', ' ', N, -1, -1, -1);
  LWKOPT = max(1, N * NB);
  WORK[1] = LWKOPT.toDouble();

  LQUERY = (LWORK == -1);
  if (N < 0) {
    INFO.value = -1;
  } else if (LDA < max(1, N)) {
    INFO.value = -3;
  } else if (LWORK < max(1, N) && !LQUERY) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DGETRI', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
  // and the inverse is not computed.

  dtrtri('Upper', 'Non-unit', N, A, LDA, INFO);
  if (INFO.value > 0) return;

  NBMIN = 2;
  LDWORK = N;
  if (NB > 1 && NB < N) {
    IWS = max(LDWORK * NB, 1);
    if (LWORK < IWS) {
      NB = LWORK ~/ LDWORK;
      NBMIN = max(2, ilaenv(2, 'DGETRI', ' ', N, -1, -1, -1));
    }
  } else {
    IWS = N;
  }

  // Solve the equation inv(A)*L = inv(U) for inv(A).

  if (NB < NBMIN || NB >= N) {
    // Use unblocked code.

    for (J = N; J >= 1; J--) {
      // Copy current column of L to WORK and replace with zeros.

      for (I = J + 1; I <= N; I++) {
        WORK[I] = A[I][J];
        A[I][J] = ZERO;
      }

      // Compute current column of inv(A).

      if (J < N) {
        dgemv('No transpose', N, N - J, -ONE, A(1, J + 1), LDA, WORK(J + 1), 1,
            ONE, A(1, J).asArray(), 1);
      }
    }
  } else {
    // Use blocked code.

    NN = ((N - 1) ~/ NB) * NB + 1;
    for (J = NN; J >= 1; J -= NB) {
      JB = min(NB, N - J + 1);

      // Copy current block column of L to WORK and replace with
      // zeros.

      for (JJ = J; JJ <= J + JB - 1; JJ++) {
        for (I = JJ + 1; I <= N; I++) {
          WORK[I + (JJ - J) * LDWORK] = A[I][JJ];
          A[I][JJ] = ZERO;
        }
      }

      // Compute current block column of inv(A).

      if (J + JB <= N) {
        dgemm(
            'No transpose',
            'No transpose',
            N,
            JB,
            N - J - JB + 1,
            -ONE,
            A(1, J + JB),
            LDA,
            WORK(J + JB).asMatrix(LDWORK),
            LDWORK,
            ONE,
            A(1, J),
            LDA);
      }
      dtrsm('Right', 'Lower', 'No transpose', 'Unit', N, JB, ONE,
          WORK(J).asMatrix(LDWORK), LDWORK, A(1, J), LDA);
    }
  }

  // Apply column interchanges.

  for (J = N - 1; J >= 1; J--) {
    JP = IPIV[J];
    if (JP != J) dswap(N, A(1, J).asArray(), 1, A(1, JP).asArray(), 1);
  }

  WORK[1] = IWS.toDouble();
}
