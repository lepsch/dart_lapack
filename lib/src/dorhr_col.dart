import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaorhr_col_getrfnp.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dorhr_col(
  final int M,
  final int N,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> T_,
  final int LDT,
  final Array<double> D_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final D = D_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, J, JB, JBTEMP1, JBTEMP2, JNB, NPLUSONE;
  final IINFO = Box(0);
  // Test the input parameters

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || N > M) {
    INFO.value = -2;
  } else if (NB < 1) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDT < max(1, min(NB, N))) {
    INFO.value = -7;
  }

  // Handle error in the input parameters.

  if (INFO.value != 0) {
    xerbla('DORHR_COL', -INFO.value);
    return;
  }

  // Quick return if possible

  if (min(M, N) == 0) {
    return;
  }

  // On input, the M-by-N matrix A contains the orthogonal
  // M-by-N matrix Q_in.

  // (1) Compute the unit lower-trapezoidal V (ones on the diagonal
  // are not stored) by performing the "modified" LU-decomposition.

  // Q_in - ( S ) = V * U = ( V1 ) * U,
  //        ( 0 )           ( V2 )

  // where 0 is an (M-N)-by-N zero matrix.

  // (1-1) Factor V1 and U.

  dlaorhr_col_getrfnp(N, N, A, LDA, D, IINFO);

  // (1-2) Solve for V2.

  if (M > N) {
    dtrsm('R', 'U', 'N', 'N', M - N, N, ONE, A, LDA, A(N + 1, 1), LDA);
  }

  // (2) Reconstruct the block reflector T stored in T(1:NB, 1:N)
  // as a sequence of upper-triangular blocks with NB-size column
  // blocking.

  // Loop over the column blocks of size NB of the array A(1:M,1:N)
  // and the array T(1:NB,1:N), JB is the column index of a column
  // block, JNB is the column block size at each step JB.

  NPLUSONE = N + 1;
  for (JB = 1; NB < 0 ? JB >= N : JB <= N; JB += NB) {
    // (2-0) Determine the column block size JNB.

    JNB = min(NPLUSONE - JB, NB);

    // (2-1) Copy the upper-triangular part of the current JNB-by-JNB
    // diagonal block U(JB) (of the N-by-N matrix U) stored
    // in A(JB:JB+JNB-1,JB:JB+JNB-1) into the upper-triangular part
    // of the current JNB-by-JNB block T(1:JNB,JB:JB+JNB-1)
    // column-by-column, total JNB*(JNB+1)/2 elements.

    JBTEMP1 = JB - 1;
    for (J = JB; J <= JB + JNB - 1; J++) {
      dcopy(J - JBTEMP1, A(JB, J).asArray(), 1, T(1, J).asArray(), 1);
    }

    // (2-2) Perform on the upper-triangular part of the current
    // JNB-by-JNB diagonal block U(JB) (of the N-by-N matrix U) stored
    // in T(1:JNB,JB:JB+JNB-1) the following operation in place:
    // (-1)*U(JB)*S(JB), i.e the result will be stored in the upper-
    // triangular part of T(1:JNB,JB:JB+JNB-1). This multiplication
    // of the JNB-by-JNB diagonal block U(JB) by the JNB-by-JNB
    // diagonal block S(JB) of the N-by-N sign matrix S from the
    // right means changing the sign of each J-th column of the block
    // U(JB) according to the sign of the diagonal element of the block
    // S(JB), i.e. S(J,J) that is stored in the array element D(J).

    for (J = JB; J <= JB + JNB - 1; J++) {
      if (D[J] == ONE) {
        dscal(J - JBTEMP1, -ONE, T(1, J).asArray(), 1);
      }
    }

    // (2-3) Perform the triangular solve for the current block
    // matrix X(JB):

    // X(JB) * (A(JB)**T) = B(JB), where:

    // A(JB)**T  is a JNB-by-JNB unit upper-triangular
    //           coefficient block, and A(JB)=V1(JB), which
    //           is a JNB-by-JNB unit lower-triangular block
    //           stored in A(JB:JB+JNB-1,JB:JB+JNB-1).
    //           The N-by-N matrix V1 is the upper part
    //           of the M-by-N lower-trapezoidal matrix V
    //           stored in A(1:M,1:N);

    // B(JB)     is a JNB-by-JNB  upper-triangular right-hand
    //           side block, B(JB) = (-1)*U(JB)*S(JB), and
    //           B(JB) is stored in T(1:JNB,JB:JB+JNB-1);

    // X(JB)     is a JNB-by-JNB upper-triangular solution
    //           block, X(JB) is the upper-triangular block
    //           reflector T(JB), and X(JB) is stored
    //           in T(1:JNB,JB:JB+JNB-1).

    // In other words, we perform the triangular solve for the
    // upper-triangular block T(JB):

    // T(JB) * (V1(JB)**T) = (-1)*U(JB)*S(JB).

    // Even though the blocks X(JB) and B(JB) are upper-
    // triangular, the routine DTRSM will access all JNB**2
    // elements of the square T(1:JNB,JB:JB+JNB-1). Therefore,
    // we need to set to zero the elements of the block
    // T(1:JNB,JB:JB+JNB-1) below the diagonal before the call
    // to DTRSM.

    // (2-3a) Set the elements to zero.

    JBTEMP2 = JB - 2;
    for (J = JB; J <= JB + JNB - 2; J++) {
      for (I = J - JBTEMP2; I <= NB; I++) {
        T[I][J] = ZERO;
      }
    }

    // (2-3b) Perform the triangular solve.

    dtrsm('R', 'L', 'T', 'U', JNB, JNB, ONE, A(JB, JB), LDA, T(1, JB), LDT);
  }
}
