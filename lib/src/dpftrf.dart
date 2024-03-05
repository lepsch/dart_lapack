import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dpotrf.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpftrf(
  final String TRANSR,
  final String UPLO,
  final int N,
  final Array<double> A_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having()(1, offset: zeroIndexedArrayOffset);
  const ONE = 1.0;
  bool LOWER, NISODD, NORMALTRANSR;
  int N1, N2, K = 0;

  // Test the input parameters.

  INFO.value = 0;
  NORMALTRANSR = lsame(TRANSR, 'N');
  LOWER = lsame(UPLO, 'L');
  if (!NORMALTRANSR && !lsame(TRANSR, 'T')) {
    INFO.value = -1;
  } else if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('DPFTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // If N is odd, set NISODD = true;
  // If N is even, set K = N/2 and NISODD = false;

  if ((N % 2) == 0) {
    K = N ~/ 2;
    NISODD = false;
  } else {
    NISODD = true;
  }

  // Set N1 and N2 depending on LOWER

  if (LOWER) {
    N2 = N ~/ 2;
    N1 = N - N2;
  } else {
    N1 = N ~/ 2;
    N2 = N - N1;
  }

  // start execution: there are eight cases

  if (NISODD) {
    // N is odd

    if (NORMALTRANSR) {
      // N is odd and TRANSR = 'N'

      if (LOWER) {
        // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
        // T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
        // T1 -> a(0), T2 -> a(n), S -> a(n1)

        dpotrf('L', N1, A(0).asMatrix(N), N, INFO);
        if (INFO.value > 0) return;
        dtrsm('R', 'L', 'T', 'N', N2, N1, ONE, A(0).asMatrix(N), N,
            A(N1).asMatrix(N), N);
        dsyrk('U', 'N', N2, N1, -ONE, A(N1).asMatrix(N), N, ONE,
            A(N).asMatrix(N), N);
        dpotrf('U', N2, A(N).asMatrix(N), N, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + N1;
      } else {
        // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
        // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
        // T1 -> a(n2), T2 -> a(n1), S -> a(0)

        dpotrf('L', N1, A(N2).asMatrix(N), N, INFO);
        if (INFO.value > 0) return;
        dtrsm('L', 'L', 'N', 'N', N1, N2, ONE, A(N2).asMatrix(N), N,
            A(0).asMatrix(N), N);
        dsyrk('U', 'T', N2, N1, -ONE, A(0).asMatrix(N), N, ONE,
            A(N1).asMatrix(N), N);
        dpotrf('U', N2, A(N1).asMatrix(N), N, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + N1;
      }
    } else {
      // N is odd and TRANSR = 'T'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE and N is odd
        // T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
        // T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1

        dpotrf('U', N1, A(0).asMatrix(N1), N1, INFO);
        if (INFO.value > 0) return;
        dtrsm('L', 'U', 'T', 'N', N1, N2, ONE, A(0).asMatrix(N1), N1,
            A(N1 * N1).asMatrix(N1), N1);
        dsyrk('L', 'T', N2, N1, -ONE, A(N1 * N1).asMatrix(N1), N1, ONE,
            A(1).asMatrix(N1), N1);
        dpotrf('L', N2, A(1).asMatrix(N1), N1, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + N1;
      } else {
        // SRPA for UPPER, TRANSPOSE and N is odd
        // T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
        // T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2

        dpotrf('U', N1, A(N2 * N2).asMatrix(N2), N2, INFO);
        if (INFO.value > 0) return;
        dtrsm('R', 'U', 'N', 'N', N2, N1, ONE, A(N2 * N2).asMatrix(N2), N2,
            A(0).asMatrix(N2), N2);
        dsyrk('L', 'N', N2, N1, -ONE, A(0).asMatrix(N2), N2, ONE,
            A(N1 * N2).asMatrix(N2), N2);
        dpotrf('L', N2, A(N1 * N2).asMatrix(N2), N2, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + N1;
      }
    }
  } else {
    // N is even

    if (NORMALTRANSR) {
      // N is even and TRANSR = 'N'

      if (LOWER) {
        // SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
        // T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
        // T1 -> a(1), T2 -> a(0), S -> a(k+1)

        dpotrf('L', K, A(1).asMatrix(N + 1), N + 1, INFO);
        if (INFO.value > 0) return;
        dtrsm('R', 'L', 'T', 'N', K, K, ONE, A(1).asMatrix(N + 1), N + 1,
            A(K + 1).asMatrix(N + 1), N + 1);
        dsyrk('U', 'N', K, K, -ONE, A(K + 1).asMatrix(N + 1), N + 1, ONE,
            A(0).asMatrix(N + 1), N + 1);
        dpotrf('U', K, A(0).asMatrix(N + 1), N + 1, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + K;
      } else {
        // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
        // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
        // T1 -> a(k+1), T2 -> a(k), S -> a(0)

        dpotrf('L', K, A(K + 1).asMatrix(N + 1), N + 1, INFO);
        if (INFO.value > 0) return;
        dtrsm('L', 'L', 'N', 'N', K, K, ONE, A(K + 1).asMatrix(N + 1), N + 1,
            A(0).asMatrix(N + 1), N + 1);
        dsyrk('U', 'T', K, K, -ONE, A(0).asMatrix(N + 1), N + 1, ONE,
            A(K).asMatrix(N + 1), N + 1);
        dpotrf('U', K, A(K).asMatrix(N + 1), N + 1, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + K;
      }
    } else {
      // N is even and TRANSR = 'T'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE and N is even (see paper)
        // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
        // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

        dpotrf('U', K, A(0 + K).asMatrix(K), K, INFO);
        if (INFO.value > 0) return;
        dtrsm('L', 'U', 'T', 'N', K, K, ONE, A(K).asMatrix(K), N1,
            A(K * (K + 1)).asMatrix(K), K);
        dsyrk('L', 'T', K, K, -ONE, A(K * (K + 1)).asMatrix(K), K, ONE,
            A(0).asMatrix(K), K);
        dpotrf('L', K, A(0).asMatrix(K), K, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + K;
      } else {
        // SRPA for UPPER, TRANSPOSE and N is even (see paper)
        // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
        // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

        dpotrf('U', K, A(K * (K + 1)).asMatrix(K), K, INFO);
        if (INFO.value > 0) return;
        dtrsm('R', 'U', 'N', 'N', K, K, ONE, A(K * (K + 1)).asMatrix(K), K,
            A(0).asMatrix(K), K);
        dsyrk('L', 'N', K, K, -ONE, A(0).asMatrix(K), K, ONE,
            A(K * K).asMatrix(K), K);
        dpotrf('L', K, A(K * K).asMatrix(K), K, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + K;
      }
    }
  }
}
