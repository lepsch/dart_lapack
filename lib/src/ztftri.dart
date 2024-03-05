import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/ztrtri.dart';

void ztftri(
  final String TRANSR,
  final String UPLO,
  final String DIAG,
  final int N,
  final Array<Complex> A_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(offset: zeroIndexedArrayOffset);
  bool LOWER, NISODD, NORMALTRANSR;
  int N1, N2, K = 0;

  // Test the input parameters.

  INFO.value = 0;
  NORMALTRANSR = lsame(TRANSR, 'N');
  LOWER = lsame(UPLO, 'L');
  if (!NORMALTRANSR && !lsame(TRANSR, 'C')) {
    INFO.value = -1;
  } else if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -2;
  } else if (!lsame(DIAG, 'N') && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZTFTRI', -INFO.value);
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

        ztrtri('L', DIAG, N1, A(0).asMatrix(), N, INFO);
        if (INFO.value > 0) return;
        ztrmm('R', 'L', 'N', DIAG, N2, N1, -Complex.one, A(0).asMatrix(), N,
            A(N1).asMatrix(), N);
        ztrtri('U', DIAG, N2, A(N).asMatrix(), N, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + N1;
        if (INFO.value > 0) return;
        ztrmm('L', 'U', 'C', DIAG, N2, N1, Complex.one, A(N).asMatrix(), N,
            A(N1).asMatrix(), N);
      } else {
        // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
        // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
        // T1 -> a(n2), T2 -> a(n1), S -> a(0)

        ztrtri('L', DIAG, N1, A(N2).asMatrix(), N, INFO);
        if (INFO.value > 0) return;
        ztrmm('L', 'L', 'C', DIAG, N1, N2, -Complex.one, A(N2).asMatrix(), N,
            A(0).asMatrix(), N);
        ztrtri('U', DIAG, N2, A(N1).asMatrix(), N, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + N1;
        if (INFO.value > 0) return;
        ztrmm('R', 'U', 'N', DIAG, N1, N2, Complex.one, A(N1).asMatrix(), N,
            A(0).asMatrix(), N);
      }
    } else {
      // N is odd and TRANSR = 'C'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE and N is odd
        // T1 -> a(0), T2 -> a(1), S -> a(0+n1*n1)

        ztrtri('U', DIAG, N1, A(0).asMatrix(), N1, INFO);
        if (INFO.value > 0) return;
        ztrmm('L', 'U', 'N', DIAG, N1, N2, -Complex.one, A(0).asMatrix(), N1,
            A(N1 * N1).asMatrix(), N1);
        ztrtri('L', DIAG, N2, A(1).asMatrix(), N1, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + N1;
        if (INFO.value > 0) return;
        ztrmm('R', 'L', 'C', DIAG, N1, N2, Complex.one, A(1).asMatrix(), N1,
            A(N1 * N1).asMatrix(), N1);
      } else {
        // SRPA for UPPER, TRANSPOSE and N is odd
        // T1 -> a(0+n2*n2), T2 -> a(0+n1*n2), S -> a(0)

        ztrtri('U', DIAG, N1, A(N2 * N2).asMatrix(), N2, INFO);
        if (INFO.value > 0) return;
        ztrmm('R', 'U', 'C', DIAG, N2, N1, -Complex.one, A(N2 * N2).asMatrix(),
            N2, A(0).asMatrix(), N2);
        ztrtri('L', DIAG, N2, A(N1 * N2).asMatrix(), N2, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + N1;
        if (INFO.value > 0) return;
        ztrmm('L', 'L', 'N', DIAG, N2, N1, Complex.one, A(N1 * N2).asMatrix(),
            N2, A(0).asMatrix(), N2);
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

        ztrtri('L', DIAG, K, A(1).asMatrix(), N + 1, INFO);
        if (INFO.value > 0) return;
        ztrmm('R', 'L', 'N', DIAG, K, K, -Complex.one, A(1).asMatrix(), N + 1,
            A(K + 1).asMatrix(), N + 1);
        ztrtri('U', DIAG, K, A(0).asMatrix(), N + 1, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + K;
        if (INFO.value > 0) return;
        ztrmm('L', 'U', 'C', DIAG, K, K, Complex.one, A(0).asMatrix(), N + 1,
            A(K + 1).asMatrix(), N + 1);
      } else {
        // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
        // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
        // T1 -> a(k+1), T2 -> a(k), S -> a(0)

        ztrtri('L', DIAG, K, A(K + 1).asMatrix(), N + 1, INFO);
        if (INFO.value > 0) return;
        ztrmm('L', 'L', 'C', DIAG, K, K, -Complex.one, A(K + 1).asMatrix(),
            N + 1, A(0).asMatrix(), N + 1);
        ztrtri('U', DIAG, K, A(K).asMatrix(), N + 1, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + K;
        if (INFO.value > 0) return;
        ztrmm('R', 'U', 'N', DIAG, K, K, Complex.one, A(K).asMatrix(), N + 1,
            A(0).asMatrix(), N + 1);
      }
    } else {
      // N is even and TRANSR = 'C'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE and N is even (see paper)
        // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
        // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

        ztrtri('U', DIAG, K, A(K).asMatrix(), K, INFO);
        if (INFO.value > 0) return;
        ztrmm('L', 'U', 'N', DIAG, K, K, -Complex.one, A(K).asMatrix(), K,
            A(K * (K + 1)).asMatrix(), K);
        ztrtri('L', DIAG, K, A(0).asMatrix(), K, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + K;
        if (INFO.value > 0) return;
        ztrmm('R', 'L', 'C', DIAG, K, K, Complex.one, A(0).asMatrix(), K,
            A(K * (K + 1)).asMatrix(), K);
      } else {
        // SRPA for UPPER, TRANSPOSE and N is even (see paper)
        // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
        // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

        ztrtri('U', DIAG, K, A(K * (K + 1)).asMatrix(), K, INFO);
        if (INFO.value > 0) return;
        ztrmm('R', 'U', 'C', DIAG, K, K, -Complex.one,
            A(K * (K + 1)).asMatrix(), K, A(0).asMatrix(), K);
        ztrtri('L', DIAG, K, A(K * K).asMatrix(), K, INFO);
        if (INFO.value > 0) INFO.value = INFO.value + K;
        if (INFO.value > 0) return;
        ztrmm('L', 'L', 'N', DIAG, K, K, Complex.one, A(K * K).asMatrix(), K,
            A(0).asMatrix(), K);
      }
    }
  }
}
