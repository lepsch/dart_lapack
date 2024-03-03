import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlauum.dart';
import 'package:lapack/src/ztftri.dart';

void zpftri(
  final String TRANSR,
  final String UPLO,
  final int N,
  final Array<Complex> A_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim()(1, offset: zeroIndexedArrayOffset);
  const ONE = 1.0;
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
  } else if (N < 0) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('ZPFTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Invert the triangular Cholesky factor U or L.

  ztftri(TRANSR, UPLO, 'N', N, A, INFO);
  if (INFO.value > 0) return;

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

  // Start execution of triangular matrix multiply: inv(U)*inv(U)^C or
  // inv(L)^C*inv(L). There are eight cases.

  if (NISODD) {
    // N is odd

    if (NORMALTRANSR) {
      // N is odd and TRANSR = 'N'

      if (LOWER) {
        // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:N1-1) )
        // T1 -> a(0,0), T2 -> a(0,1), S -> a(N1,0)
        // T1 -> a(0), T2 -> a(n), S -> a(N1)

        zlauum('L', N1, A(0).asMatrix(), N, INFO);
        zherk('L', 'C', N1, N2, ONE, A(N1).asMatrix(), N, ONE, A(0).asMatrix(),
            N);
        ztrmm('L', 'U', 'N', 'N', N2, N1, Complex.one, A(N).asMatrix(), N,
            A(N1).asMatrix(), N);
        zlauum('U', N2, A(N).asMatrix(), N, INFO);
      } else {
        // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:N2-1)
        // T1 -> a(N1+1,0), T2 -> a(N1,0), S -> a(0,0)
        // T1 -> a(N2), T2 -> a(N1), S -> a(0)

        zlauum('L', N1, A(N2).asMatrix(), N, INFO);
        zherk('L', 'N', N1, N2, ONE, A(0).asMatrix(), N, ONE, A(N2).asMatrix(),
            N);
        ztrmm('R', 'U', 'C', 'N', N1, N2, Complex.one, A(N1).asMatrix(), N,
            A(0).asMatrix(), N);
        zlauum('U', N2, A(N1).asMatrix(), N, INFO);
      }
    } else {
      // N is odd and TRANSR = 'C'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE, and N is odd
        // T1 -> a(0), T2 -> a(1), S -> a(0+N1*N1)

        zlauum('U', N1, A(0).asMatrix(), N1, INFO);
        zherk('U', 'N', N1, N2, ONE, A(N1 * N1).asMatrix(), N1, ONE,
            A(0).asMatrix(), N1);
        ztrmm('R', 'L', 'N', 'N', N1, N2, Complex.one, A(1).asMatrix(), N1,
            A(N1 * N1).asMatrix(), N1);
        zlauum('L', N2, A(1).asMatrix(), N1, INFO);
      } else {
        // SRPA for UPPER, TRANSPOSE, and N is odd
        // T1 -> a(0+N2*N2), T2 -> a(0+N1*N2), S -> a(0)

        zlauum('U', N1, A(N2 * N2).asMatrix(), N2, INFO);
        zherk('U', 'C', N1, N2, ONE, A(0).asMatrix(), N2, ONE,
            A(N2 * N2).asMatrix(), N2);
        ztrmm('L', 'L', 'C', 'N', N2, N1, Complex.one, A(N1 * N2).asMatrix(),
            N2, A(0).asMatrix(), N2);
        zlauum('L', N2, A(N1 * N2).asMatrix(), N2, INFO);
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

        zlauum('L', K, A(1).asMatrix(), N + 1, INFO);
        zherk('L', 'C', K, K, ONE, A(K + 1).asMatrix(), N + 1, ONE,
            A(1).asMatrix(), N + 1);
        ztrmm('L', 'U', 'N', 'N', K, K, Complex.one, A(0).asMatrix(), N + 1,
            A(K + 1).asMatrix(), N + 1);
        zlauum('U', K, A(0).asMatrix(), N + 1, INFO);
      } else {
        // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
        // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
        // T1 -> a(k+1), T2 -> a(k), S -> a(0)

        zlauum('L', K, A(K + 1).asMatrix(), N + 1, INFO);
        zherk('L', 'N', K, K, ONE, A(0).asMatrix(), N + 1, ONE,
            A(K + 1).asMatrix(), N + 1);
        ztrmm('R', 'U', 'C', 'N', K, K, Complex.one, A(K).asMatrix(), N + 1,
            A(0).asMatrix(), N + 1);
        zlauum('U', K, A(K).asMatrix(), N + 1, INFO);
      }
    } else {
      // N is even and TRANSR = 'C'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE, and N is even (see paper)
        // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1),
        // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

        zlauum('U', K, A(K).asMatrix(), K, INFO);
        zherk('U', 'N', K, K, ONE, A(K * (K + 1)).asMatrix(), K, ONE,
            A(K).asMatrix(), K);
        ztrmm('R', 'L', 'N', 'N', K, K, Complex.one, A(0).asMatrix(), K,
            A(K * (K + 1)).asMatrix(), K);
        zlauum('L', K, A(0).asMatrix(), K, INFO);
      } else {
        // SRPA for UPPER, TRANSPOSE, and N is even (see paper)
        // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0),
        // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

        zlauum('U', K, A(K * (K + 1)).asMatrix(), K, INFO);
        zherk('U', 'C', K, K, ONE, A(0).asMatrix(), K, ONE,
            A(K * (K + 1)).asMatrix(), K);
        ztrmm('L', 'L', 'C', 'N', K, K, Complex.one, A(K * K).asMatrix(), K,
            A(0).asMatrix(), K);
        zlauum('L', K, A(K * K).asMatrix(), K, INFO);
      }
    }
  }
}
