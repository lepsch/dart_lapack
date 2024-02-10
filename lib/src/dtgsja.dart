import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlags2.dart';
import 'package:lapack/src/dlapll.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtgsja(
  final String JOBU,
  final String JOBV,
  final String JOBQ,
  final int M,
  final int P,
  final int N,
  final int K,
  final int L,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> B,
  final int LDB,
  final double TOLA,
  final double TOLB,
  final Array<double> ALPHA,
  final Array<double> BETA,
  final Matrix<double> U,
  final int LDU,
  final Matrix<double> V,
  final int LDV,
  final Matrix<double> Q,
  final int LDQ,
  final Array<double> WORK,
  final Box<int> NCYCLE,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const MAXIT = 40;
  const ZERO = 0.0, ONE = 1.0;
  bool INITQ, INITU, INITV, UPPER, WANTQ, WANTU, WANTV;
  int I, J, KCYCLE;
  double A1, A2, A3, B1, B2, B3, ERROR, GAMMA, SSMIN = 0;
  const HUGENUM = double.maxFinite;
  final CSU = Box(0.0),
      SNU = Box(0.0),
      CSV = Box(0.0),
      SNV = Box(0.0),
      CSQ = Box(0.0),
      SNQ = Box(0.0),
      RWK = Box(0.0);

  // Decode and test the input parameters

  INITU = lsame(JOBU, 'I');
  WANTU = INITU || lsame(JOBU, 'U');

  INITV = lsame(JOBV, 'I');
  WANTV = INITV || lsame(JOBV, 'V');

  INITQ = lsame(JOBQ, 'I');
  WANTQ = INITQ || lsame(JOBQ, 'Q');

  INFO.value = 0;
  if (!(INITU || WANTU || lsame(JOBU, 'N'))) {
    INFO.value = -1;
  } else if (!(INITV || WANTV || lsame(JOBV, 'N'))) {
    INFO.value = -2;
  } else if (!(INITQ || WANTQ || lsame(JOBQ, 'N'))) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if (P < 0) {
    INFO.value = -5;
  } else if (N < 0) {
    INFO.value = -6;
  } else if (LDA < max(1, M)) {
    INFO.value = -10;
  } else if (LDB < max(1, P)) {
    INFO.value = -12;
  } else if (LDU < 1 || (WANTU && LDU < M)) {
    INFO.value = -18;
  } else if (LDV < 1 || (WANTV && LDV < P)) {
    INFO.value = -20;
  } else if (LDQ < 1 || (WANTQ && LDQ < N)) {
    INFO.value = -22;
  }
  if (INFO.value != 0) {
    xerbla('DTGSJA', -INFO.value);
    return;
  }

  // Initialize U, V and Q, if necessary

  if (INITU) dlaset('Full', M, M, ZERO, ONE, U, LDU);
  if (INITV) dlaset('Full', P, P, ZERO, ONE, V, LDV);
  if (INITQ) dlaset('Full', N, N, ZERO, ONE, Q, LDQ);

  // Loop until convergence

  UPPER = false;
  var hasConverged = false;
  for (KCYCLE = 1; KCYCLE <= MAXIT; KCYCLE++) {
    UPPER = !UPPER;

    for (I = 1; I <= L - 1; I++) {
      for (J = I + 1; J <= L; J++) {
        A1 = ZERO;
        A2 = ZERO;
        A3 = ZERO;
        if (K + I <= M) A1 = A[K + I][N - L + I];
        if (K + J <= M) A3 = A[K + J][N - L + J];

        B1 = B[I][N - L + I];
        B3 = B[J][N - L + J];

        if (UPPER) {
          if (K + I <= M) A2 = A[K + I][N - L + J];
          B2 = B[I][N - L + J];
        } else {
          if (K + J <= M) A2 = A[K + J][N - L + I];
          B2 = B[J][N - L + I];
        }

        dlags2(UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ);

        // Update (K+I)-th and (K+J)-th rows of matrix A: U**T *A

        if (K + J <= M) {
          drot(L, A(K + J, N - L + 1).asArray(), LDA,
              A(K + I, N - L + 1).asArray(), LDA, CSU.value, SNU.value);
        }

        // Update I-th and J-th rows of matrix B: V**T *B

        drot(L, B(J, N - L + 1).asArray(), LDB, B(I, N - L + 1).asArray(), LDB,
            CSV.value, SNV.value);

        // Update (N-L+I)-th and (N-L+J)-th columns of matrices
        // A and B: A*Q and B*Q

        drot(min(K + L, M), A(1, N - L + J).asArray(), 1,
            A(1, N - L + I).asArray(), 1, CSQ.value, SNQ.value);

        drot(L, B(1, N - L + J).asArray(), 1, B(1, N - L + I).asArray(), 1,
            CSQ.value, SNQ.value);

        if (UPPER) {
          if (K + I <= M) A[K + I][N - L + J] = ZERO;
          B[I][N - L + J] = ZERO;
        } else {
          if (K + J <= M) A[K + J][N - L + I] = ZERO;
          B[J][N - L + I] = ZERO;
        }

        // Update orthogonal matrices U, V, Q, if desired.

        if (WANTU && K + J <= M) {
          drot(M, U(1, K + J).asArray(), 1, U(1, K + I).asArray(), 1, CSU.value,
              SNU.value);
        }

        if (WANTV) {
          drot(P, V(1, J).asArray(), 1, V(1, I).asArray(), 1, CSV.value,
              SNV.value);
        }

        if (WANTQ) {
          drot(N, Q(1, N - L + J).asArray(), 1, Q(1, N - L + I).asArray(), 1,
              CSQ.value, SNQ.value);
        }
      }
    }

    if (!UPPER) {
      // The matrices A13 and B13 were lower triangular at the start
      // of the cycle, and are now upper triangular.

      // Convergence test: test the parallelism of the corresponding
      // rows of A and B.

      ERROR = ZERO;
      for (I = 1; I <= min(L, M - K); I++) {
        dcopy(L - I + 1, A(K + I, N - L + I).asArray(), LDA, WORK, 1);
        dcopy(L - I + 1, B(I, N - L + I).asArray(), LDB, WORK(L + 1), 1);
        dlapll(L - I + 1, WORK, 1, WORK[L + 1], 1, SSMIN);
        ERROR = max(ERROR, SSMIN);
      }

      if ((ERROR).abs() <= min(TOLA, TOLB)) {
        hasConverged = true;
        break;
      }
    }
    // End of cycle loop
  }

  if (!hasConverged) {
    // The algorithm has not converged after MAXIT cycles.
    INFO.value = 1;
  } else {
    // If ERROR <= min(TOLA,TOLB), then the algorithm has converged.
    // Compute the generalized singular value pairs (ALPHA, BETA), and
    // set the triangular matrix R to array A.

    for (I = 1; I <= K; I++) {
      ALPHA[I] = ONE;
      BETA[I] = ZERO;
    }

    for (I = 1; I <= min(L, M - K); I++) {
      A1 = A[K + I][N - L + I];
      B1 = B[I][N - L + I];
      GAMMA = B1 / A1;

      if ((GAMMA <= HUGENUM) && (GAMMA >= -HUGENUM)) {
        // change sign if necessary

        if (GAMMA < ZERO) {
          dscal(L - I + 1, -ONE, B(I, N - L + I).asArray(), LDB);
          if (WANTV) dscal(P, -ONE, V(1, I).asArray(), 1);
        }

        dlartg((GAMMA).abs(), ONE, BETA.box(K + I), ALPHA.box(K + I), RWK);

        if (ALPHA[K + I] >= BETA[K + I]) {
          dscal(L - I + 1, ONE / ALPHA[K + I], A(K + I, N - L + I).asArray(),
              LDA);
        } else {
          dscal(L - I + 1, ONE / BETA[K + I], B(I, N - L + I).asArray(), LDB);
          dcopy(L - I + 1, B(I, N - L + I).asArray(), LDB,
              A(K + I, N - L + I).asArray(), LDA);
        }
      } else {
        ALPHA[K + I] = ZERO;
        BETA[K + I] = ONE;
        dcopy(L - I + 1, B(I, N - L + I).asArray(), LDB,
            A(K + I, N - L + I).asArray(), LDA);
      }
    }

    // Post-assignment

    for (I = M + 1; I <= K + L; I++) {
      ALPHA[I] = ZERO;
      BETA[I] = ONE;
    }

    if (K + L < N) {
      for (I = K + L + 1; I <= N; I++) {
        ALPHA[I] = ZERO;
        BETA[I] = ZERO;
      }
    }
  }
  NCYCLE.value = KCYCLE;
}
