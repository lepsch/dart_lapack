// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dlartg.dart';
import 'package:dart_lapack/src/intrinsics/huge.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlags2.dart';
import 'package:dart_lapack/src/zlapll.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zrot.dart';

void ztgsja(
  final String JOBU,
  final String JOBV,
  final String JOBQ,
  final int M,
  final int P,
  final int N,
  final int K,
  final int L,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final double TOLA,
  final double TOLB,
  final Array<double> ALPHA_,
  final Array<double> BETA_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<Complex> WORK_,
  final Box<int> NCYCLE,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final Q = Q_.having(ld: LDQ);
  final WORK = WORK_.having();
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  const MAXIT = 40;
  const ZERO = 0.0, ONE = 1.0;
  bool INITQ, INITU, INITV, UPPER, WANTQ, WANTU, WANTV;
  int I, J, KCYCLE;
  double A1, A3, B1, B3, ERROR, GAMMA;
  Complex A2, B2;
  final CSU = Box(0.0),
      CSV = Box(0.0),
      CSQ = Box(0.0),
      RWK = Box(0.0),
      SSMIN = Box(0.0);
  final SNU = Box(Complex.zero),
      SNV = Box(Complex.zero),
      SNQ = Box(Complex.zero);
  final HUGENUM = huge(ZERO);

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
    xerbla('ZTGSJA', -INFO.value);
    return;
  }

  // Initialize U, V and Q, if necessary

  if (INITU) zlaset('Full', M, M, Complex.zero, Complex.one, U, LDU);
  if (INITV) zlaset('Full', P, P, Complex.zero, Complex.one, V, LDV);
  if (INITQ) zlaset('Full', N, N, Complex.zero, Complex.one, Q, LDQ);

  // Loop until convergence
  var converged = false;
  UPPER = false;
  for (KCYCLE = 1; KCYCLE <= MAXIT; KCYCLE++) {
    UPPER = !UPPER;

    for (I = 1; I <= L - 1; I++) {
      for (J = I + 1; J <= L; J++) {
        A1 = ZERO;
        A2 = Complex.zero;
        A3 = ZERO;
        if (K + I <= M) A1 = A[K + I][N - L + I].real;
        if (K + J <= M) A3 = A[K + J][N - L + J].real;

        B1 = B[I][N - L + I].real;
        B3 = B[J][N - L + J].real;

        if (UPPER) {
          if (K + I <= M) A2 = A[K + I][N - L + J];
          B2 = B[I][N - L + J];
        } else {
          if (K + J <= M) A2 = A[K + J][N - L + I];
          B2 = B[J][N - L + I];
        }

        zlags2(UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ);

        // Update (K+I)-th and (K+J)-th rows of matrix A: U**H *A

        if (K + J <= M) {
          zrot(
              L,
              A(K + J, N - L + 1).asArray(),
              LDA,
              A(K + I, N - L + 1).asArray(),
              LDA,
              CSU.value,
              SNU.value.conjugate());
        }

        // Update I-th and J-th rows of matrix B: V**H *B

        zrot(L, B(J, N - L + 1).asArray(), LDB, B(I, N - L + 1).asArray(), LDB,
            CSV.value, SNV.value.conjugate());

        // Update (N-L+I)-th and (N-L+J)-th columns of matrices
        // A and B: A*Q and B*Q

        zrot(min(K + L, M), A(1, N - L + J).asArray(), 1,
            A(1, N - L + I).asArray(), 1, CSQ.value, SNQ.value);

        zrot(L, B(1, N - L + J).asArray(), 1, B(1, N - L + I).asArray(), 1,
            CSQ.value, SNQ.value);

        if (UPPER) {
          if (K + I <= M) A[K + I][N - L + J] = Complex.zero;
          B[I][N - L + J] = Complex.zero;
        } else {
          if (K + J <= M) A[K + J][N - L + I] = Complex.zero;
          B[J][N - L + I] = Complex.zero;
        }

        // Ensure that the diagonal elements of A and B are real.

        if (K + I <= M) {
          A[K + I][N - L + I] = A[K + I][N - L + I].real.toComplex();
        }
        if (K + J <= M) {
          A[K + J][N - L + J] = A[K + J][N - L + J].real.toComplex();
        }
        B[I][N - L + I] = B[I][N - L + I].real.toComplex();
        B[J][N - L + J] = B[J][N - L + J].real.toComplex();

        // Update unitary matrices U, V, Q, if desired.

        if (WANTU && K + J <= M) {
          zrot(M, U(1, K + J).asArray(), 1, U(1, K + I).asArray(), 1, CSU.value,
              SNU.value);
        }

        if (WANTV) {
          zrot(P, V(1, J).asArray(), 1, V(1, I).asArray(), 1, CSV.value,
              SNV.value);
        }

        if (WANTQ) {
          zrot(N, Q(1, N - L + J).asArray(), 1, Q(1, N - L + I).asArray(), 1,
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
        zcopy(L - I + 1, A(K + I, N - L + I).asArray(), LDA, WORK, 1);
        zcopy(L - I + 1, B(I, N - L + I).asArray(), LDB, WORK(L + 1), 1);
        zlapll(L - I + 1, WORK, 1, WORK(L + 1), 1, SSMIN);
        ERROR = max(ERROR, SSMIN.value);
      }

      if (ERROR.abs() <= min(TOLA, TOLB)) {
        converged = true;
        break;
      }
    }

    // End of cycle loop
  }

  if (!converged) {
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
      A1 = A[K + I][N - L + I].real;
      B1 = B[I][N - L + I].real;
      GAMMA = B1 / A1;

      if ((GAMMA <= HUGENUM) && (GAMMA >= -HUGENUM)) {
        if (GAMMA < ZERO) {
          zdscal(L - I + 1, -ONE, B(I, N - L + I).asArray(), LDB);
          if (WANTV) zdscal(P, -ONE, V(1, I).asArray(), 1);
        }

        dlartg(GAMMA.abs(), ONE, BETA(K + I), ALPHA(K + I), RWK);

        if (ALPHA[K + I] >= BETA[K + I]) {
          zdscal(L - I + 1, ONE / ALPHA[K + I], A(K + I, N - L + I).asArray(),
              LDA);
        } else {
          zdscal(L - I + 1, ONE / BETA[K + I], B(I, N - L + I).asArray(), LDB);
          zcopy(L - I + 1, B(I, N - L + I).asArray(), LDB,
              A(K + I, N - L + I).asArray(), LDA);
        }
      } else {
        ALPHA[K + I] = ZERO;
        BETA[K + I] = ONE;
        zcopy(L - I + 1, B(I, N - L + I).asArray(), LDB,
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
