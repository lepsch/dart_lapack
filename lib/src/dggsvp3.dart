import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeqp3.dart';
import 'package:lapack/src/dgeqr2.dart';
import 'package:lapack/src/dgerq2.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlapmt.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorg2r.dart';
import 'package:lapack/src/dorm2r.dart';
import 'package:lapack/src/dormr2.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dggsvp3(
  final String JOBU,
  final String JOBV,
  final String JOBQ,
  final int M,
  final int P,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> B,
  final int LDB,
  final double TOLA,
  final double TOLB,
  final Box<int> K,
  final Box<int> L,
  final Matrix<double> U,
  final int LDU,
  final Matrix<double> V,
  final int LDV,
  final Matrix<double> Q,
  final int LDQ,
  final Array<int> IWORK,
  final Array<double> TAU,
  final Array<double> WORK,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool FORWRD, WANTQ, WANTU, WANTV, LQUERY;
  int I, J, LWKOPT;

  // Test the input parameters

  WANTU = lsame(JOBU, 'U');
  WANTV = lsame(JOBV, 'V');
  WANTQ = lsame(JOBQ, 'Q');
  FORWRD = true;
  LQUERY = (LWORK == -1);
  LWKOPT = 1;

  // Test the input arguments

  INFO.value = 0;
  if (!(WANTU || lsame(JOBU, 'N'))) {
    INFO.value = -1;
  } else if (!(WANTV || lsame(JOBV, 'N'))) {
    INFO.value = -2;
  } else if (!(WANTQ || lsame(JOBQ, 'N'))) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if (P < 0) {
    INFO.value = -5;
  } else if (N < 0) {
    INFO.value = -6;
  } else if (LDA < max(1, M)) {
    INFO.value = -8;
  } else if (LDB < max(1, P)) {
    INFO.value = -10;
  } else if (LDU < 1 || (WANTU && LDU < M)) {
    INFO.value = -16;
  } else if (LDV < 1 || (WANTV && LDV < P)) {
    INFO.value = -18;
  } else if (LDQ < 1 || (WANTQ && LDQ < N)) {
    INFO.value = -20;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -24;
  }

  // Compute workspace

  if (INFO.value == 0) {
    dgeqp3(P, N, B, LDB, IWORK, TAU, WORK, -1, INFO);
    LWKOPT = WORK[1].toInt();
    if (WANTV) {
      LWKOPT = max(LWKOPT, P);
    }
    LWKOPT = max(LWKOPT, min(N, P));
    LWKOPT = max(LWKOPT, M);
    if (WANTQ) {
      LWKOPT = max(LWKOPT, N);
    }
    dgeqp3(M, N, A, LDA, IWORK, TAU, WORK, -1, INFO);
    LWKOPT = max(LWKOPT, WORK[1].toInt());
    LWKOPT = max(1, LWKOPT);
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DGGSVP3', -INFO.value);
    return;
  }
  if (LQUERY) {
    return;
  }

  // QR with column pivoting of B: B*P = V*( S11 S12 )
  //                                       (  0   0  )

  for (I = 1; I <= N; I++) {
    IWORK[I] = 0;
  }
  dgeqp3(P, N, B, LDB, IWORK, TAU, WORK, LWORK, INFO);

  // Update A := A*P

  dlapmt(FORWRD, M, N, A, LDA, IWORK);

  // Determine the effective rank of matrix B.

  L.value = 0;
  for (I = 1; I <= min(P, N); I++) {
    if ((B[I][I]).abs() > TOLB) L.value = L.value + 1;
  }

  if (WANTV) {
    // Copy the details of V, and form V.

    dlaset('Full', P, P, ZERO, ZERO, V, LDV);
    if (P > 1) dlacpy('Lower', P - 1, N, B(2, 1), LDB, V(2, 1), LDV);
    dorg2r(P, P, min(P, N), V, LDV, TAU, WORK, INFO);
  }

  // Clean up B

  for (J = 1; J <= L.value - 1; J++) {
    for (I = J + 1; I <= L.value; I++) {
      B[I][J] = ZERO;
    }
  }
  if (P > L.value) {
    dlaset('Full', P - L.value, N, ZERO, ZERO, B(L.value + 1, 1), LDB);
  }

  if (WANTQ) {
    // Set Q = I and Update Q := Q*P

    dlaset('Full', N, N, ZERO, ONE, Q, LDQ);
    dlapmt(FORWRD, N, N, Q, LDQ, IWORK);
  }

  if (P >= L.value && N != L.value) {
    // RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z

    dgerq2(L.value, N, B, LDB, TAU, WORK, INFO);

    // Update A := A*Z**T

    dormr2(
        'Right', 'Transpose', M, N, L.value, B, LDB, TAU, A, LDA, WORK, INFO);

    if (WANTQ) {
      // Update Q := Q*Z**T

      dormr2(
          'Right', 'Transpose', N, N, L.value, B, LDB, TAU, Q, LDQ, WORK, INFO);
    }

    // Clean up B

    dlaset('Full', L.value, N - L.value, ZERO, ZERO, B, LDB);
    for (J = N - L.value + 1; J <= N; J++) {
      for (I = J - N + L.value + 1; I <= L.value; I++) {
        B[I][J] = ZERO;
      }
    }
  }

  // Let              N-L     L
  //            A = ( A11    A12 ) M,
  //
  // then the following does the complete QR decomposition of A11:
  //
  //          A11 = U*(  0  T12 )*P1**T
  //                  (  0   0  )

  for (I = 1; I <= N - L.value; I++) {
    IWORK[I] = 0;
  }
  dgeqp3(M, N - L.value, A, LDA, IWORK, TAU, WORK, LWORK, INFO);

  // Determine the effective rank of A11

  K.value = 0;
  for (I = 1; I <= min(M, N - L.value); I++) {
    if ((A[I][I]).abs() > TOLA) K.value = K.value + 1;
  }

  // Update A12 := U**T*A12, where A12 = A[ 1:M][ N-L.value+1:N ]

  dorm2r('Left', 'Transpose', M, L.value, min(M, N - L.value), A, LDA, TAU,
      A[1][N - L.value + 1], LDA, WORK, INFO);

  if (WANTU) {
    // Copy the details of U, and form U

    dlaset('Full', M, M, ZERO, ZERO, U, LDU);
    if (M > 1) dlacpy('Lower', M - 1, N - L.value, A(2, 1), LDA, U(2, 1), LDU);
    dorg2r(M, M, min(M, N - L.value), U, LDU, TAU, WORK, INFO);
  }

  if (WANTQ) {
    // Update Q[ 1:N][ 1:N-L.value ]  = Q[ 1:N][ 1:N-L.value ]*P1

    dlapmt(FORWRD, N, N - L.value, Q, LDQ, IWORK);
  }

  // Clean up A: set the strictly lower triangular part of
  // A[1:K.value][ 1:K.value] = 0, and A[ K.value+1:M][ 1:N-L.value ] = 0.

  for (J = 1; J <= K.value - 1; J++) {
    for (I = J + 1; I <= K.value; I++) {
      A[I][J] = ZERO;
    }
  }
  if (M > K.value) {
    dlaset(
        'Full', M - K.value, N - L.value, ZERO, ZERO, A(K.value + 1, 1), LDA);
  }

  if (N - L.value > K.value) {
    // RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1

    dgerq2(K.value, N - L.value, A, LDA, TAU, WORK, INFO);

    if (WANTQ) {
      // Update Q[ 1:N][1:N-L.value ] = Q[ 1:N][1:N-L.value ]*Z1**T

      dormr2('Right', 'Transpose', N, N - L.value, K.value, A, LDA, TAU, Q, LDQ,
          WORK, INFO);
    }

    // Clean up A

    dlaset('Full', K.value, N - L.value - K.value, ZERO, ZERO, A, LDA);
    for (J = N - L.value - K.value + 1; J <= N - L.value; J++) {
      for (I = J - N + L.value + K.value + 1; I <= K.value; I++) {
        A[I][J] = ZERO;
      }
    }
  }

  if (M > K.value) {
    // QR factorization of A[ K.value+1:M][N-L.value+1:N ]

    dgeqr2(M - K.value, L.value, A(K.value + 1, N - L.value + 1), LDA, TAU,
        WORK, INFO);

    if (WANTU) {
      // Update U[:][K.value+1:M] := U[:][K.value+1:M]*U1

      dorm2r(
          'Right',
          'No transpose',
          M,
          M - K.value,
          min(M - K.value, L.value),
          A[K.value + 1][N - L.value + 1],
          LDA,
          TAU,
          U[1][K.value + 1],
          LDU,
          WORK,
          INFO);
    }

    // Clean up

    for (J = N - L.value + 1; J <= N; J++) {
      for (I = J - N + K.value + L.value + 1; I <= M; I++) {
        A[I][J] = ZERO;
      }
    }
  }

  WORK[1] = LWKOPT.toDouble();
}
