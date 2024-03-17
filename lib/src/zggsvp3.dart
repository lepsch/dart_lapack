import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgeqp3.dart';
import 'package:lapack/src/zgeqr2.dart';
import 'package:lapack/src/zgerq2.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlapmt.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zung2r.dart';
import 'package:lapack/src/zunm2r.dart';
import 'package:lapack/src/zunmr2.dart';

void zggsvp3(
  final String JOBU,
  final String JOBV,
  final String JOBQ,
  final int M,
  final int P,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final double TOLA,
  final double TOLB,
  final Box<int> K,
  final Box<int> L,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<int> IWORK_,
  final Array<double> RWORK_,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
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
  final TAU = TAU_.having();
  final IWORK = IWORK_.having();
  final RWORK = RWORK_.having();
  final WORK = WORK_.having();
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
    zgeqp3(P, N, B, LDB, IWORK, TAU, WORK, -1, RWORK, INFO);
    LWKOPT = WORK[1].toInt();
    if (WANTV) {
      LWKOPT = max(LWKOPT, P);
    }
    LWKOPT = max(LWKOPT, min(N, P));
    LWKOPT = max(LWKOPT, M);
    if (WANTQ) {
      LWKOPT = max(LWKOPT, N);
    }
    zgeqp3(M, N, A, LDA, IWORK, TAU, WORK, -1, RWORK, INFO);
    LWKOPT = max(LWKOPT, WORK[1].toInt());
    LWKOPT = max(1, LWKOPT);
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZGGSVP3', -INFO.value);
    return;
  }
  if (LQUERY) {
    return;
  }

  // QR with column pivoting of B: B*P = V*( S11 S12 )
  //                                       (  0   0  )

  for (I = 1; I <= N; I++) {
    // 10
    IWORK[I] = 0;
  } // 10
  zgeqp3(P, N, B, LDB, IWORK, TAU, WORK, LWORK, RWORK, INFO);

  // Update A := A*P

  zlapmt(FORWRD, M, N, A, LDA, IWORK);

  // Determine the effective rank of matrix B.

  L.value = 0;
  for (I = 1; I <= min(P, N); I++) {
    // 20
    if ((B[I][I]).abs() > TOLB) L.value++;
  } // 20

  if (WANTV) {
    // Copy the details of V, and form V.

    zlaset('Full', P, P, Complex.zero, Complex.zero, V, LDV);
    if (P > 1) zlacpy('Lower', P - 1, N, B(2, 1), LDB, V(2, 1), LDV);
    zung2r(P, P, min(P, N), V, LDV, TAU, WORK, INFO);
  }

  // Clean up B

  for (J = 1; J <= L.value - 1; J++) {
    // 40
    for (I = J + 1; I <= L.value; I++) {
      // 30
      B[I][J] = Complex.zero;
    } // 30
  } // 40
  if (P > L.value) {
    zlaset('Full', P - L.value, N, Complex.zero, Complex.zero,
        B(L.value + 1, 1), LDB);
  }

  if (WANTQ) {
    // Set Q = I and Update Q := Q*P

    zlaset('Full', N, N, Complex.zero, Complex.one, Q, LDQ);
    zlapmt(FORWRD, N, N, Q, LDQ, IWORK);
  }

  if (P >= L.value && N != L.value) {
    // RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z

    zgerq2(L.value, N, B, LDB, TAU, WORK, INFO);

    // Update A := A*Z**H

    zunmr2('Right', 'Conjugate transpose', M, N, L.value, B, LDB, TAU, A, LDA,
        WORK, INFO);
    if (WANTQ) {
      // Update Q := Q*Z**H

      zunmr2('Right', 'Conjugate transpose', N, N, L.value, B, LDB, TAU, Q, LDQ,
          WORK, INFO);
    }

    // Clean up B

    zlaset('Full', L.value, N - L.value, Complex.zero, Complex.zero, B, LDB);
    for (J = N - L.value + 1; J <= N; J++) {
      // 60
      for (I = J - N + L.value + 1; I <= L.value; I++) {
        // 50
        B[I][J] = Complex.zero;
      } // 50
    } // 60
  }

  // Let              N-L.value     L.value
  //            A = ( A11    A12 ) M,

  // then the following does the complete QR decomposition of A11:

  // A11 = U*(  0  T12 )*P1**H
  //         (  0   0  )

  for (I = 1; I <= N - L.value; I++) {
    // 70
    IWORK[I] = 0;
  } // 70
  zgeqp3(M, N - L.value, A, LDA, IWORK, TAU, WORK, LWORK, RWORK, INFO);

  // Determine the effective rank of A11

  K.value = 0;
  for (I = 1; I <= min(M, N - L.value); I++) {
    // 80
    if ((A[I][I]).abs() > TOLA) K.value++;
  } // 80

  // Update A12 := U**H*A12, where A12 = A( 1:M, N-L.value+1:N )

  zunm2r('Left', 'Conjugate transpose', M, L.value, min(M, N - L.value), A, LDA,
      TAU, A(1, N - L.value + 1), LDA, WORK, INFO);

  if (WANTU) {
    // Copy the details of U, and form U

    zlaset('Full', M, M, Complex.zero, Complex.zero, U, LDU);
    if (M > 1) zlacpy('Lower', M - 1, N - L.value, A(2, 1), LDA, U(2, 1), LDU);
    zung2r(M, M, min(M, N - L.value), U, LDU, TAU, WORK, INFO);
  }

  if (WANTQ) {
    // Update Q( 1:N, 1:N-L.value )  = Q( 1:N, 1:N-L.value )*P1

    zlapmt(FORWRD, N, N - L.value, Q, LDQ, IWORK);
  }

  // Clean up A: set the strictly lower triangular part of
  // A(1:K.value, 1:K.value) = 0, and A( K.value+1:M, 1:N-L.value ) = 0.

  for (J = 1; J <= K.value - 1; J++) {
    // 100
    for (I = J + 1; I <= K.value; I++) {
      // 90
      A[I][J] = Complex.zero;
    } // 90
  } // 100
  if (M > K.value) {
    zlaset('Full', M - K.value, N - L.value, Complex.zero, Complex.zero,
        A(K.value + 1, 1), LDA);
  }

  if (N - L.value > K.value) {
    // RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1

    zgerq2(K.value, N - L.value, A, LDA, TAU, WORK, INFO);

    if (WANTQ) {
      // Update Q( 1:N,1:N-L.value ) = Q( 1:N,1:N-L.value )*Z1**H

      zunmr2('Right', 'Conjugate transpose', N, N - L.value, K.value, A, LDA,
          TAU, Q, LDQ, WORK, INFO);
    }

    // Clean up A

    zlaset('Full', K.value, N - L.value - K.value, Complex.zero, Complex.zero,
        A, LDA);
    for (J = N - L.value - K.value + 1; J <= N - L.value; J++) {
      // 120
      for (I = J - N + L.value + K.value + 1; I <= K.value; I++) {
        // 110
        A[I][J] = Complex.zero;
      } // 110
    } // 120
  }

  if (M > K.value) {
    // QR factorization of A( K.value+1:M,N-L.value+1:N )

    zgeqr2(M - K.value, L.value, A(K.value + 1, N - L.value + 1), LDA, TAU,
        WORK, INFO);

    if (WANTU) {
      // Update U(:,K.value+1:M) := U(:,K.value+1:M)*U1

      zunm2r(
          'Right',
          'No transpose',
          M,
          M - K.value,
          min(M - K.value, L.value),
          A(K.value + 1, N - L.value + 1),
          LDA,
          TAU,
          U(1, K.value + 1),
          LDU,
          WORK,
          INFO);
    }

    // Clean up

    for (J = N - L.value + 1; J <= N; J++) {
      // 140
      for (I = J - N + K.value + L.value + 1; I <= M; I++) {
        // 130
        A[I][J] = Complex.zero;
      } // 130
    } // 140
  }

  WORK[1] = LWKOPT.toComplex();
}
