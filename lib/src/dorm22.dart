import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dtrmm.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dorm22(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int N1,
  final int N2,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Q = Q_.having(ld: LDQ);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  const ONE = 1.0;
  bool LEFT, LQUERY, NOTRAN;
  int I, LDWORK, LEN, LWKOPT = 0, NB, NQ, NW;

  // Test the input arguments

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  NOTRAN = lsame(TRANS, 'N');
  LQUERY = (LWORK == -1);

  // NQ is the order of Q;
  // NW is the minimum dimension of WORK.

  if (LEFT) {
    NQ = M;
  } else {
    NQ = N;
  }
  NW = NQ;
  if (N1 == 0 || N2 == 0) NW = 1;
  if (!LEFT && !lsame(SIDE, 'R')) {
    INFO.value = -1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T')) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (N1 < 0 || N1 + N2 != NQ) {
    INFO.value = -5;
  } else if (N2 < 0) {
    INFO.value = -6;
  } else if (LDQ < max(1, NQ)) {
    INFO.value = -8;
  } else if (LDC < max(1, M)) {
    INFO.value = -10;
  } else if (LWORK < NW && !LQUERY) {
    INFO.value = -12;
  }

  if (INFO.value == 0) {
    LWKOPT = M * N;
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DORM22', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    WORK[1] = 1;
    return;
  }

  // Degenerate cases (N1 = 0 or N2 = 0) are handled using DTRMM.

  if (N1 == 0) {
    dtrmm(SIDE, 'Upper', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC);
    WORK[1] = ONE;
    return;
  } else if (N2 == 0) {
    dtrmm(SIDE, 'Lower', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC);
    WORK[1] = ONE;
    return;
  }

  // Compute the largest chunk size available from the workspace.

  NB = max(1, min(LWORK, LWKOPT) ~/ NQ);

  if (LEFT) {
    if (NOTRAN) {
      for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) {
        LEN = min(NB, N - I + 1);
        LDWORK = M;

        // Multiply bottom part of C by Q12.

        dlacpy(
            'All', N1, LEN, C(N2 + 1, I), LDC, WORK.asMatrix(LDWORK), LDWORK);
        dtrmm('Left', 'Lower', 'No Transpose', 'Non-Unit', N1, LEN, ONE,
            Q(1, N2 + 1), LDQ, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply top part of C by Q11.

        dgemm('No Transpose', 'No Transpose', N1, LEN, N2, ONE, Q, LDQ, C(1, I),
            LDC, ONE, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply top part of C by Q21.

        dlacpy('All', N2, LEN, C(1, I), LDC, WORK(N1 + 1).asMatrix(LDWORK),
            LDWORK);
        dtrmm('Left', 'Upper', 'No Transpose', 'Non-Unit', N2, LEN, ONE,
            Q(N1 + 1, 1), LDQ, WORK(N1 + 1).asMatrix(LDWORK), LDWORK);

        // Multiply bottom part of C by Q22.

        dgemm(
            'No Transpose',
            'No Transpose',
            N2,
            LEN,
            N1,
            ONE,
            Q(N1 + 1, N2 + 1),
            LDQ,
            C(N2 + 1, I),
            LDC,
            ONE,
            WORK(N1 + 1).asMatrix(LDWORK),
            LDWORK);

        // Copy everything back.

        dlacpy('All', M, LEN, WORK.asMatrix(LDWORK), LDWORK, C(1, I), LDC);
      }
    } else {
      for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) {
        LEN = min(NB, N - I + 1);
        LDWORK = M;

        // Multiply bottom part of C by Q21**T.

        dlacpy(
            'All', N2, LEN, C(N1 + 1, I), LDC, WORK.asMatrix(LDWORK), LDWORK);
        dtrmm('Left', 'Upper', 'Transpose', 'Non-Unit', N2, LEN, ONE,
            Q(N1 + 1, 1), LDQ, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply top part of C by Q11**T.

        dgemm('Transpose', 'No Transpose', N2, LEN, N1, ONE, Q, LDQ, C(1, I),
            LDC, ONE, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply top part of C by Q12**T.

        dlacpy('All', N1, LEN, C(1, I), LDC, WORK(N2 + 1).asMatrix(LDWORK),
            LDWORK);
        dtrmm('Left', 'Lower', 'Transpose', 'Non-Unit', N1, LEN, ONE,
            Q(1, N2 + 1), LDQ, WORK(N2 + 1).asMatrix(LDWORK), LDWORK);

        // Multiply bottom part of C by Q22**T.

        dgemm('Transpose', 'No Transpose', N1, LEN, N2, ONE, Q(N1 + 1, N2 + 1),
            LDQ, C(N1 + 1, I), LDC, ONE, WORK(N2 + 1).asMatrix(LDWORK), LDWORK);

        // Copy everything back.

        dlacpy('All', M, LEN, WORK.asMatrix(LDWORK), LDWORK, C(1, I), LDC);
      }
    }
  } else {
    if (NOTRAN) {
      for (I = 1; NB < 0 ? I >= M : I <= M; I += NB) {
        LEN = min(NB, M - I + 1);
        LDWORK = LEN;

        // Multiply right part of C by Q21.

        dlacpy(
            'All', LEN, N2, C(I, N1 + 1), LDC, WORK.asMatrix(LDWORK), LDWORK);
        dtrmm('Right', 'Upper', 'No Transpose', 'Non-Unit', LEN, N2, ONE,
            Q(N1 + 1, 1), LDQ, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply left part of C by Q11.

        dgemm('No Transpose', 'No Transpose', LEN, N2, N1, ONE, C(I, 1), LDC, Q,
            LDQ, ONE, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply left part of C by Q12.

        dlacpy('All', LEN, N1, C(I, 1), LDC,
            WORK(1 + N2 * LDWORK).asMatrix(LDWORK), LDWORK);
        dtrmm('Right', 'Lower', 'No Transpose', 'Non-Unit', LEN, N1, ONE,
            Q(1, N2 + 1), LDQ, WORK(1 + N2 * LDWORK).asMatrix(LDWORK), LDWORK);

        // Multiply right part of C by Q22.

        dgemm(
            'No Transpose',
            'No Transpose',
            LEN,
            N1,
            N2,
            ONE,
            C(I, N1 + 1),
            LDC,
            Q(N1 + 1, N2 + 1),
            LDQ,
            ONE,
            WORK(1 + N2 * LDWORK).asMatrix(LDWORK),
            LDWORK);

        // Copy everything back.

        dlacpy('All', LEN, N, WORK.asMatrix(LDWORK), LDWORK, C(I, 1), LDC);
      }
    } else {
      for (I = 1; NB < 0 ? I >= M : I <= M; I += NB) {
        LEN = min(NB, M - I + 1);
        LDWORK = LEN;

        // Multiply right part of C by Q12**T.

        dlacpy(
            'All', LEN, N1, C(I, N2 + 1), LDC, WORK.asMatrix(LDWORK), LDWORK);
        dtrmm('Right', 'Lower', 'Transpose', 'Non-Unit', LEN, N1, ONE,
            Q(1, N2 + 1), LDQ, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply left part of C by Q11**T.

        dgemm('No Transpose', 'Transpose', LEN, N1, N2, ONE, C(I, 1), LDC, Q,
            LDQ, ONE, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply left part of C by Q21**T.

        dlacpy('All', LEN, N2, C(I, 1), LDC,
            WORK(1 + N1 * LDWORK).asMatrix(LDWORK), LDWORK);
        dtrmm('Right', 'Upper', 'Transpose', 'Non-Unit', LEN, N2, ONE,
            Q(N1 + 1, 1), LDQ, WORK(1 + N1 * LDWORK).asMatrix(LDWORK), LDWORK);

        // Multiply right part of C by Q22**T.

        dgemm(
            'No Transpose',
            'Transpose',
            LEN,
            N2,
            N1,
            ONE,
            C(I, N2 + 1),
            LDC,
            Q(N1 + 1, N2 + 1),
            LDQ,
            ONE,
            WORK(1 + N1 * LDWORK).asMatrix(LDWORK),
            LDWORK);

        // Copy everything back.

        dlacpy('All', LEN, N, WORK.asMatrix(LDWORK), LDWORK, C(I, 1), LDC);
      }
    }
  }

  WORK[1] = LWKOPT.toDouble();
}
