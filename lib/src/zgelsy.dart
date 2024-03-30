import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgeqp3.dart';
import 'package:lapack/src/zlaic1.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztzrzf.dart';
import 'package:lapack/src/zunmqr.dart';
import 'package:lapack/src/zunmrz.dart';

void zgelsy(
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<int> JPVT_,
  final double RCOND,
  final Box<int> RANK,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final JPVT = JPVT_.having();
  const IMAX = 1, IMIN = 2;
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY;
  int I, IASCL, IBSCL, ISMAX, ISMIN, J, LWKOPT, MN, NB, NB1, NB2, NB3, NB4;
  double ANRM, BIGNUM, BNRM, SMAX, SMIN, SMLNUM, WSIZE;
  final C1 = Box(Complex.zero),
      C2 = Box(Complex.zero),
      S1 = Box(Complex.zero),
      S2 = Box(Complex.zero);
  final SMINPR = Box(0.0), SMAXPR = Box(0.0);

  MN = min(M, N);
  ISMIN = MN + 1;
  ISMAX = 2 * MN + 1;

  // Test the input arguments.

  INFO.value = 0;
  NB1 = ilaenv(1, 'ZGEQRF', ' ', M, N, -1, -1);
  NB2 = ilaenv(1, 'ZGERQF', ' ', M, N, -1, -1);
  NB3 = ilaenv(1, 'ZUNMQR', ' ', M, N, NRHS, -1);
  NB4 = ilaenv(1, 'ZUNMRQ', ' ', M, N, NRHS, -1);
  NB = max(max(NB1, NB2), max(NB3, NB4));
  LWKOPT = max(1, max(MN + 2 * N + NB * (N + 1), 2 * MN + NB * NRHS));
  WORK[1] = LWKOPT.toComplex();
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDB < max(1, max(M, N))) {
    INFO.value = -7;
  } else if (LWORK < (MN + max(2 * MN, max(N + 1, MN + NRHS))) && !LQUERY) {
    INFO.value = -12;
  }

  if (INFO.value != 0) {
    xerbla('ZGELSY', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (min(M, min(N, NRHS)) == 0) {
    RANK.value = 0;
    return;
  }

  // Get machine parameters

  SMLNUM = dlamch('S') / dlamch('P');
  BIGNUM = ONE / SMLNUM;

  // Scale A, B if max entries outside range [SMLNUM,BIGNUM]

  ANRM = zlange('M', M, N, A, LDA, RWORK);
  IASCL = 0;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO);
    IASCL = 1;
  } else if (ANRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO);
    IASCL = 2;
  } else if (ANRM == ZERO) {
    // Matrix all zero. Return zero solution.

    zlaset('F', max(M, N), NRHS, Complex.zero, Complex.zero, B, LDB);
    RANK.value = 0;
    WORK[1] = LWKOPT.toComplex();
    return;
  }

  BNRM = zlange('M', M, NRHS, B, LDB, RWORK);
  IBSCL = 0;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    zlascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO);
    IBSCL = 1;
  } else if (BNRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    zlascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO);
    IBSCL = 2;
  }

  // Compute QR factorization with column pivoting of A:
  //    A * P = Q * R

  zgeqp3(M, N, A, LDA, JPVT, WORK(1), WORK(MN + 1), LWORK - MN, RWORK, INFO);
  WSIZE = MN + WORK[MN + 1].real;

  // complex workspace: MN+NB*(N+1). real workspace 2*N.
  // Details of Householder rotations stored in WORK(1:MN).

  // Determine RANK using incremental condition estimation

  WORK[ISMIN] = Complex.one;
  WORK[ISMAX] = Complex.one;
  SMAX = A[1][1].abs();
  SMIN = SMAX;
  if (A[1][1].abs() == ZERO) {
    RANK.value = 0;
    zlaset('F', max(M, N), NRHS, Complex.zero, Complex.zero, B, LDB);
    WORK[1] = LWKOPT.toComplex();
    return;
  } else {
    RANK.value = 1;
  }
  while (RANK.value < MN) {
    I = RANK.value + 1;
    zlaic1(IMIN, RANK.value, WORK(ISMIN), SMIN, A(1, I).asArray(), A[I][I],
        SMINPR, S1, C1);
    zlaic1(IMAX, RANK.value, WORK(ISMAX), SMAX, A(1, I).asArray(), A[I][I],
        SMAXPR, S2, C2);

    if (SMAXPR.value * RCOND <= SMINPR.value) {
      for (I = 1; I <= RANK.value; I++) {
        WORK[ISMIN + I - 1] = S1.value * WORK[ISMIN + I - 1];
        WORK[ISMAX + I - 1] = S2.value * WORK[ISMAX + I - 1];
      }
      WORK[ISMIN + RANK.value] = C1.value;
      WORK[ISMAX + RANK.value] = C2.value;
      SMIN = SMINPR.value;
      SMAX = SMAXPR.value;
      RANK.value++;
      continue;
    }
    break;
  }

  // complex workspace: 3*MN.

  // Logically partition R = [ R11 R12 ]
  //                         [  0  R22 ]
  // where R11 = R(1:RANK,1:RANK)

  // [R11,R12] = [ T11, 0 ] * Y

  if (RANK.value < N) {
    ztzrzf(RANK.value, N, A, LDA, WORK(MN + 1), WORK(2 * MN + 1),
        LWORK - 2 * MN, INFO);
  }

  // complex workspace: 2*MN.
  // Details of Householder rotations stored in WORK(MN+1:2*MN)

  // B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)

  zunmqr('Left', 'Conjugate transpose', M, NRHS, MN, A, LDA, WORK(1), B, LDB,
      WORK(2 * MN + 1), LWORK - 2 * MN, INFO);
  WSIZE = max(WSIZE, 2 * MN + WORK[2 * MN + 1].real);

  // complex workspace: 2*MN+NB*NRHS.

  // B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)

  ztrsm('Left', 'Upper', 'No transpose', 'Non-unit', RANK.value, NRHS,
      Complex.one, A, LDA, B, LDB);

  for (J = 1; J <= NRHS; J++) {
    for (I = RANK.value + 1; I <= N; I++) {
      B[I][J] = Complex.zero;
    }
  }

  // B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS)

  if (RANK.value < N) {
    zunmrz('Left', 'Conjugate transpose', N, NRHS, RANK.value, N - RANK.value,
        A, LDA, WORK(MN + 1), B, LDB, WORK(2 * MN + 1), LWORK - 2 * MN, INFO);
  }

  // complex workspace: 2*MN+NRHS.

  // B(1:N,1:NRHS) := P * B(1:N,1:NRHS)

  for (J = 1; J <= NRHS; J++) {
    for (I = 1; I <= N; I++) {
      WORK[JPVT[I]] = B[I][J];
    }
    zcopy(N, WORK(1), 1, B(1, J).asArray(), 1);
  }

  // complex workspace: N.

  // Undo scaling

  if (IASCL == 1) {
    zlascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO);
    zlascl('U', 0, 0, SMLNUM, ANRM, RANK.value, RANK.value, A, LDA, INFO);
  } else if (IASCL == 2) {
    zlascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO);
    zlascl('U', 0, 0, BIGNUM, ANRM, RANK.value, RANK.value, A, LDA, INFO);
  }
  if (IBSCL == 1) {
    zlascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO);
  } else if (IBSCL == 2) {
    zlascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO);
  }

  WORK[1] = LWKOPT.toComplex();
}
