import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dbdsqr.dart';
import 'package:lapack/src/dgebd2.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

double dqrt12(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> S_,
  final Array<double> WORK_,
  final int LWORK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final S = S_.having();
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, ONE = 1.0;
  final DUMMY = Array<double>(1);
  final INFO = Box(0);

  // Test that enough workspace is supplied

  if (LWORK <
      max(M * N + 4 * min(M, N) + max(M, N), M * N + 2 * min(M, N) + 4 * N)) {
    xerbla('DQRT12', 7);
    return ZERO;
  }

  // Quick return if possible

  final MN = min(M, N);
  if (MN <= ZERO) return ZERO;

  final NRMSVL = dnrm2(MN, S, 1);

  // Copy upper triangle of A into work

  dlaset('Full', M, N, ZERO, ZERO, WORK.asMatrix(), M);
  for (var J = 1; J <= N; J++) {
    for (var I = 1; I <= min(J, M); I++) {
      WORK[(J - 1) * M + I] = A[I][J];
    }
  }

  // Get machine parameters

  final SMLNUM = dlamch('S') / dlamch('P');
  final BIGNUM = ONE / SMLNUM;

  // Scale work if max entry outside range [SMLNUM,BIGNUM]

  final ANRM = dlange('M', M, N, WORK.asMatrix(), M, DUMMY);
  final int ISCL;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    dlascl('G', 0, 0, ANRM, SMLNUM, M, N, WORK.asMatrix(), M, INFO);
    ISCL = 1;
  } else if (ANRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    dlascl('G', 0, 0, ANRM, BIGNUM, M, N, WORK.asMatrix(), M, INFO);
    ISCL = 1;
  } else {
    ISCL = 0;
  }

  if (ANRM != ZERO) {
    // Compute SVD of work

    dgebd2(
        M,
        N,
        WORK.asMatrix(),
        M,
        WORK(M * N + 1),
        WORK(M * N + MN + 1),
        WORK(M * N + 2 * MN + 1),
        WORK(M * N + 3 * MN + 1),
        WORK(M * N + 4 * MN + 1),
        INFO);
    dbdsqr(
        'Upper',
        MN,
        0,
        0,
        0,
        WORK(M * N + 1),
        WORK(M * N + MN + 1),
        DUMMY.asMatrix(),
        MN,
        DUMMY.asMatrix(),
        1,
        DUMMY.asMatrix(),
        MN,
        WORK(M * N + 2 * MN + 1),
        INFO);

    if (ISCL == 1) {
      if (ANRM > BIGNUM) {
        dlascl('G', 0, 0, BIGNUM, ANRM, MN, 1, WORK(M * N + 1).asMatrix(), MN,
            INFO);
      }
      if (ANRM < SMLNUM) {
        dlascl('G', 0, 0, SMLNUM, ANRM, MN, 1, WORK(M * N + 1).asMatrix(), MN,
            INFO);
      }
    }
  } else {
    for (var I = 1; I <= MN; I++) {
      WORK[M * N + I] = ZERO;
    }
  }

  // Compare s and singular values of work

  daxpy(MN, -ONE, S, 1, WORK(M * N + 1), 1);

  final result = dasum(MN, WORK(M * N + 1), 1) /
      (dlamch('Epsilon') * (max(M, N)).toDouble());

  return NRMSVL != ZERO ? result / NRMSVL : result;
}
