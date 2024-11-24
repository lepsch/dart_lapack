// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dpbstf.dart';
import 'package:lapack/src/dsbgst.dart';
import 'package:lapack/src/dsbtrd.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsbgv(
  final String JOBZ,
  final String UPLO,
  final int N,
  final int KA,
  final int KB,
  final Matrix<double> AB_,
  final int LDAB,
  final Matrix<double> BB_,
  final int LDBB,
  final Array<double> W_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final BB = BB_.having(ld: LDBB);
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  bool UPPER, WANTZ;
  String VECT;
  int INDE, INDWRK;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(UPPER || lsame(UPLO, 'L'))) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (KA < 0) {
    INFO.value = -4;
  } else if (KB < 0 || KB > KA) {
    INFO.value = -5;
  } else if (LDAB < KA + 1) {
    INFO.value = -7;
  } else if (LDBB < KB + 1) {
    INFO.value = -9;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -12;
  }
  if (INFO.value != 0) {
    xerbla('DSBGV', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a split Cholesky factorization of B.

  dpbstf(UPLO, N, KB, BB, LDBB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem.

  INDE = 1;
  INDWRK = INDE + N;
  dsbgst(
      JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK(INDWRK), IINFO);

  // Reduce to tridiagonal form.

  if (WANTZ) {
    VECT = 'U';
  } else {
    VECT = 'N';
  }
  dsbtrd(
      VECT, UPLO, N, KA, AB, LDAB, W, WORK(INDE), Z, LDZ, WORK(INDWRK), IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.

  if (!WANTZ) {
    dsterf(N, W, WORK(INDE), INFO);
  } else {
    dsteqr(JOBZ, N, W, WORK(INDE), Z, LDZ, WORK(INDWRK), INFO);
  }
}
