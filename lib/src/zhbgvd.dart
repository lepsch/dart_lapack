// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhbgst.dart';
import 'package:lapack/src/zhbtrd.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zpbstf.dart';
import 'package:lapack/src/zstedc.dart';

void zhbgvd(
  final String JOBZ,
  final String UPLO,
  final int N,
  final int KA,
  final int KB,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> BB_,
  final int LDBB,
  final Array<double> W_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final BB = BB_.having(ld: LDBB);
  final Z = Z_.having(ld: LDZ);
  final W = W_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  bool LQUERY, UPPER, WANTZ;
  String VECT;
  int INDE, INDWK2, INDWRK, LIWMIN, LLRWK, LLWK2, LRWMIN, LWMIN;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1 || LRWORK == -1 || LIWORK == -1);

  INFO.value = 0;
  if (N <= 1) {
    LWMIN = 1 + N;
    LRWMIN = 1 + N;
    LIWMIN = 1;
  } else if (WANTZ) {
    LWMIN = 2 * pow(N, 2).toInt();
    LRWMIN = 1 + 5 * N + 2 * pow(N, 2).toInt();
    LIWMIN = 3 + 5 * N;
  } else {
    LWMIN = N;
    LRWMIN = N;
    LIWMIN = 1;
  }
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

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toComplex();
    RWORK[1] = LRWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -14;
    } else if (LRWORK < LRWMIN && !LQUERY) {
      INFO.value = -16;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -18;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHBGVD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a split Cholesky factorization of B.

  zpbstf(UPLO, N, KB, BB, LDBB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem.

  INDE = 1;
  INDWRK = INDE + N;
  INDWK2 = 1 + N * N;
  LLWK2 = LWORK - INDWK2 + 2;
  LLRWK = LRWORK - INDWRK + 2;
  zhbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, RWORK, IINFO);

  // Reduce Hermitian band matrix to tridiagonal form.

  if (WANTZ) {
    VECT = 'U';
  } else {
    VECT = 'N';
  }
  zhbtrd(VECT, UPLO, N, KA, AB, LDAB, W, RWORK(INDE), Z, LDZ, WORK, IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEDC.

  if (!WANTZ) {
    dsterf(N, W, RWORK(INDE), INFO);
  } else {
    zstedc('I', N, W, RWORK(INDE), WORK.asMatrix(N), N, WORK(INDWK2), LLWK2,
        RWORK(INDWRK), LLRWK, IWORK, LIWORK, INFO);
    zgemm('N', 'N', N, N, N, Complex.one, Z, LDZ, WORK.asMatrix(N), N,
        Complex.zero, WORK(INDWK2).asMatrix(N), N);
    zlacpy('A', N, N, WORK(INDWK2).asMatrix(N), N, Z, LDZ);
  }

  WORK[1] = LWMIN.toComplex();
  RWORK[1] = LRWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
