// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dtpmv.dart';
import 'package:dart_lapack/src/blas/dtpsv.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dpptrf.dart';
import 'package:dart_lapack/src/dspevd.dart';
import 'package:dart_lapack/src/dspgst.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dspgvd(
  final int ITYPE,
  final String JOBZ,
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Array<double> BP_,
  final Array<double> W_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final BP = BP_.having();
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  bool LQUERY, UPPER, WANTZ;
  String TRANS;
  int J, LIWMIN = 0, LWMIN = 0, NEIG;

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1 || LIWORK == -1);

  INFO.value = 0;
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -2;
  } else if (!(UPPER || lsame(UPLO, 'L'))) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -9;
  }

  if (INFO.value == 0) {
    if (N <= 1) {
      LIWMIN = 1;
      LWMIN = 1;
    } else {
      if (WANTZ) {
        LIWMIN = 3 + 5 * N;
        LWMIN = 1 + 6 * N + 2 * pow(N, 2).toInt();
      } else {
        LIWMIN = 1;
        LWMIN = 2 * N;
      }
    }
    WORK[1] = LWMIN.toDouble();
    IWORK[1] = LIWMIN;
    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -11;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -13;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSPGVD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a Cholesky factorization of BP.

  dpptrf(UPLO, N, BP, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem and solve.

  dspgst(ITYPE, UPLO, N, AP, BP, INFO);
  dspevd(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO);
  LWMIN = max(LWMIN, WORK[1]).toInt();
  LIWMIN = max(LIWMIN, IWORK[1]);

  if (WANTZ) {
    // Backtransform eigenvectors to the original problem.

    NEIG = N;
    if (INFO.value > 0) NEIG = INFO.value - 1;
    if (ITYPE == 1 || ITYPE == 2) {
      // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
      // backtransform eigenvectors: x = inv(L)**T *y or inv(U)*y

      if (UPPER) {
        TRANS = 'N';
      } else {
        TRANS = 'T';
      }

      for (J = 1; J <= NEIG; J++) {
        dtpsv(UPLO, TRANS, 'Non-unit', N, BP, Z(1, J).asArray(), 1);
      }
    } else if (ITYPE == 3) {
      // For B*A*x=(lambda)*x;
      // backtransform eigenvectors: x = L*y or U**T *y

      if (UPPER) {
        TRANS = 'T';
      } else {
        TRANS = 'N';
      }

      for (J = 1; J <= NEIG; J++) {
        dtpmv(UPLO, TRANS, 'Non-unit', N, BP, Z(1, J).asArray(), 1);
      }
    }
  }

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
