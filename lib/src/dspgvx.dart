// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dtpmv.dart';
import 'package:dart_lapack/src/blas/dtpsv.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dpptrf.dart';
import 'package:dart_lapack/src/dspevx.dart';
import 'package:dart_lapack/src/dspgst.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dspgvx(
  final int ITYPE,
  final String JOBZ,
  final String RANGE,
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Array<double> BP_,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final double ABSTOL,
  final Box<int> M,
  final Array<double> W_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Array<int> IFAIL_,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final BP = BP_.having();
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final IFAIL = IFAIL_.having();
  bool ALLEIG, INDEIG, UPPER, VALEIG, WANTZ;
  String TRANS;
  int J;

  // Test the input parameters.

  UPPER = lsame(UPLO, 'U');
  WANTZ = lsame(JOBZ, 'V');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');

  INFO.value = 0;
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -2;
  } else if (!(ALLEIG || VALEIG || INDEIG)) {
    INFO.value = -3;
  } else if (!(UPPER || lsame(UPLO, 'L'))) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else {
    if (VALEIG) {
      if (N > 0 && VU <= VL) {
        INFO.value = -9;
      }
    } else if (INDEIG) {
      if (IL < 1) {
        INFO.value = -10;
      } else if (IU < min(N, IL) || IU > N) {
        INFO.value = -11;
      }
    }
  }
  if (INFO.value == 0) {
    if (LDZ < 1 || (WANTZ && LDZ < N)) {
      INFO.value = -16;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSPGVX', -INFO.value);
    return;
  }

  // Quick return if possible

  M.value = 0;
  if (N == 0) return;

  // Form a Cholesky factorization of B.

  dpptrf(UPLO, N, BP, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem and solve.

  dspgst(ITYPE, UPLO, N, AP, BP, INFO);
  dspevx(JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
      IWORK, IFAIL, INFO);

  if (WANTZ) {
    // Backtransform eigenvectors to the original problem.

    if (INFO.value > 0) M.value = INFO.value - 1;
    if (ITYPE == 1 || ITYPE == 2) {
      // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
      // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

      if (UPPER) {
        TRANS = 'N';
      } else {
        TRANS = 'T';
      }

      for (J = 1; J <= M.value; J++) {
        dtpsv(UPLO, TRANS, 'Non-unit', N, BP, Z(1, J).asArray(), 1);
      }
    } else if (ITYPE == 3) {
      // For B*A*x=(lambda)*x;
      // backtransform eigenvectors: x = L*y or U**T*y

      if (UPPER) {
        TRANS = 'T';
      } else {
        TRANS = 'N';
      }

      for (J = 1; J <= M.value; J++) {
        dtpmv(UPLO, TRANS, 'Non-unit', N, BP, Z(1, J).asArray(), 1);
      }
    }
  }
}
