// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dstebz.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zhbgst.dart';
import 'package:dart_lapack/src/zhbtrd.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zpbstf.dart';
import 'package:dart_lapack/src/zstein.dart';
import 'package:dart_lapack/src/zsteqr.dart';

void zhbgvx(
  final String JOBZ,
  final String RANGE,
  final String UPLO,
  final int N,
  final int KA,
  final int KB,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> BB_,
  final int LDBB,
  final Matrix<Complex> Q_,
  final int LDQ,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final double ABSTOL,
  final Box<int> M,
  final Array<double> W_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Array<int> IFAIL_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final BB = BB_.having(ld: LDBB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final W = W_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final IFAIL = IFAIL_.having();
  const ZERO = 0.0;
  bool ALLEIG, INDEIG, TEST, UPPER, VALEIG, WANTZ;
  String ORDER, VECT;
  int I, INDD, INDE, INDEE, INDISP, INDIWK, INDRWK, INDWRK, ITMP1, J, JJ;
  double TMP1;
  final IINFO = Box(0), NSPLIT = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(ALLEIG || VALEIG || INDEIG)) {
    INFO.value = -2;
  } else if (!(UPPER || lsame(UPLO, 'L'))) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (KA < 0) {
    INFO.value = -5;
  } else if (KB < 0 || KB > KA) {
    INFO.value = -6;
  } else if (LDAB < KA + 1) {
    INFO.value = -8;
  } else if (LDBB < KB + 1) {
    INFO.value = -10;
  } else if (LDQ < 1 || (WANTZ && LDQ < N)) {
    INFO.value = -12;
  } else {
    if (VALEIG) {
      if (N > 0 && VU <= VL) INFO.value = -14;
    } else if (INDEIG) {
      if (IL < 1 || IL > max(1, N)) {
        INFO.value = -15;
      } else if (IU < min(N, IL) || IU > N) {
        INFO.value = -16;
      }
    }
  }
  if (INFO.value == 0) {
    if (LDZ < 1 || (WANTZ && LDZ < N)) {
      INFO.value = -21;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHBGVX', -INFO.value);
    return;
  }

  // Quick return if possible

  M.value = 0;
  if (N == 0) return;

  // Form a split Cholesky factorization of B.

  zpbstf(UPLO, N, KB, BB, LDBB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem.

  zhbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, WORK, RWORK, IINFO);

  // Solve the standard eigenvalue problem.
  // Reduce Hermitian band matrix to tridiagonal form.

  INDD = 1;
  INDE = INDD + N;
  INDRWK = INDE + N;
  INDWRK = 1;
  if (WANTZ) {
    VECT = 'U';
  } else {
    VECT = 'N';
  }
  zhbtrd(VECT, UPLO, N, KA, AB, LDAB, RWORK(INDD), RWORK(INDE), Q, LDQ,
      WORK(INDWRK), IINFO);

  // If all eigenvalues are desired and ABSTOL is less than or equal
  // to zero, then call DSTERF or ZSTEQR.  If this fails for some
  // eigenvalue, then try DSTEBZ.

  TEST = false;
  if (INDEIG) {
    if (IL == 1 && IU == N) {
      TEST = true;
    }
  }
  var success = false;
  if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
    dcopy(N, RWORK(INDD), 1, W, 1);
    INDEE = INDRWK + 2 * N;
    dcopy(N - 1, RWORK(INDE), 1, RWORK(INDEE), 1);
    if (!WANTZ) {
      dsterf(N, W, RWORK(INDEE), INFO);
    } else {
      zlacpy('A', N, N, Q, LDQ, Z, LDZ);
      zsteqr(JOBZ, N, W, RWORK(INDEE), Z, LDZ, RWORK(INDRWK), INFO);
      if (INFO.value == 0) {
        for (I = 1; I <= N; I++) {
          IFAIL[I] = 0;
        }
      }
    }
    if (INFO.value == 0) {
      M.value = N;
      success = true;
    }
    INFO.value = 0;
  }

  if (!success) {
    // Otherwise, call DSTEBZ and, if eigenvectors are desired,
    // call ZSTEIN.

    if (WANTZ) {
      ORDER = 'B';
    } else {
      ORDER = 'E';
    }
    INDISP = 1 + N;
    INDIWK = INDISP + N;
    dstebz(RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, RWORK(INDD), RWORK(INDE), M,
        NSPLIT, W, IWORK(1), IWORK(INDISP), RWORK(INDRWK), IWORK(INDIWK), INFO);

    if (WANTZ) {
      zstein(N, RWORK(INDD), RWORK(INDE), M.value, W, IWORK(1), IWORK(INDISP),
          Z, LDZ, RWORK(INDRWK), IWORK(INDIWK), IFAIL, INFO);

      // Apply unitary matrix used in reduction to tridiagonal
      // form to eigenvectors returned by ZSTEIN.

      for (J = 1; J <= M.value; J++) {
        zcopy(N, Z(1, J).asArray(), 1, WORK(1), 1);
        zgemv('N', N, N, Complex.one, Q, LDQ, WORK, 1, Complex.zero,
            Z(1, J).asArray(), 1);
      }
    }
  }

  // If eigenvalues are not in order, then sort them, along with
  // eigenvectors.

  if (WANTZ) {
    for (J = 1; J <= M.value - 1; J++) {
      I = 0;
      TMP1 = W[J];
      for (JJ = J + 1; JJ <= M.value; JJ++) {
        if (W[JJ] < TMP1) {
          I = JJ;
          TMP1 = W[JJ];
        }
      }

      if (I != 0) {
        ITMP1 = IWORK[1 + I - 1];
        W[I] = W[J];
        IWORK[1 + I - 1] = IWORK[1 + J - 1];
        W[J] = TMP1;
        IWORK[1 + J - 1] = ITMP1;
        zswap(N, Z(1, I).asArray(), 1, Z(1, J).asArray(), 1);
        if (INFO.value != 0) {
          ITMP1 = IFAIL[I];
          IFAIL[I] = IFAIL[J];
          IFAIL[J] = ITMP1;
        }
      }
    }
  }
}
