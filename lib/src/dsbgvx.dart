// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/blas/dswap.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dpbstf.dart';
import 'package:dart_lapack/src/dsbgst.dart';
import 'package:dart_lapack/src/dsbtrd.dart';
import 'package:dart_lapack/src/dstebz.dart';
import 'package:dart_lapack/src/dstein.dart';
import 'package:dart_lapack/src/dsteqr.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsbgvx(
  final String JOBZ,
  final String RANGE,
  final String UPLO,
  final int N,
  final int KA,
  final int KB,
  final Matrix<double> AB_,
  final int LDAB,
  final Matrix<double> BB_,
  final int LDBB,
  final Matrix<double> Q_,
  final int LDQ,
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
  final AB = AB_.having(ld: LDAB);
  final BB = BB_.having(ld: LDBB);
  final Q = Q_.having(ld: LDQ);
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final IFAIL = IFAIL_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool ALLEIG, INDEIG, TEST, UPPER, VALEIG, WANTZ;
  String ORDER, VECT;
  int I, INDD, INDE, INDEE, INDISP, INDIWO, INDWRK, ITMP1, J, JJ;
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
    xerbla('DSBGVX', -INFO.value);
    return;
  }

  // Quick return if possible

  M.value = 0;
  if (N == 0) return;

  // Form a split Cholesky factorization of B.

  dpbstf(UPLO, N, KB, BB, LDBB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem.

  dsbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, WORK, IINFO);

  // Reduce symmetric band matrix to tridiagonal form.

  INDD = 1;
  INDE = INDD + N;
  INDWRK = INDE + N;
  if (WANTZ) {
    VECT = 'U';
  } else {
    VECT = 'N';
  }
  dsbtrd(VECT, UPLO, N, KA, AB, LDAB, WORK(INDD), WORK(INDE), Q, LDQ,
      WORK(INDWRK), IINFO);

  // If all eigenvalues are desired and ABSTOL is less than or equal
  // to zero, then call DSTERF or SSTEQR.  If this fails for some
  // eigenvalue, then try DSTEBZ.

  TEST = false;
  if (INDEIG) {
    if (IL == 1 && IU == N) {
      TEST = true;
    }
  }
  while (true) {
    if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
      dcopy(N, WORK(INDD), 1, W, 1);
      INDEE = INDWRK + 2 * N;
      dcopy(N - 1, WORK(INDE), 1, WORK(INDEE), 1);
      if (!WANTZ) {
        dsterf(N, W, WORK(INDEE), INFO);
      } else {
        dlacpy('A', N, N, Q, LDQ, Z, LDZ);
        dsteqr(JOBZ, N, W, WORK(INDEE), Z, LDZ, WORK(INDWRK), INFO);
        if (INFO.value == 0) {
          for (I = 1; I <= N; I++) {
            IFAIL[I] = 0;
          }
        }
      }
      if (INFO.value == 0) {
        M.value = N;
        break;
      }
      INFO.value = 0;
    }

    // Otherwise, call DSTEBZ and, if eigenvectors are desired,
    // call DSTEIN.

    if (WANTZ) {
      ORDER = 'B';
    } else {
      ORDER = 'E';
    }
    INDISP = 1 + N;
    INDIWO = INDISP + N;
    dstebz(RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, WORK(INDD), WORK(INDE), M,
        NSPLIT, W, IWORK(1), IWORK(INDISP), WORK(INDWRK), IWORK(INDIWO), INFO);

    if (WANTZ) {
      dstein(N, WORK(INDD), WORK(INDE), M.value, W, IWORK(1), IWORK(INDISP), Z,
          LDZ, WORK(INDWRK), IWORK(INDIWO), IFAIL, INFO);

      // Apply transformation matrix used in reduction to tridiagonal
      // form to eigenvectors returned by DSTEIN.

      for (J = 1; J <= M.value; J++) {
        dcopy(N, Z(1, J).asArray(), 1, WORK(1), 1);
        dgemv('N', N, N, ONE, Q, LDQ, WORK, 1, ZERO, Z(1, J).asArray(), 1);
      }
    }

    break;
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
        dswap(N, Z(1, I).asArray(), 1, Z(1, J).asArray(), 1);
        if (INFO.value != 0) {
          ITMP1 = IFAIL[I];
          IFAIL[I] = IFAIL[J];
          IFAIL[J] = ITMP1;
        }
      }
    }
  }
}
