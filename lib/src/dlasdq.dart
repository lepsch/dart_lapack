// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dbdsqr.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlasr.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlasdq(
  final String UPLO,
  final int SQRE,
  final int N,
  final int NCVT,
  final int NRU,
  final int NCC,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> VT_,
  final int LDVT,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final VT = VT_.having(ld: LDVT);
  final U = U_.having(ld: LDU);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  const ZERO = 0.0;
  bool ROTATE;
  int I, ISUB, IUPLO, J, NP1, SQRE1;
  double SMIN;
  final CS = Box(0.0), R = Box(0.0), SN = Box(0.0);

  // Test the input parameters.

  INFO.value = 0;
  IUPLO = 0;
  if (lsame(UPLO, 'U')) IUPLO = 1;
  if (lsame(UPLO, 'L')) IUPLO = 2;
  if (IUPLO == 0) {
    INFO.value = -1;
  } else if ((SQRE < 0) || (SQRE > 1)) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NCVT < 0) {
    INFO.value = -4;
  } else if (NRU < 0) {
    INFO.value = -5;
  } else if (NCC < 0) {
    INFO.value = -6;
  } else if ((NCVT == 0 && LDVT < 1) || (NCVT > 0 && LDVT < max(1, N))) {
    INFO.value = -10;
  } else if (LDU < max(1, NRU)) {
    INFO.value = -12;
  } else if ((NCC == 0 && LDC < 1) || (NCC > 0 && LDC < max(1, N))) {
    INFO.value = -14;
  }
  if (INFO.value != 0) {
    xerbla('DLASDQ', -INFO.value);
    return;
  }
  if (N == 0) return;

  // ROTATE is true if any singular vectors desired, false otherwise

  ROTATE = (NCVT > 0) || (NRU > 0) || (NCC > 0);
  NP1 = N + 1;
  SQRE1 = SQRE;

  // If matrix non-square upper bidiagonal, rotate to be lower
  // bidiagonal.  The rotations are on the right.

  if ((IUPLO == 1) && (SQRE1 == 1)) {
    for (I = 1; I <= N - 1; I++) {
      dlartg(D[I], E[I], CS, SN, R);
      D[I] = R.value;
      E[I] = SN.value * D[I + 1];
      D[I + 1] = CS.value * D[I + 1];
      if (ROTATE) {
        WORK[I] = CS.value;
        WORK[N + I] = SN.value;
      }
    }
    dlartg(D[N], E[N], CS, SN, R);
    D[N] = R.value;
    E[N] = ZERO;
    if (ROTATE) {
      WORK[N] = CS.value;
      WORK[N + N] = SN.value;
    }
    IUPLO = 2;
    SQRE1 = 0;

    // Update singular vectors if desired.

    if (NCVT > 0) dlasr('L', 'V', 'F', NP1, NCVT, WORK(1), WORK(NP1), VT, LDVT);
  }

  // If matrix lower bidiagonal, rotate to be upper bidiagonal
  // by applying Givens rotations on the left.

  if (IUPLO == 2) {
    for (I = 1; I <= N - 1; I++) {
      dlartg(D[I], E[I], CS, SN, R);
      D[I] = R.value;
      E[I] = SN.value * D[I + 1];
      D[I + 1] = CS.value * D[I + 1];
      if (ROTATE) {
        WORK[I] = CS.value;
        WORK[N + I] = SN.value;
      }
    }

    // If matrix (N+1)-by-N lower bidiagonal, one additional
    // rotation is needed.

    if (SQRE1 == 1) {
      dlartg(D[N], E[N], CS, SN, R);
      D[N] = R.value;
      if (ROTATE) {
        WORK[N] = CS.value;
        WORK[N + N] = SN.value;
      }
    }

    // Update singular vectors if desired.

    if (NRU > 0) {
      if (SQRE1 == 0) {
        dlasr('R', 'V', 'F', NRU, N, WORK(1), WORK(NP1), U, LDU);
      } else {
        dlasr('R', 'V', 'F', NRU, NP1, WORK(1), WORK(NP1), U, LDU);
      }
    }
    if (NCC > 0) {
      if (SQRE1 == 0) {
        dlasr('L', 'V', 'F', N, NCC, WORK(1), WORK(NP1), C, LDC);
      } else {
        dlasr('L', 'V', 'F', NP1, NCC, WORK(1), WORK(NP1), C, LDC);
      }
    }
  }

  // Call DBDSQR to compute the SVD of the reduced real
  // N-by-N upper bidiagonal matrix.

  dbdsqr('U', N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO);

  // Sort the singular values into ascending order (insertion sort on
  // singular values, but only one transposition per singular vector)

  for (I = 1; I <= N; I++) {
    // Scan for smallest D(I).

    ISUB = I;
    SMIN = D[I];
    for (J = I + 1; J <= N; J++) {
      if (D[J] < SMIN) {
        ISUB = J;
        SMIN = D[J];
      }
    }
    if (ISUB != I) {
      // Swap singular values and vectors.

      D[ISUB] = D[I];
      D[I] = SMIN;
      if (NCVT > 0) {
        dswap(NCVT, VT(ISUB, 1).asArray(), LDVT, VT(I, 1).asArray(), LDVT);
      }
      if (NRU > 0) dswap(NRU, U(1, ISUB).asArray(), 1, U(1, I).asArray(), 1);
      if (NCC > 0) {
        dswap(NCC, C(ISUB, 1).asArray(), LDC, C(I, 1).asArray(), LDC);
      }
    }
  }
}
