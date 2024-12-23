// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlarf.dart';

void zupmtr(
  final String SIDE,
  final String UPLO,
  final String TRANS,
  final int M,
  final int N,
  final Array<Complex> AP,
  final Array<Complex> TAU,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  bool FORWRD, LEFT, NOTRAN, UPPER;
  int I, I1, I2, I3, IC = 0, II, JC = 0, MI = 0, NI = 0, NQ;
  Complex AII, TAUI;

  // Test the input arguments

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  NOTRAN = lsame(TRANS, 'N');
  UPPER = lsame(UPLO, 'U');

  // NQ is the order of Q

  if (LEFT) {
    NQ = M;
  } else {
    NQ = N;
  }
  if (!LEFT && !lsame(SIDE, 'R')) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (!NOTRAN && !lsame(TRANS, 'C')) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDC < max(1, M)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('ZUPMTR', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  if (UPPER) {
    // Q was determined by a call to ZHPTRD with UPLO = 'U'

    FORWRD = (LEFT && NOTRAN) || (!LEFT && !NOTRAN);

    if (FORWRD) {
      I1 = 1;
      I2 = NQ - 1;
      I3 = 1;
      II = 2;
    } else {
      I1 = NQ - 1;
      I2 = 1;
      I3 = -1;
      II = NQ * (NQ + 1) ~/ 2 - 1;
    }

    if (LEFT) {
      NI = N;
    } else {
      MI = M;
    }

    for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
      if (LEFT) {
        // H(i) or H(i)**H is applied to C(1:i,1:n)

        MI = I;
      } else {
        // H(i) or H(i)**H is applied to C(1:m,1:i)

        NI = I;
      }

      // Apply H(i) or H(i)**H

      if (NOTRAN) {
        TAUI = TAU[I];
      } else {
        TAUI = TAU[I].conjugate();
      }
      AII = AP[II];
      AP[II] = Complex.one;
      zlarf(SIDE, MI, NI, AP(II - I + 1), 1, TAUI, C, LDC, WORK);
      AP[II] = AII;

      if (FORWRD) {
        II += I + 2;
      } else {
        II -= I + 1;
      }
    }
  } else {
    // Q was determined by a call to ZHPTRD with UPLO = 'L'.

    FORWRD = (LEFT && !NOTRAN) || (!LEFT && NOTRAN);

    if (FORWRD) {
      I1 = 1;
      I2 = NQ - 1;
      I3 = 1;
      II = 2;
    } else {
      I1 = NQ - 1;
      I2 = 1;
      I3 = -1;
      II = NQ * (NQ + 1) ~/ 2 - 1;
    }

    if (LEFT) {
      NI = N;
      JC = 1;
    } else {
      MI = M;
      IC = 1;
    }

    for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
      AII = AP[II];
      AP[II] = Complex.one;
      if (LEFT) {
        // H(i) or H(i)**H is applied to C(i+1:m,1:n)

        MI = M - I;
        IC = I + 1;
      } else {
        // H(i) or H(i)**H is applied to C(1:m,i+1:n)

        NI = N - I;
        JC = I + 1;
      }

      // Apply H(i) or H(i)**H

      if (NOTRAN) {
        TAUI = TAU[I];
      } else {
        TAUI = TAU[I].conjugate();
      }
      zlarf(SIDE, MI, NI, AP(II), 1, TAUI, C(IC, JC), LDC, WORK);
      AP[II] = AII;

      if (FORWRD) {
        II += NQ - I + 1;
      } else {
        II -= NQ - I + 2;
      }
    }
  }
}
