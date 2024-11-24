// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zgerc.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlarfg.dart';
import 'package:lapack/src/zlarnv.dart';
import 'package:lapack/src/zlaset.dart';

import 'dlatm1.dart';
import 'zlarge.dart';
import 'zlarnd.dart';
import 'zlatm1.dart';

void zlatme(
  final int N,
  final String DIST,
  final Array<int> ISEED_,
  final Array<Complex> D_,
  final int MODE,
  final double COND,
  final Complex DMAX,
  final String RSIGN,
  final String UPPER,
  final String SIM,
  final Array<double> DS_,
  final int MODES,
  final double CONDS,
  final int KL,
  final int KU,
  final double ANORM,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final DS = DS_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0;
  const ONE = 1.0;
  bool BADS;
  int I, IC, ICOLS, IDIST, IR, IROWS, IRSIGN, ISIM, IUPPER, J, JC, JCR;
  double RALPHA, TEMP;
  Complex ALPHA;
  final TEMPA = Array<double>(1);
  final IINFO = Box(0);
  final TAU = Box(Complex.zero), XNORMS = Box(Complex.zero);

  // 1)      Decode and Test the input parameters.
  //         Initialize flags & seed.

  INFO.value = 0;

  // Quick return if possible

  if (N == 0) return;

  // Decode DIST

  if (lsame(DIST, 'U')) {
    IDIST = 1;
  } else if (lsame(DIST, 'S')) {
    IDIST = 2;
  } else if (lsame(DIST, 'N')) {
    IDIST = 3;
  } else if (lsame(DIST, 'D')) {
    IDIST = 4;
  } else {
    IDIST = -1;
  }

  // Decode RSIGN

  if (lsame(RSIGN, 'T')) {
    IRSIGN = 1;
  } else if (lsame(RSIGN, 'F')) {
    IRSIGN = 0;
  } else {
    IRSIGN = -1;
  }

  // Decode UPPER

  if (lsame(UPPER, 'T')) {
    IUPPER = 1;
  } else if (lsame(UPPER, 'F')) {
    IUPPER = 0;
  } else {
    IUPPER = -1;
  }

  // Decode SIM

  if (lsame(SIM, 'T')) {
    ISIM = 1;
  } else if (lsame(SIM, 'F')) {
    ISIM = 0;
  } else {
    ISIM = -1;
  }

  // Check DS, if MODES=0 and ISIM=1

  BADS = false;
  if (MODES == 0 && ISIM == 1) {
    for (J = 1; J <= N; J++) {
      if (DS[J] == ZERO) BADS = true;
    }
  }

  // Set INFO if an error

  if (N < 0) {
    INFO.value = -1;
  } else if (IDIST == -1) {
    INFO.value = -2;
  } else if (MODE.abs() > 6) {
    INFO.value = -5;
  } else if ((MODE != 0 && MODE.abs() != 6) && COND < ONE) {
    INFO.value = -6;
  } else if (IRSIGN == -1) {
    INFO.value = -9;
  } else if (IUPPER == -1) {
    INFO.value = -10;
  } else if (ISIM == -1) {
    INFO.value = -11;
  } else if (BADS) {
    INFO.value = -12;
  } else if (ISIM == 1 && MODES.abs() > 5) {
    INFO.value = -13;
  } else if (ISIM == 1 && MODES != 0 && CONDS < ONE) {
    INFO.value = -14;
  } else if (KL < 1) {
    INFO.value = -15;
  } else if (KU < 1 || (KU < N - 1 && KL < N - 1)) {
    INFO.value = -16;
  } else if (LDA < max(1, N)) {
    INFO.value = -19;
  }

  if (INFO.value != 0) {
    xerbla('ZLATME', -INFO.value);
    return;
  }

  // Initialize random number generator

  for (I = 1; I <= 4; I++) {
    ISEED[I] = ISEED[I].abs() % 4096;
  }

  if ((ISEED[4] % 2) != 1) ISEED[4]++;

  // 2)      Set up diagonal of A

  // Compute D according to COND and MODE

  zlatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, N, IINFO);
  if (IINFO.value != 0) {
    INFO.value = 1;
    return;
  }
  if (MODE != 0 && MODE.abs() != 6) {
    // Scale by DMAX

    TEMP = D[1].abs();
    for (I = 2; I <= N; I++) {
      TEMP = max(TEMP, D[I].abs());
    }

    if (TEMP > ZERO) {
      ALPHA = DMAX / TEMP.toComplex();
    } else {
      INFO.value = 2;
      return;
    }

    zscal(N, ALPHA, D, 1);
  }

  zlaset('Full', N, N, Complex.zero, Complex.zero, A, LDA);
  zcopy(N, D, 1, A.asArray(), LDA + 1);

  // 3)      If UPPER='T', set upper triangle of A to random numbers.

  if (IUPPER != 0) {
    for (JC = 2; JC <= N; JC++) {
      zlarnv(IDIST, ISEED, JC - 1, A(1, JC).asArray());
    }
  }

  // 4)      If SIM='T', apply similarity transformation.

  // -1
  // Transform is  X A X  , where X = U S V, thus

  // it is  U S V A V' (1/S) U'

  if (ISIM != 0) {
    // Compute S (singular values of the eigenvector matrix)
    // according to CONDS and MODES

    dlatm1(MODES, CONDS, 0, 0, ISEED, DS, N, IINFO);
    if (IINFO.value != 0) {
      INFO.value = 3;
      return;
    }

    // Multiply by V and V'

    zlarge(N, A, LDA, ISEED, WORK, IINFO);
    if (IINFO.value != 0) {
      INFO.value = 4;
      return;
    }

    // Multiply by S and (1/S)

    for (J = 1; J <= N; J++) {
      zdscal(N, DS[J], A(J, 1).asArray(), LDA);
      if (DS[J] != ZERO) {
        zdscal(N, ONE / DS[J], A(1, J).asArray(), 1);
      } else {
        INFO.value = 5;
        return;
      }
    }

    // Multiply by U and U'

    zlarge(N, A, LDA, ISEED, WORK, IINFO);
    if (IINFO.value != 0) {
      INFO.value = 4;
      return;
    }
  }

  // 5)      Reduce the bandwidth.

  if (KL < N - 1) {
    // Reduce bandwidth -- kill column

    for (JCR = KL + 1; JCR <= N - 1; JCR++) {
      IC = JCR - KL;
      IROWS = N + 1 - JCR;
      ICOLS = N + KL - JCR;

      zcopy(IROWS, A(JCR, IC).asArray(), 1, WORK, 1);
      XNORMS.value = WORK[1];
      zlarfg(IROWS, XNORMS, WORK(2), 1, TAU);
      TAU.value = TAU.value.conjugate();
      WORK[1] = Complex.one;
      ALPHA = zlarnd(5, ISEED);

      zgemv('C', IROWS, ICOLS, Complex.one, A(JCR, IC + 1), LDA, WORK, 1,
          Complex.zero, WORK(IROWS + 1), 1);
      zgerc(IROWS, ICOLS, -TAU.value, WORK, 1, WORK(IROWS + 1), 1,
          A(JCR, IC + 1), LDA);

      zgemv('N', N, IROWS, Complex.one, A(1, JCR), LDA, WORK, 1, Complex.zero,
          WORK(IROWS + 1), 1);
      zgerc(N, IROWS, -TAU.value.conjugate(), WORK(IROWS + 1), 1, WORK, 1,
          A(1, JCR), LDA);

      A[JCR][IC] = XNORMS.value;
      zlaset('Full', IROWS - 1, 1, Complex.zero, Complex.zero, A(JCR + 1, IC),
          LDA);

      zscal(ICOLS + 1, ALPHA, A(JCR, IC).asArray(), LDA);
      zscal(N, ALPHA.conjugate(), A(1, JCR).asArray(), 1);
    }
  } else if (KU < N - 1) {
    // Reduce upper bandwidth -- kill a row at a time.

    for (JCR = KU + 1; JCR <= N - 1; JCR++) {
      IR = JCR - KU;
      IROWS = N + KU - JCR;
      ICOLS = N + 1 - JCR;

      zcopy(ICOLS, A(IR, JCR).asArray(), LDA, WORK, 1);
      XNORMS.value = WORK[1];
      zlarfg(ICOLS, XNORMS, WORK(2), 1, TAU);
      TAU.value = TAU.value.conjugate();
      WORK[1] = Complex.one;
      zlacgv(ICOLS - 1, WORK(2), 1);
      ALPHA = zlarnd(5, ISEED);

      zgemv('N', IROWS, ICOLS, Complex.one, A(IR + 1, JCR), LDA, WORK, 1,
          Complex.zero, WORK(ICOLS + 1), 1);
      zgerc(IROWS, ICOLS, -TAU.value, WORK(ICOLS + 1), 1, WORK, 1,
          A(IR + 1, JCR), LDA);

      zgemv('C', ICOLS, N, Complex.one, A(JCR, 1), LDA, WORK, 1, Complex.zero,
          WORK(ICOLS + 1), 1);
      zgerc(ICOLS, N, -TAU.value.conjugate(), WORK, 1, WORK(ICOLS + 1), 1,
          A(JCR, 1), LDA);

      A[IR][JCR] = XNORMS.value;
      zlaset('Full', 1, ICOLS - 1, Complex.zero, Complex.zero, A(IR, JCR + 1),
          LDA);

      zscal(IROWS + 1, ALPHA, A(IR, JCR).asArray(), 1);
      zscal(N, ALPHA.conjugate(), A(JCR, 1).asArray(), LDA);
    }
  }

  // Scale the matrix to have norm ANORM

  if (ANORM >= ZERO) {
    TEMP = zlange('M', N, N, A, LDA, TEMPA);
    if (TEMP > ZERO) {
      RALPHA = ANORM / TEMP;
      for (J = 1; J <= N; J++) {
        zdscal(N, RALPHA, A(1, J).asArray(), 1);
      }
    }
  }
}
