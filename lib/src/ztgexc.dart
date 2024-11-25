// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/ztgex2.dart';

void ztgexc(
  final bool WANTQ,
  final bool WANTZ,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final int IFST,
  final Box<int> ILST,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  int HERE = 0;

  // Decode and test input arguments.
  INFO.value = 0;
  if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  } else if (LDQ < 1 || WANTQ && (LDQ < max(1, N))) {
    INFO.value = -9;
  } else if (LDZ < 1 || WANTZ && (LDZ < max(1, N))) {
    INFO.value = -11;
  } else if (IFST < 1 || IFST > N) {
    INFO.value = -12;
  } else if (ILST.value < 1 || ILST.value > N) {
    INFO.value = -13;
  }
  if (INFO.value != 0) {
    xerbla('ZTGEXC', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 1) return;
  if (IFST == ILST.value) return;

  if (IFST < ILST.value) {
    HERE = IFST;

    do {
      // Swap with next one below

      ztgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, INFO);
      if (INFO.value != 0) {
        ILST.value = HERE;
        return;
      }
      HERE++;
    } while (HERE < ILST.value);
    HERE--;
  } else {
    HERE = IFST - 1;

    do {
      // Swap with next one above

      ztgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, INFO);
      if (INFO.value != 0) {
        ILST.value = HERE;
        return;
      }
      HERE--;
    } while (HERE >= ILST.value);
    HERE++;
  }
  ILST.value = HERE;
}
