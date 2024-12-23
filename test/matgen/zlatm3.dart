// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'dlaran.dart';
import 'zlarnd.dart';

Complex zlatm3(
  final int M,
  final int N,
  final int I,
  final int J,
  final Box<int> ISUB,
  final Box<int> JSUB,
  final int KL,
  final int KU,
  final int IDIST,
  final Array<int> ISEED_,
  final Array<Complex> D_,
  final int IGRADE,
  final Array<Complex> DL_,
  final Array<Complex> DR_,
  final int IPVTNG,
  final Array<int> IWORK_,
  final double SPARSE,
) {
  final ISEED = ISEED_.having();
  final IWORK = IWORK_.having();
  final D = D_.having();
  final DL = DL_.having();
  final DR = DR_.having();
  const ZERO = 0.0;
  Complex CTEMP;

  // Check for I and J in range

  if (I < 1 || I > M || J < 1 || J > N) {
    ISUB.value = I;
    JSUB.value = J;
    return Complex.zero;
  }

  // Compute subscripts depending on IPVTNG

  if (IPVTNG == 0) {
    ISUB.value = I;
    JSUB.value = J;
  } else if (IPVTNG == 1) {
    ISUB.value = IWORK[I];
    JSUB.value = J;
  } else if (IPVTNG == 2) {
    ISUB.value = I;
    JSUB.value = IWORK[J];
  } else if (IPVTNG == 3) {
    ISUB.value = IWORK[I];
    JSUB.value = IWORK[J];
  }

  // Check for banding

  if (JSUB.value > ISUB.value + KU || JSUB.value < ISUB.value - KL) {
    return Complex.zero;
  }

  // Check for sparsity

  if (SPARSE > ZERO) {
    if (dlaran(ISEED) < SPARSE) {
      return Complex.zero;
    }
  }

  // Compute entry and grade it according to IGRADE

  if (I == J) {
    CTEMP = D[I];
  } else {
    CTEMP = zlarnd(IDIST, ISEED);
  }
  if (IGRADE == 1) {
    CTEMP *= DL[I];
  } else if (IGRADE == 2) {
    CTEMP *= DR[J];
  } else if (IGRADE == 3) {
    CTEMP *= DL[I] * DR[J];
  } else if (IGRADE == 4 && I != J) {
    CTEMP *= DL[I] / DL[J];
  } else if (IGRADE == 5) {
    CTEMP *= DL[I] * DL[J].conjugate();
  } else if (IGRADE == 6) {
    CTEMP *= DL[I] * DL[J];
  }
  return CTEMP;
}
