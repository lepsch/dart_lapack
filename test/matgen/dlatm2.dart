// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

import 'dlaran.dart';
import 'dlarnd.dart';

double dlatm2(
  final int M,
  final int N,
  final int I,
  final int J,
  final int KL,
  final int KU,
  final int IDIST,
  final Array<int> ISEED_,
  final Array<double> D_,
  final int IGRADE,
  final Array<double> DL_,
  final Array<double> DR_,
  final int IPVTNG,
  final Array<int> IWORK_,
  final double SPARSE,
) {
  final ISEED = ISEED_.having();
  final D = D_.having();
  final DL = DL_.having();
  final DR = DR_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0;
  int ISUB = 0, JSUB = 0;
  double TEMP;

  // Check for I and J in range
  if (I < 1 || I > M || J < 1 || J > N) return ZERO;

  // Check for banding
  if (J > I + KU || J < I - KL) return ZERO;

  // Check for sparsity
  if (SPARSE > ZERO) {
    if (dlaran(ISEED) < SPARSE) {
      return ZERO;
    }
  }

  // Compute subscripts depending on IPVTNG
  if (IPVTNG == 0) {
    ISUB = I;
    JSUB = J;
  } else if (IPVTNG == 1) {
    ISUB = IWORK[I];
    JSUB = J;
  } else if (IPVTNG == 2) {
    ISUB = I;
    JSUB = IWORK[J];
  } else if (IPVTNG == 3) {
    ISUB = IWORK[I];
    JSUB = IWORK[J];
  }

  // Compute entry and grade it according to IGRADE
  if (ISUB == JSUB) {
    TEMP = D[ISUB];
  } else {
    TEMP = dlarnd(IDIST, ISEED);
  }
  if (IGRADE == 1) {
    TEMP *= DL[ISUB];
  } else if (IGRADE == 2) {
    TEMP *= DR[JSUB];
  } else if (IGRADE == 3) {
    TEMP *= DL[ISUB] * DR[JSUB];
  } else if (IGRADE == 4 && ISUB != JSUB) {
    TEMP *= DL[ISUB] / DL[JSUB];
  } else if (IGRADE == 5) {
    TEMP *= DL[ISUB] * DL[JSUB];
  }
  return TEMP;
}
