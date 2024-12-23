// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void dlamrg(
  final int N1,
  final int N2,
  final Array<double> A_,
  final int DTRD1,
  final int DTRD2,
  final Array<int> INDEX,
) {
  final A = A_.having();
  int I, IND1, IND2, N1SV, N2SV;

  N1SV = N1;
  N2SV = N2;
  if (DTRD1 > 0) {
    IND1 = 1;
  } else {
    IND1 = N1;
  }
  if (DTRD2 > 0) {
    IND2 = 1 + N1;
  } else {
    IND2 = N1 + N2;
  }
  I = 1;
  while (N1SV > 0 && N2SV > 0) {
    if (A[IND1] <= A[IND2]) {
      INDEX[I] = IND1;
      I++;
      IND1 += DTRD1;
      N1SV--;
    } else {
      INDEX[I] = IND2;
      I++;
      IND2 += DTRD2;
      N2SV--;
    }
  }
  // end while
  if (N1SV == 0) {
    for (N1SV = 1; N1SV <= N2SV; N1SV++) {
      INDEX[I] = IND2;
      I++;
      IND2 += DTRD2;
    }
  } else {
    // N2SV == 0
    for (N2SV = 1; N2SV <= N1SV; N2SV++) {
      INDEX[I] = IND1;
      I++;
      IND1 += DTRD1;
    }
  }
}
