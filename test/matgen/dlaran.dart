// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/matrix.dart';

double dlaran(final Array<int> ISEED_) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
  final ISEED = ISEED_.having();
  const M1 = 494, M2 = 322, M3 = 2508, M4 = 2549;
  const ONE = 1.0;
  const IPW2 = 4096, R = ONE / IPW2;
  int IT1, IT2, IT3, IT4;
  double RNDOUT;

  do {
    // multiply the seed by the multiplier modulo 2**48

    IT4 = ISEED[4] * M4;
    IT3 = IT4 ~/ IPW2;
    IT4 -= IPW2 * IT3;
    IT3 += ISEED[3] * M4 + ISEED[4] * M3;
    IT2 = IT3 ~/ IPW2;
    IT3 -= IPW2 * IT2;
    IT2 += ISEED[2] * M4 + ISEED[3] * M3 + ISEED[4] * M2;
    IT1 = IT2 ~/ IPW2;
    IT2 -= IPW2 * IT1;
    IT1 += ISEED[1] * M4 + ISEED[2] * M3 + ISEED[3] * M2 + ISEED[4] * M1;
    IT1 = IT1 % IPW2;

    // return updated seed

    ISEED[1] = IT1;
    ISEED[2] = IT2;
    ISEED[3] = IT3;
    ISEED[4] = IT4;

    // convert 48-bit integer to a real number in the interval (0,1)

    RNDOUT = R * (IT1 + R * (IT2 + R * (IT3 + R * IT4)));

    // If a real number has n bits of precision, and the first
    // n bits of the 48-bit integer above happen to be all 1 (which
    // will occur about once every 2**n calls), then DLARAN will
    // be rounded to exactly 1.0.
    // Since DLARAN is not supposed to return exactly 0.0 or 1.0
    // (and some callers of DLARAN, such as CLARND, depend on that),
    // the statistically correct thing to do in this situation is
    // simply to iterate again.
    // N.B. the case DLARAN = 0.0 should not be possible.
  } while (RNDOUT == 1.0);

  return RNDOUT;
}
