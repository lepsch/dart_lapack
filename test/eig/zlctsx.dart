// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';

import 'common.dart';

bool zlctsx(final Complex ALPHA, final Complex BETA) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  bool result;
  if (mn.FS) {
    mn.I++;
    if (mn.I <= mn.M) {
      result = false;
    } else {
      result = true;
    }
    if (mn.I == mn.MPLUSN) {
      mn.FS = false;
      mn.I = 0;
    }
  } else {
    mn.I++;
    if (mn.I <= mn.N) {
      result = true;
    } else {
      result = false;
    }
    if (mn.I == mn.MPLUSN) {
      mn.FS = true;
      mn.I = 0;
    }
  }

  // IF( BETA == CZERO ) THEN
  // ZLCTSX = ( ALPHA > ZERO )
  // ELSE
  // ZLCTSX = ( (ALPHA/BETA) > ZERO )
  // END IF

  return result;
}
