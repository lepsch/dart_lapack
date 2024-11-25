// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'common.dart';

bool dlctsx(
  double AR,
  double AI,
  double BETA,
) {
  final bool result;

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

  // IF( AR/BETA > 0.0 )THEN
  // DLCTSX = true;
  // ELSE
  // DLCTSX = false;
  // END IF

  return result;
}
