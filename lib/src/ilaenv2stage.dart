// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/iparam2stage.dart';

int Function(
  int ISPEC,
  String NAME,
  String OPTS,
  int N1,
  int N2,
  int N3,
  int N4,
) ilaenv2stage = _ilaenv2stage;

int _ilaenv2stage(
  final int ISPEC,
  final String NAME,
  final String OPTS,
  final int N1,
  final int N2,
  final int N3,
  final int N4,
) {
  int IISPEC;
  switch (ISPEC) {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:

      // 2stage eigenvalues and SVD or related subroutines.

      IISPEC = 16 + ISPEC;
      return iparam2stage(IISPEC, NAME, OPTS, N1, N2, N3, N4);
    default:
      // Invalid value for ISPEC

      return -1;
  }
}
