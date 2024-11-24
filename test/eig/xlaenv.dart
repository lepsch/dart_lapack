// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'common.dart';

void xlaenv(final int ISPEC, final int NVALUE) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  if (ISPEC >= 1 && ISPEC <= 16) {
    claenv.IPARMS[ISPEC] = NVALUE;
  }
}
