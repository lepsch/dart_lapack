// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'common.dart';

void xlaenv(final int ISPEC, final int NVALUE) {
  if (ISPEC >= 1 && ISPEC <= 9) {
    claenv.IPARMS[ISPEC] = NVALUE;
  }
}
