// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

bool lsame(final String CA, final String CB) {
  return CA.isNotEmpty &&
      CB.isNotEmpty &&
      CA.toUpperCase()[0] == CB.toUpperCase()[0];
}
