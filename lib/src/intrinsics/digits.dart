// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

const int Function(double _) digits = digits64;

int digits32(final double _) => 24; // Float 32 fraction bits
int digits64(final double _) => 53; // Float 64 fraction bits
