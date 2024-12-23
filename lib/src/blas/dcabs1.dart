// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/blas.dart';

double dcabs1(final Complex Z) {
  return Z.real.abs() + Z.imaginary.abs();
}
