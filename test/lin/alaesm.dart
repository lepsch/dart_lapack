// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/nio.dart';

void alaesm(
  final String PATH,
  final bool OK,
  final Nout NOUT,
) {
  if (OK) {
    NOUT.println(' ${PATH.a3} routines passed the tests of the error exits');
  } else {
    NOUT.println(
        ' *** ${PATH.a3} routines failed the tests of the error exits ***');
  }
}
