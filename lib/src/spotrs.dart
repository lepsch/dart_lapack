// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/dpotrs.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void spotrs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  dpotrs(UPLO, NRHS, N, A_, LDA, B_, LDB, INFO);
}
