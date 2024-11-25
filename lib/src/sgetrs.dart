// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgetrs.dart';
import 'package:dart_lapack/src/matrix.dart';

void sgetrs(
  final String TRANS,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
  dgetrs(TRANS, N, NRHS, A_, LDA, IPIV_, B_, LDB, INFO);
}
