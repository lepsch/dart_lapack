// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dtrsm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgbtrs.dart';
import 'package:dart_lapack/src/dlaswp.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsytrs_aa_2stage(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TB_,
  final int LTB,
  final Array<int> IPIV_,
  final Array<int> IPIV2_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TB = TB_.having();
  final IPIV = IPIV_.having();
  final IPIV2 = IPIV2_.having();
  final B = B_.having(ld: LDB);
  const ONE = 1.0;
  int LDTB, NB;
  bool UPPER;

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LTB < (4 * N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('DSYTRS_AA_2STAGE', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  // Read NB and compute LDTB

  NB = TB[1].toInt();
  LDTB = LTB ~/ N;

  if (UPPER) {
    // Solve A*X = B, where A = U**T*T*U.

    if (N > NB) {
      // Pivot, P**T * B -> B

      dlaswp(NRHS, B, LDB, NB + 1, N, IPIV, 1);

      // Compute (U**T \ B) -> B    [ (U**T \P**T * B) ]

      dtrsm('L', 'U', 'T', 'U', N - NB, NRHS, ONE, A(1, NB + 1), LDA,
          B(NB + 1, 1), LDB);
    }

    // Compute T \ B -> B   [ T \ (U**T \P**T * B) ]

    dgbtrs('N', N, NB, NB, NRHS, TB.asMatrix(LDTB), LDTB, IPIV2, B, LDB, INFO);
    if (N > NB) {
      // Compute (U \ B) -> B   [ U \ (T \ (U**T \P**T * B) ) ]

      dtrsm('L', 'U', 'N', 'U', N - NB, NRHS, ONE, A(1, NB + 1), LDA,
          B(NB + 1, 1), LDB);

      // Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]

      dlaswp(NRHS, B, LDB, NB + 1, N, IPIV, -1);
    }
  } else {
    // Solve A*X = B, where A = L*T*L**T.

    if (N > NB) {
      // Pivot, P**T * B -> B

      dlaswp(NRHS, B, LDB, NB + 1, N, IPIV, 1);

      // Compute (L \ B) -> B    [ (L \P**T * B) ]

      dtrsm('L', 'L', 'N', 'U', N - NB, NRHS, ONE, A(NB + 1, 1), LDA,
          B(NB + 1, 1), LDB);
    }

    // Compute T \ B -> B   [ T \ (L \P**T * B) ]

    dgbtrs('N', N, NB, NB, NRHS, TB.asMatrix(LDTB), LDTB, IPIV2, B, LDB, INFO);
    if (N > NB) {
      // Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]

      dtrsm('L', 'L', 'T', 'U', N - NB, NRHS, ONE, A(NB + 1, 1), LDA,
          B(NB + 1, 1), LDB);

      // Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]

      dlaswp(NRHS, B, LDB, NB + 1, N, IPIV, -1);
    }
  }
}
