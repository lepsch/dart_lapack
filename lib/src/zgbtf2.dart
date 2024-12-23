// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/blas/zgeru.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zgbtf2(
  final int M,
  final int N,
  final int KL,
  final int KU,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
  final AB = AB_.having(ld: LDAB);
  final IPIV = IPIV_.having();
  int I, J, JP, JU, KM, KV;

  // KV is the number of superdiagonals in the factor U, allowing for
  // fill-in.

  KV = KU + KL;

  // Test the input parameters.

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KL < 0) {
    INFO.value = -3;
  } else if (KU < 0) {
    INFO.value = -4;
  } else if (LDAB < KL + KV + 1) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZGBTF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Gaussian elimination with partial pivoting

  // Set fill-in elements in columns KU+2 to KV to zero.

  for (J = KU + 2; J <= min(KV, N); J++) {
    for (I = KV - J + 2; I <= KL; I++) {
      AB[I][J] = Complex.zero;
    }
  }

  // JU is the index of the last column affected by the current stage
  // of the factorization.

  JU = 1;

  for (J = 1; J <= min(M, N); J++) {
    // Set fill-in elements in column J+KV to zero.

    if (J + KV <= N) {
      for (I = 1; I <= KL; I++) {
        AB[I][J + KV] = Complex.zero;
      }
    }

    // Find pivot and test for singularity. KM is the number of
    // subdiagonal elements in the current column.

    KM = min(KL, M - J);
    JP = izamax(KM + 1, AB(KV + 1, J).asArray(), 1);
    IPIV[J] = JP + J - 1;
    if (AB[KV + JP][J] != Complex.zero) {
      JU = max(JU, min(J + KU + JP - 1, N));

      // Apply interchange to columns J to JU.

      if (JP != 1) {
        zswap(JU - J + 1, AB(KV + JP, J).asArray(), LDAB - 1,
            AB(KV + 1, J).asArray(), LDAB - 1);
      }
      if (KM > 0) {
        // Compute multipliers.

        zscal(KM, Complex.one / AB[KV + 1][J], AB(KV + 2, J).asArray(), 1);

        // Update trailing submatrix within the band.

        if (JU > J) {
          zgeru(KM, JU - J, -Complex.one, AB(KV + 2, J).asArray(), 1,
              AB(KV, J + 1).asArray(), LDAB - 1, AB(KV + 1, J + 1), LDAB - 1);
        }
      }
    } else {
      // If pivot is zero, set INFO to the index of the pivot
      // unless a zero pivot has already been found.

      if (INFO.value == 0) INFO.value = J;
    }
  }
}
