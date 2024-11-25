// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

double zla_porpvgrw(
  final String UPLO,
  final int NCOLS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
  final Array<double> WORK_,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final WORK = WORK_.having();
  int I, J;
  double AMAX, UMAX, RPVGRW;
  bool UPPER;

  UPPER = lsame('Upper', UPLO);

  // DPOTRF will have factored only the NCOLSxNCOLS leading submatrix,
  // so we restrict the growth search to that submatrix and use only
  // the first 2*NCOLS workspace entries.

  RPVGRW = 1.0;
  for (I = 1; I <= 2 * NCOLS; I++) {
    WORK[I] = 0.0;
  }

  // Find the max magnitude entry of each column.

  if (UPPER) {
    for (J = 1; J <= NCOLS; J++) {
      for (I = 1; I <= J; I++) {
        WORK[NCOLS + J] = max(A[I][J].cabs1(), WORK[NCOLS + J]);
      }
    }
  } else {
    for (J = 1; J <= NCOLS; J++) {
      for (I = J; I <= NCOLS; I++) {
        WORK[NCOLS + J] = max(A[I][J].cabs1(), WORK[NCOLS + J]);
      }
    }
  }

  // Now find the max magnitude entry of each column of the factor in
  // AF.  No pivoting, so no permutations.

  if (lsame('Upper', UPLO)) {
    for (J = 1; J <= NCOLS; J++) {
      for (I = 1; I <= J; I++) {
        WORK[J] = max(AF[I][J].cabs1(), WORK[J]);
      }
    }
  } else {
    for (J = 1; J <= NCOLS; J++) {
      for (I = J; I <= NCOLS; I++) {
        WORK[J] = max(AF[I][J].cabs1(), WORK[J]);
      }
    }
  }

  // Compute the *inverse* of the max element growth factor.  Dividing
  // by zero would imply the largest entry of the factor's column is
  // zero.  Than can happen when either the column of A is zero or
  // massive pivots made the factor underflow to zero.  Neither counts
  // as growth in itself, so simply ignore terms with zero
  // denominators.

  if (lsame('Upper', UPLO)) {
    for (I = 1; I <= NCOLS; I++) {
      UMAX = WORK[I];
      AMAX = WORK[NCOLS + I];
      if (UMAX != 0.0) {
        RPVGRW = min(AMAX / UMAX, RPVGRW);
      }
    }
  } else {
    for (I = 1; I <= NCOLS; I++) {
      UMAX = WORK[I];
      AMAX = WORK[NCOLS + I];
      if (UMAX != 0.0) {
        RPVGRW = min(AMAX / UMAX, RPVGRW);
      }
    }
  }

  return RPVGRW;
}
