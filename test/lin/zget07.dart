// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';

void zget07(
  final String TRANS,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> XACT_,
  final int LDXACT,
  final Array<double> FERR_,
  final bool CHKFERR,
  final Array<double> BERR_,
  final Array<double> RESLTS_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final XACT = XACT_.having(ld: LDXACT);
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final RESLTS = RESLTS_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0 or NRHS = 0.

  if (N <= 0 || NRHS <= 0) {
    RESLTS[1] = ZERO;
    RESLTS[2] = ZERO;
    return;
  }

  final EPS = dlamch('Epsilon');
  final UNFL = dlamch('Safe minimum');
  final OVFL = ONE / UNFL;
  final NOTRAN = lsame(TRANS, 'N');

  // Test 1:  Compute the maximum of
  //    norm(X - XACT) / ( norm(X) * FERR )
  // over all the vectors X and XACT using the infinity-norm.

  var ERRBND = ZERO;
  if (CHKFERR) {
    for (var J = 1; J <= NRHS; J++) {
      final IMAX = izamax(N, X(1, J).asArray(), 1);
      final XNORM = max(X[IMAX][J].cabs1(), UNFL);
      var DIFF = ZERO;
      for (var I = 1; I <= N; I++) {
        DIFF = max(DIFF, (X[I][J] - XACT[I][J]).cabs1());
      }

      if (XNORM > ONE) {
        //
      } else if (DIFF <= OVFL * XNORM) {
        //
      } else {
        ERRBND = ONE / EPS;
        continue;
      }

      if (DIFF / XNORM <= FERR[J]) {
        ERRBND = max(ERRBND, (DIFF / XNORM) / FERR[J]);
      } else {
        ERRBND = ONE / EPS;
      }
    }
  }
  RESLTS[1] = ERRBND;

  // Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
  // (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )

  for (var K = 1; K <= NRHS; K++) {
    var AXBI = ZERO;
    for (var I = 1; I <= N; I++) {
      var TMP = B[I][K].cabs1();
      if (NOTRAN) {
        for (var J = 1; J <= N; J++) {
          TMP += A[I][J].cabs1() * X[J][K].cabs1();
        }
      } else {
        for (var J = 1; J <= N; J++) {
          TMP += A[J][I].cabs1() * X[J][K].cabs1();
        }
      }
      if (I == 1) {
        AXBI = TMP;
      } else {
        AXBI = min(AXBI, TMP);
      }
    }
    final TMP =
        BERR[K] / ((N + 1) * EPS + (N + 1) * UNFL / max(AXBI, (N + 1) * UNFL));
    if (K == 1) {
      RESLTS[2] = TMP;
    } else {
      RESLTS[2] = max(RESLTS[2], TMP);
    }
  }
}
