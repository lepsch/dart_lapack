// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';

void ztbt05(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int KD,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> XACT_,
  final int LDXACT,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<double> RESLTS_,
) {
  final AB = AB_.having(ld: LDAB);
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
  final UPPER = lsame(UPLO, 'U');
  final NOTRAN = lsame(TRANS, 'N');
  final UNIT = lsame(DIAG, 'U');
  final NZ = min(KD, N - 1) + 1;

  // Test 1:  Compute the maximum of
  //    norm(X - XACT) / ( norm(X) * FERR )
  // over all the vectors X and XACT using the infinity-norm.

  var ERRBND = ZERO;
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
  RESLTS[1] = ERRBND;

  // Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
  // (*) = NZ*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

  final IFU = UNIT ? 1 : 0;
  for (var K = 1; K <= NRHS; K++) {
    var AXBI = ZERO;
    for (var I = 1; I <= N; I++) {
      var TMP = B[I][K].cabs1();
      if (UPPER) {
        if (!NOTRAN) {
          for (var J = max(I - KD, 1); J <= I - IFU; J++) {
            TMP += AB[KD + 1 - I + J][I].cabs1() * X[J][K].cabs1();
          }
          if (UNIT) TMP += X[I][K].cabs1();
        } else {
          if (UNIT) TMP += X[I][K].cabs1();
          for (var J = I + IFU; J <= min(I + KD, N); J++) {
            TMP += AB[KD + 1 + I - J][J].cabs1() * X[J][K].cabs1();
          }
        }
      } else {
        if (NOTRAN) {
          for (var J = max(I - KD, 1); J <= I - IFU; J++) {
            TMP += AB[1 + I - J][J].cabs1() * X[J][K].cabs1();
          }
          if (UNIT) TMP += X[I][K].cabs1();
        } else {
          if (UNIT) TMP += X[I][K].cabs1();
          for (var J = I + IFU; J <= min(I + KD, N); J++) {
            TMP += AB[1 + J - I][I].cabs1() * X[J][K].cabs1();
          }
        }
      }
      if (I == 1) {
        AXBI = TMP;
      } else {
        AXBI = min(AXBI, TMP);
      }
    }
    final TMP = BERR[K] / (NZ * EPS + NZ * UNFL / max(AXBI, NZ * UNFL));
    if (K == 1) {
      RESLTS[2] = TMP;
    } else {
      RESLTS[2] = max(RESLTS[2], TMP);
    }
  }
}
