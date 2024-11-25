// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zdotc.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgecon.dart';
import 'package:dart_lapack/src/zgesc2.dart';
import 'package:dart_lapack/src/zlassq.dart';
import 'package:dart_lapack/src/zlaswp.dart';

void zlatdf(
  final int IJOB,
  final int N,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> RHS_,
  final Box<double> RDSUM,
  final Box<double> RDSCAL,
  final Array<int> IPIV_,
  final Array<int> JPIV_,
) {
  final Z = Z_.having(ld: LDZ);
  final RHS = RHS_.having();
  final IPIV = IPIV_.having();
  final JPIV = JPIV_.having();

  const MAXDIM = 2;
  const ZERO = 0.0, ONE = 1.0;
  int I, J, K;
  double SMINU, SPLUS;
  Complex BM, BP, PMONE, TEMP;
  final RWORK = Array<double>(MAXDIM);
  final WORK = Array<Complex>(4 * MAXDIM),
      XM = Array<Complex>(MAXDIM),
      XP = Array<Complex>(MAXDIM);
  final INFO = Box(0);
  final RTEMP = Box(0.0), SCALE = Box(0.0);

  if (IJOB != 2) {
    // Apply permutations IPIV to RHS

    zlaswp(1, RHS.asMatrix(LDZ), LDZ, 1, N - 1, IPIV, 1);

    // Solve for L-part choosing RHS either to +1 or -1.

    PMONE = -Complex.one;
    for (J = 1; J <= N - 1; J++) {
      BP = RHS[J] + Complex.one;
      BM = RHS[J] - Complex.one;
      SPLUS = ONE;

      // Look-ahead for L- part RHS(1:N-1) = +-1
      // SPLUS and SMIN computed more efficiently than in BSOLVE[1].

      SPLUS +=
          zdotc(N - J, Z(J + 1, J).asArray(), 1, Z(J + 1, J).asArray(), 1).real;
      SMINU = zdotc(N - J, Z(J + 1, J).asArray(), 1, RHS(J + 1), 1).real;
      SPLUS *= RHS[J].real;
      if (SPLUS > SMINU) {
        RHS[J] = BP;
      } else if (SMINU > SPLUS) {
        RHS[J] = BM;
      } else {
        // In this case the updating sums are equal and we can
        // choose RHS(J) +1 or -1. The first time this happens we
        // choose -1, thereafter +1. This is a simple way to get
        // good estimates of matrices like Byers well-known example
        // (see [1]). (Not done in BSOLVE.)

        RHS[J] += PMONE;
        PMONE = Complex.one;
      }

      // Compute the remaining r.h.s.

      TEMP = -RHS[J];
      zaxpy(N - J, TEMP, Z(J + 1, J).asArray(), 1, RHS(J + 1), 1);
    }

    // Solve for U- part, lockahead for RHS(N) = +-1. This is not done
    // In BSOLVE and will hopefully give us a better estimate because
    // any ill-conditioning of the original matrix is transferred to U
    // and not to L. U(N, N) is an approximation to sigma_min(LU).

    zcopy(N - 1, RHS, 1, WORK, 1);
    WORK[N] = RHS[N] + Complex.one;
    RHS[N] -= Complex.one;
    SPLUS = ZERO;
    SMINU = ZERO;
    for (I = N; I >= 1; I--) {
      TEMP = Complex.one / Z[I][I];
      WORK[I] *= TEMP;
      RHS[I] *= TEMP;
      for (K = I + 1; K <= N; K++) {
        WORK[I] -= WORK[K] * (Z[I][K] * TEMP);
        RHS[I] -= RHS[K] * (Z[I][K] * TEMP);
      }
      SPLUS += WORK[I].abs();
      SMINU += RHS[I].abs();
    }
    if (SPLUS > SMINU) zcopy(N, WORK, 1, RHS, 1);

    // Apply the permutations JPIV to the computed solution (RHS)

    zlaswp(1, RHS.asMatrix(LDZ), LDZ, 1, N - 1, JPIV, -1);

    // Compute the sum of squares

    zlassq(N, RHS, 1, RDSCAL, RDSUM);
    return;
  }

  // ENTRY IJOB = 2

  // Compute approximate nullvector XM of Z

  zgecon('I', N, Z, LDZ, ONE, RTEMP, WORK, RWORK, INFO);
  zcopy(N, WORK(N + 1), 1, XM, 1);

  // Compute RHS

  zlaswp(1, XM.asMatrix(LDZ), LDZ, 1, N - 1, IPIV, -1);
  TEMP = Complex.one / zdotc(N, XM, 1, XM, 1).sqrt();
  zscal(N, TEMP, XM, 1);
  zcopy(N, XM, 1, XP, 1);
  zaxpy(N, Complex.one, RHS, 1, XP, 1);
  zaxpy(N, -Complex.one, XM, 1, RHS, 1);
  zgesc2(N, Z, LDZ, RHS, IPIV, JPIV, SCALE);
  zgesc2(N, Z, LDZ, XP, IPIV, JPIV, SCALE);
  if (dzasum(N, XP, 1) > dzasum(N, RHS, 1)) zcopy(N, XP, 1, RHS, 1);

  // Compute the sum of squares

  zlassq(N, RHS, 1, RDSCAL, RDSUM);
}
