// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/ztrmv.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlarft(
  final String DIRECT,
  final String STOREV,
  final int N,
  final int K,
  final Matrix<Complex> V_,
  final int LDV,
  final Array<Complex> TAU_,
  final Matrix<Complex> T_,
  final int LDT,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.having(ld: LDV);
  final TAU = TAU_.having();
  final T = T_.having(ld: LDT);
  int I, J, PREVLASTV, LASTV;

  // Quick return if possible

  if (N == 0) return;

  if (lsame(DIRECT, 'F')) {
    PREVLASTV = N;
    for (I = 1; I <= K; I++) {
      PREVLASTV = max(PREVLASTV, I);
      if (TAU[I] == Complex.zero) {
        // H(i)  =  I

        for (J = 1; J <= I; J++) {
          T[J][I] = Complex.zero;
        }
      } else {
        // general case

        if (lsame(STOREV, 'C')) {
          // Skip any trailing zeros.
          for (LASTV = N; LASTV >= I + 1; LASTV--) {
            if (V[LASTV][I] != Complex.zero) break;
          }
          for (J = 1; J <= I - 1; J++) {
            T[J][I] = -TAU[I] * V[I][J].conjugate();
          }
          J = min(LASTV, PREVLASTV);

          // T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)

          zgemv('Conjugate transpose', J - I, I - 1, -TAU[I], V(I + 1, 1), LDV,
              V(I + 1, I).asArray(), 1, Complex.one, T(1, I).asArray(), 1);
        } else {
          // Skip any trailing zeros.
          for (LASTV = N; LASTV >= I + 1; LASTV--) {
            if (V[I][LASTV] != Complex.zero) break;
          }
          for (J = 1; J <= I - 1; J++) {
            T[J][I] = -TAU[I] * V[J][I];
          }
          J = min(LASTV, PREVLASTV);

          // T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H

          zgemm('N', 'C', I - 1, 1, J - I, -TAU[I], V(1, I + 1), LDV,
              V(I, I + 1), LDV, Complex.one, T(1, I), LDT);
        }

        // T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)

        ztrmv('Upper', 'No transpose', 'Non-unit', I - 1, T, LDT,
            T(1, I).asArray(), 1);
        T[I][I] = TAU[I];
        if (I > 1) {
          PREVLASTV = max(PREVLASTV, LASTV);
        } else {
          PREVLASTV = LASTV;
        }
      }
    }
  } else {
    PREVLASTV = 1;
    for (I = K; I >= 1; I--) {
      if (TAU[I] == Complex.zero) {
        // H(i)  =  I

        for (J = I; J <= K; J++) {
          T[J][I] = Complex.zero;
        }
      } else {
        // general case

        if (I < K) {
          if (lsame(STOREV, 'C')) {
            // Skip any leading zeros.
            for (LASTV = 1; LASTV <= I - 1; LASTV++) {
              if (V[LASTV][I] != Complex.zero) break;
            }
            for (J = I + 1; J <= K; J++) {
              T[J][I] = -TAU[I] * V[N - K + I][J].conjugate();
            }
            J = max(LASTV, PREVLASTV);

            // T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)

            zgemv(
                'Conjugate transpose',
                N - K + I - J,
                K - I,
                -TAU[I],
                V(J, I + 1),
                LDV,
                V(J, I).asArray(),
                1,
                Complex.one,
                T(I + 1, I).asArray(),
                1);
          } else {
            // Skip any leading zeros.
            for (LASTV = 1; LASTV <= I - 1; LASTV++) {
              if (V[I][LASTV] != Complex.zero) break;
            }
            for (J = I + 1; J <= K; J++) {
              T[J][I] = -TAU[I] * V[J][N - K + I];
            }
            J = max(LASTV, PREVLASTV);

            // T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H

            zgemm('N', 'C', K - I, 1, N - K + I - J, -TAU[I], V(I + 1, J), LDV,
                V(I, J), LDV, Complex.one, T(I + 1, I), LDT);
          }

          // T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)

          ztrmv('Lower', 'No transpose', 'Non-unit', K - I, T(I + 1, I + 1),
              LDT, T(I + 1, I).asArray(), 1);
          if (I > 1) {
            PREVLASTV = min(PREVLASTV, LASTV);
          } else {
            PREVLASTV = LASTV;
          }
        }
        T[I][I] = TAU[I];
      }
    }
  }
}
