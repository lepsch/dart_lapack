// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zhpr.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dlapy2.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zhptrf(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final IPIV = IPIV_.having();
  const ZERO = 0.0, ONE = 1.0;
  const EIGHT = 8.0, SEVTEN = 17.0;
  bool UPPER;
  int I, IMAX = 0, J, JMAX, K, KC, KK, KNC, KP, KPC = 0, KSTEP, KX, NPP;
  double ABSAKK, ALPHA, COLMAX, D, D11, D22, R1, ROWMAX, TT;
  Complex D12, D21, T, WK, WKM1, WKP1;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('ZHPTRF', -INFO.value);
    return;
  }

  // Initialize ALPHA for use in choosing pivot block size.

  ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

  if (UPPER) {
    // Factorize A as U*D*U**H using the upper triangle of A

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2

    K = N;
    KC = (N - 1) * N ~/ 2 + 1;
    while (true) {
      KNC = KC;

      // If K < 1, exit from loop

      if (K < 1) return;
      KSTEP = 1;

      // Determine rows and columns to be interchanged and whether
      // a 1-by-1 or 2-by-2 pivot block will be used

      ABSAKK = AP[KC + K - 1].real.abs();

      // IMAX is the row-index of the largest off-diagonal element in
      // column K, and COLMAX is its absolute value

      if (K > 1) {
        IMAX = izamax(K - 1, AP(KC), 1);
        COLMAX = AP[KC + IMAX - 1].cabs1();
      } else {
        COLMAX = ZERO;
      }

      if (max(ABSAKK, COLMAX) == ZERO) {
        // Column K is zero: set INFO and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
        AP[KC + K - 1] = AP[KC + K - 1].real.toComplex();
      } else {
        if (ABSAKK >= ALPHA * COLMAX) {
          // no interchange, use 1-by-1 pivot block

          KP = K;
        } else {
          // JMAX is the column-index of the largest off-diagonal
          // element in row IMAX, and ROWMAX is its absolute value

          ROWMAX = ZERO;
          JMAX = IMAX;
          KX = IMAX * (IMAX + 1) ~/ 2 + IMAX;
          for (J = IMAX + 1; J <= K; J++) {
            if (AP[KX].cabs1() > ROWMAX) {
              ROWMAX = AP[KX].cabs1();
              JMAX = J;
            }
            KX += J;
          }
          KPC = (IMAX - 1) * IMAX ~/ 2 + 1;
          if (IMAX > 1) {
            JMAX = izamax(IMAX - 1, AP(KPC), 1);
            ROWMAX = max(ROWMAX, AP[KPC + JMAX - 1].cabs1());
          }

          if (ABSAKK >= ALPHA * COLMAX * (COLMAX / ROWMAX)) {
            // no interchange, use 1-by-1 pivot block

            KP = K;
          } else if (AP[KPC + IMAX - 1].real.abs() >= ALPHA * ROWMAX) {
            // interchange rows and columns K and IMAX, use 1-by-1
            // pivot block

            KP = IMAX;
          } else {
            // interchange rows and columns K-1 and IMAX, use 2-by-2
            // pivot block

            KP = IMAX;
            KSTEP = 2;
          }
        }

        KK = K - KSTEP + 1;
        if (KSTEP == 2) KNC -= K - 1;
        if (KP != KK) {
          // Interchange rows and columns KK and KP in the leading
          // submatrix A(1:k,1:k)

          zswap(KP - 1, AP(KNC), 1, AP(KPC), 1);
          KX = KPC + KP - 1;
          for (J = KP + 1; J <= KK - 1; J++) {
            KX += J - 1;
            T = AP[KNC + J - 1].conjugate();
            AP[KNC + J - 1] = AP[KX].conjugate();
            AP[KX] = T;
          }
          AP[KX + KK - 1] = AP[KX + KK - 1].conjugate();
          R1 = AP[KNC + KK - 1].real;
          AP[KNC + KK - 1] = AP[KPC + KP - 1].real.toComplex();
          AP[KPC + KP - 1] = R1.toComplex();
          if (KSTEP == 2) {
            AP[KC + K - 1] = AP[KC + K - 1].real.toComplex();
            T = AP[KC + K - 2];
            AP[KC + K - 2] = AP[KC + KP - 1];
            AP[KC + KP - 1] = T;
          }
        } else {
          AP[KC + K - 1] = AP[KC + K - 1].real.toComplex();
          if (KSTEP == 2) AP[KC - 1] = AP[KC - 1].real.toComplex();
        }

        // Update the leading submatrix

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column k now holds
          //
          // W(k) = U(k)*D(k)
          //
          // where U(k) is the k-th column of U
          //
          // Perform a rank-1 update of A(1:k-1,1:k-1) as
          //
          // A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H

          R1 = ONE / AP[KC + K - 1].real;
          zhpr(UPLO, K - 1, -R1, AP(KC), 1, AP);

          // Store U(k) in column k

          zdscal(K - 1, R1, AP(KC), 1);
        } else {
          // 2-by-2 pivot block D(k): columns k and k-1 now hold
          //
          // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
          //
          // where U(k) and U(k-1) are the k-th and (k-1)-th columns
          // of U
          //
          // Perform a rank-2 update of A(1:k-2,1:k-2) as
          //
          // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H
          //    = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H

          if (K > 2) {
            D = dlapy2(AP[K - 1 + (K - 1) * K ~/ 2].real,
                AP[K - 1 + (K - 1) * K ~/ 2].imaginary);
            D22 = AP[K - 1 + (K - 2) * (K - 1) ~/ 2].real / D;
            D11 = AP[K + (K - 1) * K ~/ 2].real / D;
            TT = ONE / (D11 * D22 - ONE);
            D12 = AP[K - 1 + (K - 1) * K ~/ 2] / D.toComplex();
            D = TT / D;

            for (J = K - 2; J >= 1; J--) {
              WKM1 = D.toComplex() *
                  (D11.toComplex() * AP[J + (K - 2) * (K - 1) ~/ 2] -
                      D12.conjugate() * AP[J + (K - 1) * K ~/ 2]);
              WK = D.toComplex() *
                  (D22.toComplex() * AP[J + (K - 1) * K ~/ 2] -
                      D12 * AP[J + (K - 2) * (K - 1) ~/ 2]);
              for (I = J; I >= 1; I--) {
                AP[I + (J - 1) * J ~/ 2] -=
                    AP[I + (K - 1) * K ~/ 2] * WK.conjugate() +
                        AP[I + (K - 2) * (K - 1) ~/ 2] * WKM1.conjugate();
              }
              AP[J + (K - 1) * K ~/ 2] = WK;
              AP[J + (K - 2) * (K - 1) ~/ 2] = WKM1;
              AP[J + (J - 1) * J ~/ 2] =
                  AP[J + (J - 1) * J ~/ 2].real.toComplex();
            }
          }
        }
      }

      // Store details of the interchanges in IPIV

      if (KSTEP == 1) {
        IPIV[K] = KP;
      } else {
        IPIV[K] = -KP;
        IPIV[K - 1] = -KP;
      }

      // Decrease K and return to the start of the main loop

      K -= KSTEP;
      KC = KNC - K;
    }
  } else {
    // Factorize A as L*D*L**H using the lower triangle of A

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2

    K = 1;
    KC = 1;
    NPP = N * (N + 1) ~/ 2;
    while (true) {
      KNC = KC;

      // If K > N, exit from loop

      if (K > N) return;
      KSTEP = 1;

      // Determine rows and columns to be interchanged and whether
      // a 1-by-1 or 2-by-2 pivot block will be used

      ABSAKK = AP[KC].real.abs();

      // IMAX is the row-index of the largest off-diagonal element in
      // column K, and COLMAX is its absolute value

      if (K < N) {
        IMAX = K + izamax(N - K, AP(KC + 1), 1);
        COLMAX = AP[KC + IMAX - K].cabs1();
      } else {
        COLMAX = ZERO;
      }

      if (max(ABSAKK, COLMAX) == ZERO) {
        // Column K is zero: set INFO and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
        AP[KC] = AP[KC].real.toComplex();
      } else {
        if (ABSAKK >= ALPHA * COLMAX) {
          // no interchange, use 1-by-1 pivot block

          KP = K;
        } else {
          // JMAX is the column-index of the largest off-diagonal
          // element in row IMAX, and ROWMAX is its absolute value

          ROWMAX = ZERO;
          KX = KC + IMAX - K;
          for (J = K; J <= IMAX - 1; J++) {
            if (AP[KX].cabs1() > ROWMAX) {
              ROWMAX = AP[KX].cabs1();
              JMAX = J;
            }
            KX += N - J;
          }
          KPC = NPP - (N - IMAX + 1) * (N - IMAX + 2) ~/ 2 + 1;
          if (IMAX < N) {
            JMAX = IMAX + izamax(N - IMAX, AP(KPC + 1), 1);
            ROWMAX = max(ROWMAX, AP[KPC + JMAX - IMAX].cabs1());
          }

          if (ABSAKK >= ALPHA * COLMAX * (COLMAX / ROWMAX)) {
            // no interchange, use 1-by-1 pivot block

            KP = K;
          } else if (AP[KPC].real.abs() >= ALPHA * ROWMAX) {
            // interchange rows and columns K and IMAX, use 1-by-1
            // pivot block

            KP = IMAX;
          } else {
            // interchange rows and columns K+1 and IMAX, use 2-by-2
            // pivot block

            KP = IMAX;
            KSTEP = 2;
          }
        }

        KK = K + KSTEP - 1;
        if (KSTEP == 2) KNC += N - K + 1;
        if (KP != KK) {
          // Interchange rows and columns KK and KP in the trailing
          // submatrix A(k:n,k:n)

          if (KP < N) zswap(N - KP, AP(KNC + KP - KK + 1), 1, AP(KPC + 1), 1);
          KX = KNC + KP - KK;
          for (J = KK + 1; J <= KP - 1; J++) {
            KX += N - J + 1;
            T = AP[KNC + J - KK].conjugate();
            AP[KNC + J - KK] = AP[KX].conjugate();
            AP[KX] = T;
          }
          AP[KNC + KP - KK] = AP[KNC + KP - KK].conjugate();
          R1 = AP[KNC].real;
          AP[KNC] = AP[KPC].real.toComplex();
          AP[KPC] = R1.toComplex();
          if (KSTEP == 2) {
            AP[KC] = AP[KC].real.toComplex();
            T = AP[KC + 1];
            AP[KC + 1] = AP[KC + KP - K];
            AP[KC + KP - K] = T;
          }
        } else {
          AP[KC] = AP[KC].real.toComplex();
          if (KSTEP == 2) AP[KNC] = AP[KNC].real.toComplex();
        }

        // Update the trailing submatrix

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column k now holds

          // W(k) = L(k)*D(k)

          // where L(k) is the k-th column of L

          if (K < N) {
            // Perform a rank-1 update of A(k+1:n,k+1:n) as

            // A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H

            R1 = ONE / AP[KC].real;
            zhpr(UPLO, N - K, -R1, AP(KC + 1), 1, AP(KC + N - K + 1));

            // Store L(k) in column K

            zdscal(N - K, R1, AP(KC + 1), 1);
          }
        } else {
          // 2-by-2 pivot block D(k): columns K and K+1 now hold

          // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

          // where L(k) and L(k+1) are the k-th and (k+1)-th columns
          // of L

          if (K < N - 1) {
            // Perform a rank-2 update of A(k+2:n,k+2:n) as

            // A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H
            //    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H

            // where L(k) and L(k+1) are the k-th and (k+1)-th
            // columns of L

            D = dlapy2(AP[K + 1 + (K - 1) * (2 * N - K) ~/ 2].real,
                AP[K + 1 + (K - 1) * (2 * N - K) ~/ 2].imaginary);
            D11 = AP[K + 1 + K * (2 * N - K - 1) ~/ 2].real / D;
            D22 = AP[K + (K - 1) * (2 * N - K) ~/ 2].real / D;
            TT = ONE / (D11 * D22 - ONE);
            D21 = AP[K + 1 + (K - 1) * (2 * N - K) ~/ 2] / D.toComplex();
            D = TT / D;

            for (J = K + 2; J <= N; J++) {
              WK = D.toComplex() *
                  (D11.toComplex() * AP[J + (K - 1) * (2 * N - K) ~/ 2] -
                      D21 * AP[J + K * (2 * N - K - 1) ~/ 2]);
              WKP1 = D.toComplex() *
                  (D22.toComplex() * AP[J + K * (2 * N - K - 1) ~/ 2] -
                      D21.conjugate() * AP[J + (K - 1) * (2 * N - K) ~/ 2]);
              for (I = J; I <= N; I++) {
                AP[I + (J - 1) * (2 * N - J) ~/ 2] =
                    AP[I + (J - 1) * (2 * N - J) ~/ 2] -
                        AP[I + (K - 1) * (2 * N - K) ~/ 2] * WK.conjugate() -
                        AP[I + K * (2 * N - K - 1) ~/ 2] * WKP1.conjugate();
              }
              AP[J + (K - 1) * (2 * N - K) ~/ 2] = WK;
              AP[J + K * (2 * N - K - 1) ~/ 2] = WKP1;
              AP[J + (J - 1) * (2 * N - J) ~/ 2] =
                  AP[J + (J - 1) * (2 * N - J) ~/ 2].real.toComplex();
            }
          }
        }
      }

      // Store details of the interchanges in IPIV

      if (KSTEP == 1) {
        IPIV[K] = KP;
      } else {
        IPIV[K] = -KP;
        IPIV[K + 1] = -KP;
      }

      // Increase K and return to the start of the main loop

      K += KSTEP;
      KC = KNC + N - K + 2;
    }
  }
}
