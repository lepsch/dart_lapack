import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zher.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zhetf2(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();

  const ZERO = 0.0, ONE = 1.0;
  const EIGHT = 8.0, SEVTEN = 17.0;
  bool UPPER;
  int I, IMAX = 0, J, JMAX, K, KK, KP, KSTEP;
  double ABSAKK, ALPHA, COLMAX, D, D11, D22, R1, ROWMAX, TT;
  Complex D12, D21, T, WK, WKM1, WKP1;

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZHETF2', -INFO.value);
    return;
  }

  // Initialize ALPHA for use in choosing pivot block size.

  ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

  if (UPPER) {
    // Factorize A as U*D*U**H using the upper triangle of A

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2

    K = N;
    while (K >= 1) {
      KSTEP = 1;

      // Determine rows and columns to be interchanged and whether
      // a 1-by-1 or 2-by-2 pivot block will be used

      ABSAKK = A[K][K].toDouble().abs();

      // IMAX is the row-index of the largest off-diagonal element in
      // column K, and COLMAX is its absolute value.
      // Determine both COLMAX and IMAX.

      if (K > 1) {
        IMAX = izamax(K - 1, A(1, K).asArray(), 1);
        COLMAX = CABS1(A[IMAX][K]);
      } else {
        COLMAX = ZERO;
      }

      if ((max(ABSAKK, COLMAX) == ZERO) || disnan(ABSAKK)) {
        // Column K is zero or underflow, or contains a NaN:
        // set INFO.value and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
        A[K][K] = A[K][K].real.toComplex();
      } else {
        // ============================================================

        // Test for interchange

        if (ABSAKK >= ALPHA * COLMAX) {
          // no interchange, use 1-by-1 pivot block

          KP = K;
        } else {
          // JMAX is the column-index of the largest off-diagonal
          // element in row IMAX, and ROWMAX is its absolute value.
          // Determine only ROWMAX.

          JMAX = IMAX + izamax(K - IMAX, A(IMAX, IMAX + 1).asArray(), LDA);
          ROWMAX = CABS1(A[IMAX][JMAX]);
          if (IMAX > 1) {
            JMAX = izamax(IMAX - 1, A(1, IMAX).asArray(), 1);
            ROWMAX = max(ROWMAX, CABS1(A[JMAX][IMAX]));
          }

          if (ABSAKK >= ALPHA * COLMAX * (COLMAX / ROWMAX)) {
            // no interchange, use 1-by-1 pivot block

            KP = K;
          } else if (A[IMAX][IMAX].toDouble().abs() >= ALPHA * ROWMAX) {
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

        // ============================================================

        KK = K - KSTEP + 1;
        if (KP != KK) {
          // Interchange rows and columns KK and KP in the leading
          // submatrix A(1:k,1:k)

          zswap(KP - 1, A(1, KK).asArray(), 1, A(1, KP).asArray(), 1);
          for (J = KP + 1; J <= KK - 1; J++) {
            T = A[J][KK].conjugate();
            A[J][KK] = A[KP][J].conjugate();
            A[KP][J] = T;
          }
          A[KP][KK] = A[KP][KK].conjugate();
          R1 = (A[KK][KK]).toDouble();
          A[KK][KK] = (A[KP][KP]).real.toComplex();
          A[KP][KP] = R1.toComplex();
          if (KSTEP == 2) {
            A[K][K] = (A[K][K]).real.toComplex();
            T = A[K - 1][K];
            A[K - 1][K] = A[KP][K];
            A[KP][K] = T;
          }
        } else {
          A[K][K] = A[K][K].real.toComplex();
          if (KSTEP == 2) A[K - 1][K - 1] = A[K - 1][K - 1].real.toComplex();
        }

        // Update the leading submatrix

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column k now holds

          // W(k) = U(k)*D(k)

          // where U(k) is the k-th column of U

          // Perform a rank-1 update of A(1:k-1,1:k-1) as

          // A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H

          R1 = ONE / (A[K][K]).toDouble();
          zher(UPLO, K - 1, -R1, A(1, K).asArray(), 1, A, LDA);

          // Store U(k) in column k

          zdscal(K - 1, R1, A(1, K).asArray(), 1);
        } else {
          // 2-by-2 pivot block D(k): columns k and k-1 now hold

          // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

          // where U(k) and U(k-1) are the k-th and (k-1)-th columns
          // of U

          // Perform a rank-2 update of A(1:k-2,1:k-2) as

          // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H
          //    = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H

          if (K > 2) {
            D = dlapy2(A[K - 1][K].toDouble(), A[K - 1][K].imaginary);
            D22 = (A[K - 1][K - 1]).toDouble() / D;
            D11 = (A[K][K]).toDouble() / D;
            TT = ONE / (D11 * D22 - ONE);
            D12 = A[K - 1][K] / D.toComplex();
            D = TT / D;

            for (J = K - 2; J >= 1; J--) {
              WKM1 = D.toComplex() *
                  (D11.toComplex() * A[J][K - 1] - D12.conjugate() * A[J][K]);
              WK = D.toComplex() *
                  (D22.toComplex() * A[J][K] - D12 * A[J][K - 1]);
              for (I = J; I >= 1; I--) {
                A[I][J] -=
                    A[I][K] * WK.conjugate() - A[I][K - 1] * WKM1.conjugate();
              }
              A[J][K] = WK;
              A[J][K - 1] = WKM1;
              A[J][J] = A[J][J].real.toComplex();
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
    }
  } else {
    // Factorize A as L*D*L**H using the lower triangle of A

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2

    K = 1;
    while (K <= N) {
      KSTEP = 1;

      // Determine rows and columns to be interchanged and whether
      // a 1-by-1 or 2-by-2 pivot block will be used

      ABSAKK = A[K][K].toDouble().abs();

      // IMAX is the row-index of the largest off-diagonal element in
      // column K, and COLMAX is its absolute value.
      // Determine both COLMAX and IMAX.

      if (K < N) {
        IMAX = K + izamax(N - K, A(K + 1, K).asArray(), 1);
        COLMAX = CABS1(A[IMAX][K]);
      } else {
        COLMAX = ZERO;
      }

      if ((max(ABSAKK, COLMAX) == ZERO) || disnan(ABSAKK)) {
        // Column K is zero or underflow, or contains a NaN:
        // set INFO.value and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
        A[K][K] = A[K][K].real.toComplex();
      } else {
        // ============================================================

        // Test for interchange

        if (ABSAKK >= ALPHA * COLMAX) {
          // no interchange, use 1-by-1 pivot block

          KP = K;
        } else {
          // JMAX is the column-index of the largest off-diagonal
          // element in row IMAX, and ROWMAX is its absolute value.
          // Determine only ROWMAX.

          JMAX = K - 1 + izamax(IMAX - K, A(IMAX, K).asArray(), LDA);
          ROWMAX = CABS1(A[IMAX][JMAX]);
          if (IMAX < N) {
            JMAX = IMAX + izamax(N - IMAX, A(IMAX + 1, IMAX).asArray(), 1);
            ROWMAX = max(ROWMAX, CABS1(A[JMAX][IMAX]));
          }

          if (ABSAKK >= ALPHA * COLMAX * (COLMAX / ROWMAX)) {
            // no interchange, use 1-by-1 pivot block

            KP = K;
          } else if (A[IMAX][IMAX].toDouble().abs() >= ALPHA * ROWMAX) {
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

        // ============================================================

        KK = K + KSTEP - 1;
        if (KP != KK) {
          // Interchange rows and columns KK and KP in the trailing
          // submatrix A(k:n,k:n)

          if (KP < N) {
            zswap(
                N - KP, A(KP + 1, KK).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }
          for (J = KK + 1; J <= KP - 1; J++) {
            T = A[J][KK].conjugate();
            A[J][KK] = A[KP][J].conjugate();
            A[KP][J] = T;
          }
          A[KP][KK] = A[KP][KK].conjugate();
          R1 = (A[KK][KK]).toDouble();
          A[KK][KK] = (A[KP][KP]).real.toComplex();
          A[KP][KP] = R1.toComplex();
          if (KSTEP == 2) {
            A[K][K] = A[K][K].real.toComplex();
            T = A[K + 1][K];
            A[K + 1][K] = A[KP][K];
            A[KP][K] = T;
          }
        } else {
          A[K][K] = (A[K][K]).real.toComplex();
          if (KSTEP == 2) A[K + 1][K + 1] = A[K + 1][K + 1].real.toComplex();
        }

        // Update the trailing submatrix

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column k now holds

          // W(k) = L(k)*D(k)

          // where L(k) is the k-th column of L

          if (K < N) {
            // Perform a rank-1 update of A(k+1:n,k+1:n) as

            // A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H

            R1 = ONE / (A[K][K]).toDouble();
            zher(UPLO, N - K, -R1, A(K + 1, K).asArray(), 1, A(K + 1, K + 1),
                LDA);

            // Store L(k) in column K

            zdscal(N - K, R1, A(K + 1, K).asArray(), 1);
          }
        } else {
          // 2-by-2 pivot block D(k)

          if (K < N - 1) {
            // Perform a rank-2 update of A(k+2:n,k+2:n) as

            // A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H
            //    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H

            // where L(k) and L(k+1) are the k-th and (k+1)-th
            // columns of L

            D = dlapy2(A[K + 1][K].real, A[K + 1][K].imaginary);
            D11 = (A[K + 1][K + 1]).toDouble() / D;
            D22 = (A[K][K]).toDouble() / D;
            TT = ONE / (D11 * D22 - ONE);
            D21 = A[K + 1][K] / D.toComplex();
            D = TT / D;

            for (J = K + 2; J <= N; J++) {
              WK = D.toComplex() *
                  (D11.toComplex() * A[J][K] - D21 * A[J][K + 1]);
              WKP1 = D.toComplex() *
                  (D22.toComplex() * A[J][K + 1] - D21.conjugate() * A[J][K]);
              for (I = J; I <= N; I++) {
                A[I][J] -=
                    A[I][K] * WK.conjugate() - A[I][K + 1] * WKP1.conjugate();
              }
              A[J][K] = WK;
              A[J][K + 1] = WKP1;
              A[J][J] = A[J][J].real.toComplex();
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
    }
  }
}
