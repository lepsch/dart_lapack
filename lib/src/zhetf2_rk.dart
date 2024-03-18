import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zher.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zhetf2_rk(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> E_,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final IPIV = IPIV_.having();
  const ZERO = 0.0, ONE = 1.0;
  const EIGHT = 8.0, SEVTEN = 17.0;
  bool DONE, UPPER;
  int I, II, IMAX = 0, ITEMP, J, JMAX = 0, K, KK, KP = 0, KSTEP, P;
  double ABSAKK, ALPHA, COLMAX, D, D11, D22, R1, DTEMP, ROWMAX, TT, SFMIN;
  Complex D12, D21, T, WK, WKM1, WKP1;

  double CABS1(Complex Z) => Z.toDouble().abs() + Z.imaginary.abs();

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
    xerbla('ZHETF2_RK', -INFO.value);
    return;
  }

  // Initialize ALPHA for use in choosing pivot block size.

  ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

  // Compute machine safe minimum

  SFMIN = dlamch('S');

  if (UPPER) {
    // Factorize A as U*D*U**H using the upper triangle of A

    // Initialize the first entry of array E, where superdiagonal
    // elements of D are stored

    E[1] = Complex.zero;

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2

    K = N;
    if (K >= 1) {
      KSTEP = 1;
      P = K;

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

      if ((max(ABSAKK, COLMAX) == ZERO)) {
        // Column K is zero or underflow: set INFO.value and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
        A[K][K] = A[K][K]..real.toComplex();

        // Set E[K] to zero

        if (K > 1) E[K] = Complex.zero;
      } else {
        // ============================================================

        // BEGIN pivot search

        // Case(1)
        // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
        // (used to handle NaN and Inf)

        if (!(ABSAKK < ALPHA * COLMAX)) {
          // no interchange, use 1-by-1 pivot block

          KP = K;
        } else {
          DONE = false;

          // Loop until pivot found

          do {
            // BEGIN pivot search loop body

            // JMAX is the column-index of the largest off-diagonal
            // element in row IMAX, and ROWMAX is its absolute value.
            // Determine both ROWMAX and JMAX.

            if (IMAX != K) {
              JMAX = IMAX + izamax(K - IMAX, A(IMAX, IMAX + 1).asArray(), LDA);
              ROWMAX = CABS1(A[IMAX][JMAX]);
            } else {
              ROWMAX = ZERO;
            }

            if (IMAX > 1) {
              ITEMP = izamax(IMAX - 1, A(1, IMAX).asArray(), 1);
              DTEMP = CABS1(A[ITEMP][IMAX]);
              if (DTEMP > ROWMAX) {
                ROWMAX = DTEMP;
                JMAX = ITEMP;
              }
            }

            // Case(2)
            // Equivalent to testing for
            // ABS( (W( IMAX,KW-1 )).toDouble() ) >= ALPHA*ROWMAX
            // (used to handle NaN and Inf)

            if (!(A[IMAX][IMAX].toDouble().abs() < ALPHA * ROWMAX)) {
              // interchange rows and columns K and IMAX,
              // use 1-by-1 pivot block

              KP = IMAX;
              DONE = true;

              // Case(3)
              // Equivalent to testing for ROWMAX == COLMAX,
              // (used to handle NaN and Inf)
            } else if ((P == JMAX) || (ROWMAX <= COLMAX)) {
              // interchange rows and columns K-1 and IMAX,
              // use 2-by-2 pivot block

              KP = IMAX;
              KSTEP = 2;
              DONE = true;

              // Case(4)
            } else {
              // Pivot not found: set params and repeat

              P = IMAX;
              COLMAX = ROWMAX;
              IMAX = JMAX;
            }

            // END pivot search loop body
          } while (!DONE);
        }

        // END pivot search

        // ============================================================

        // KK is the column of A where pivoting step stopped

        KK = K - KSTEP + 1;

        // For only a 2x2 pivot, interchange rows and columns K and P
        // in the leading submatrix A(1:k,1:k)

        if ((KSTEP == 2) && (P != K)) {
          // (1) Swap columnar parts
          if (P > 1) zswap(P - 1, A(1, K).asArray(), 1, A(1, P).asArray(), 1);
          // (2) Swap and conjugate middle parts
          for (J = P + 1; J <= K - 1; J++) {
            T = A[J][K].conjugate();
            A[J][K] = A[P][J].conjugate();
            A[P][J] = T;
          }
          // (3) Swap and conjugate corner elements at row-col intersection
          A[P][K] = A[P][K].conjugate();
          // (4) Swap diagonal elements at row-col intersection
          R1 = A[K][K].toDouble();
          A[K][K] = A[P][P].real.toComplex();
          A[P][P] = R1.toComplex();

          // Convert upper triangle of A into U form by applying
          // the interchanges in columns k+1:N.

          if (K < N) {
            zswap(
                N - K, A(K, K + 1).asArray(), LDA, A(P, K + 1).asArray(), LDA);
          }
        }

        // For both 1x1 and 2x2 pivots, interchange rows and
        // columns KK and KP in the leading submatrix A(1:k,1:k)

        if (KP != KK) {
          // (1) Swap columnar parts
          if (KP > 1) {
            zswap(KP - 1, A(1, KK).asArray(), 1, A(1, KP).asArray(), 1);
          }
          // (2) Swap and conjugate middle parts
          for (J = KP + 1; J <= KK - 1; J++) {
            T = A[J][KK].conjugate();
            A[J][KK] = A[KP][J].conjugate();
            A[KP][J] = T;
          }
          // (3) Swap and conjugate corner elements at row-col intersection
          A[KP][KK] = A[KP][KK].conjugate();
          // (4) Swap diagonal elements at row-col intersection
          R1 = A[KK][KK].toDouble();
          A[KK][KK] = A[KP][KP].real.toComplex();
          A[KP][KP] = R1.toComplex();

          if (KSTEP == 2) {
            // (*) Make sure that diagonal element of pivot is real
            A[K][K] = A[K][K].real.toComplex();
            // (5) Swap row elements
            T = A[K - 1][K];
            A[K - 1][K] = A[KP][K];
            A[KP][K] = T;
          }

          // Convert upper triangle of A into U form by applying
          // the interchanges in columns k+1:N.

          if (K < N) {
            zswap(N - K, A(KK, K + 1).asArray(), LDA, A(KP, K + 1).asArray(),
                LDA);
          }
        } else {
          // (*) Make sure that diagonal element of pivot is real
          A[K][K] = A[K][K].real.toComplex();
          if (KSTEP == 2) A[K - 1][K - 1] = A[K - 1][K - 1].real.toComplex();
        }

        // Update the leading submatrix

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column k now holds

          // W(k) = U(k)*D(k)

          // where U(k) is the k-th column of U

          if (K > 1) {
            // Perform a rank-1 update of A(1:k-1,1:k-1) and
            // store U(k) in column k

            if (A[K][K].toDouble().abs() >= SFMIN) {
              // Perform a rank-1 update of A(1:k-1,1:k-1) as
              // A := A - U(k)*D(k)*U(k)**T
              //    = A - W(k)*1/D(k)*W(k)**T

              D11 = ONE / A[K][K].toDouble();
              zher(UPLO, K - 1, -D11, A(1, K).asArray(), 1, A, LDA);

              // Store U(k) in column k

              zdscal(K - 1, D11, A(1, K).asArray(), 1);
            } else {
              // Store L(k) in column K

              D11 = A[K][K].toDouble();
              for (II = 1; II <= K - 1; II++) {
                A[II][K] /= D11.toComplex();
              }

              // Perform a rank-1 update of A(k+1:n,k+1:n) as
              // A := A - U(k)*D(k)*U(k)**T
              //    = A - W(k)*(1/D(k))*W(k)**T
              //    = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

              zher(UPLO, K - 1, -D11, A(1, K).asArray(), 1, A, LDA);
            }

            // Store the superdiagonal element of D in array E

            E[K] = Complex.zero;
          }
        } else {
          // 2-by-2 pivot block D(k): columns k and k-1 now hold

          // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

          // where U(k) and U(k-1) are the k-th and (k-1)-th columns
          // of U

          // Perform a rank-2 update of A(1:k-2,1:k-2) as

          // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
          //    = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T

          // and store L(k) and L(k+1) in columns k and k+1

          if (K > 2) {
            // D = |A12|
            D = dlapy2(A[K - 1][K].toDouble(), A[K - 1][K].imaginary);
            D11 = (A[K][K] / D.toComplex()).toDouble();
            D22 = (A[K - 1][K - 1] / D.toComplex()).toDouble();
            D12 = A[K - 1][K] / D.toComplex();
            TT = ONE / (D11 * D22 - ONE);

            for (J = K - 2; J >= 1; J--) {
              // Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J

              WKM1 = TT.toComplex() *
                  (D11.toComplex() * A[J][K - 1] - D12.conjugate() * A[J][K]);
              WK = TT.toComplex() *
                  (D22.toComplex() * A[J][K] - D12 * A[J][K - 1]);

              // Perform a rank-2 update of A(1:k-2,1:k-2)

              for (I = J; I >= 1; I--) {
                A[I][J] -= (A[I][K] / D.toComplex()) * WK.conjugate() +
                    (A[I][K - 1] / D.toComplex()) * WKM1.conjugate();
              }

              // Store U(k) and U(k-1) in cols k and k-1 for row J

              A[J][K] = WK / D.toComplex();
              A[J][K - 1] = WKM1 / D.toComplex();
              // (*) Make sure that diagonal element of pivot is real
              A[J][J] = A[J][J].real.toComplex();
            }
          }

          // Copy superdiagonal elements of D(K) to E(K) and
          // ZERO out superdiagonal entry of A

          E[K] = A[K - 1][K];
          E[K - 1] = Complex.zero;
          A[K - 1][K] = Complex.zero;
        }

        // End column K is nonsingular
      }

      // Store details of the interchanges in IPIV

      if (KSTEP == 1) {
        IPIV[K] = KP;
      } else {
        IPIV[K] = -P;
        IPIV[K - 1] = -KP;
      }

      // Decrease K and return to the start of the main loop

      K -= KSTEP;
    }
  } else {
    // Factorize A as L*D*L**H using the lower triangle of A

    // Initialize the unused last entry of the subdiagonal array E.

    E[N] = Complex.zero;

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2

    K = 1;
    if (K <= N) {
      KSTEP = 1;
      P = K;

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

      if (max(ABSAKK, COLMAX) == ZERO) {
        // Column K is zero or underflow: set INFO.value and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
        A[K][K] = A[K][K].real.toComplex();

        // Set E[K] to zero

        if (K < N) E[K] = Complex.zero;
      } else {
        // ============================================================

        // BEGIN pivot search

        // Case(1)
        // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
        // (used to handle NaN and Inf)

        if (!(ABSAKK < ALPHA * COLMAX)) {
          // no interchange, use 1-by-1 pivot block

          KP = K;
        } else {
          DONE = false;

          // Loop until pivot found

          do {
            // BEGIN pivot search loop body

            // JMAX is the column-index of the largest off-diagonal
            // element in row IMAX, and ROWMAX is its absolute value.
            // Determine both ROWMAX and JMAX.

            if (IMAX != K) {
              JMAX = K - 1 + izamax(IMAX - K, A(IMAX, K).asArray(), LDA);
              ROWMAX = CABS1(A[IMAX][JMAX]);
            } else {
              ROWMAX = ZERO;
            }

            if (IMAX < N) {
              ITEMP = IMAX + izamax(N - IMAX, A(IMAX + 1, IMAX).asArray(), 1);
              DTEMP = CABS1(A[ITEMP][IMAX]);
              if (DTEMP > ROWMAX) {
                ROWMAX = DTEMP;
                JMAX = ITEMP;
              }
            }

            // Case(2)
            // Equivalent to testing for
            // ABS( (W( IMAX,KW-1 )).toDouble() ) >= ALPHA*ROWMAX
            // (used to handle NaN and Inf)

            if (!(A[IMAX][IMAX].toDouble().abs() < ALPHA * ROWMAX)) {
              // interchange rows and columns K and IMAX,
              // use 1-by-1 pivot block

              KP = IMAX;
              DONE = true;

              // Case(3)
              // Equivalent to testing for ROWMAX == COLMAX,
              // (used to handle NaN and Inf)
            } else if ((P == JMAX) || (ROWMAX <= COLMAX)) {
              // interchange rows and columns K+1 and IMAX,
              // use 2-by-2 pivot block

              KP = IMAX;
              KSTEP = 2;
              DONE = true;

              // Case(4)
            } else {
              // Pivot not found: set params and repeat

              P = IMAX;
              COLMAX = ROWMAX;
              IMAX = JMAX;
            }

            // END pivot search loop body
          } while (!DONE);
        }

        // END pivot search

        // ============================================================

        // KK is the column of A where pivoting step stopped

        KK = K + KSTEP - 1;

        // For only a 2x2 pivot, interchange rows and columns K and P
        // in the trailing submatrix A(k:n,k:n)

        if ((KSTEP == 2) && (P != K)) {
          // (1) Swap columnar parts
          if (P < N) {
            zswap(N - P, A(P + 1, K).asArray(), 1, A(P + 1, P).asArray(), 1);
          }
          // (2) Swap and conjugate middle parts
          for (J = K + 1; J <= P - 1; J++) {
            T = A[J][K].conjugate();
            A[J][K] = A[P][J].conjugate();
            A[P][J] = T;
          }
          // (3) Swap and conjugate corner elements at row-col intersection
          A[P][K] = A[P][K].conjugate();
          // (4) Swap diagonal elements at row-col intersection
          R1 = A[K][K].toDouble();
          A[K][K] = A[P][P].real.toComplex();
          A[P][P] = R1.toComplex();

          // Convert lower triangle of A into L form by applying
          // the interchanges in columns 1:k-1.

          if (K > 1) {
            zswap(K - 1, A(K, 1).asArray(), LDA, A(P, 1).asArray(), LDA);
          }
        }

        // For both 1x1 and 2x2 pivots, interchange rows and
        // columns KK and KP in the trailing submatrix A(k:n,k:n)

        if (KP != KK) {
          // (1) Swap columnar parts
          if (KP < N) {
            zswap(
                N - KP, A(KP + 1, KK).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }
          // (2) Swap and conjugate middle parts
          for (J = KK + 1; J <= KP - 1; J++) {
            T = A[J][KK].conjugate();
            A[J][KK] = A[KP][J].conjugate();
            A[KP][J] = T;
          }
          // (3) Swap and conjugate corner elements at row-col intersection
          A[KP][KK] = A[KP][KK].conjugate();
          // (4) Swap diagonal elements at row-col intersection
          R1 = A[KK][KK].toDouble();
          A[KK][KK] = A[KP][KP].real.toComplex();
          A[KP][KP] = R1.toComplex();

          if (KSTEP == 2) {
            // (*) Make sure that diagonal element of pivot is real
            A[K][K] = A[K][K].real.toComplex();
            // (5) Swap row elements
            T = A[K + 1][K];
            A[K + 1][K] = A[KP][K];
            A[KP][K] = T;
          }

          // Convert lower triangle of A into L form by applying
          // the interchanges in columns 1:k-1.

          if (K > 1) {
            zswap(K - 1, A(KK, 1).asArray(), LDA, A(KP, 1).asArray(), LDA);
          }
        } else {
          // (*) Make sure that diagonal element of pivot is real
          A[K][K] = A[K][K].real.toComplex();
          if (KSTEP == 2) A[K + 1][K + 1] = A[K + 1][K + 1].real.toComplex();
        }

        // Update the trailing submatrix

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column k of A now holds

          // W(k) = L(k)*D(k),

          // where L(k) is the k-th column of L

          if (K < N) {
            // Perform a rank-1 update of A(k+1:n,k+1:n) and
            // store L(k) in column k

            // Handle division by a small number

            if (A[K][K].toDouble().abs() >= SFMIN) {
              // Perform a rank-1 update of A(k+1:n,k+1:n) as
              // A := A - L(k)*D(k)*L(k)**T
              //    = A - W(k)*(1/D(k))*W(k)**T

              D11 = ONE / A[K][K].toDouble();
              zher(UPLO, N - K, -D11, A(K + 1, K).asArray(), 1, A(K + 1, K + 1),
                  LDA);

              // Store L(k) in column k

              zdscal(N - K, D11, A(K + 1, K).asArray(), 1);
            } else {
              // Store L(k) in column k

              D11 = A[K][K].toDouble();
              for (II = K + 1; II <= N; II++) {
                A[II][K] /= D11.toComplex();
              }

              // Perform a rank-1 update of A(k+1:n,k+1:n) as
              // A := A - L(k)*D(k)*L(k)**T
              //    = A - W(k)*(1/D(k))*W(k)**T
              //    = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

              zher(UPLO, N - K, -D11, A(K + 1, K).asArray(), 1, A(K + 1, K + 1),
                  LDA);
            }

            // Store the subdiagonal element of D in array E

            E[K] = Complex.zero;
          }
        } else {
          // 2-by-2 pivot block D(k): columns k and k+1 now hold

          // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

          // where L(k) and L(k+1) are the k-th and (k+1)-th columns
          // of L

          // Perform a rank-2 update of A(k+2:n,k+2:n) as

          // A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T
          //    = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T

          // and store L(k) and L(k+1) in columns k and k+1

          if (K < N - 1) {
            // D = |A21|
            D = dlapy2(A[K + 1][K].toDouble(), A[K + 1][K].imaginary);
            D11 = A[K + 1][K + 1].toDouble() / D;
            D22 = A[K][K].toDouble() / D;
            D21 = A[K + 1][K] / D.toComplex();
            TT = ONE / (D11 * D22 - ONE);

            for (J = K + 2; J <= N; J++) {
              // Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J

              WK = TT.toComplex() *
                  (D11.toComplex() * A[J][K] - D21 * A[J][K + 1]);
              WKP1 = TT.toComplex() *
                  (D22.toComplex() * A[J][K + 1] - D21.conjugate() * A[J][K]);

              // Perform a rank-2 update of A(k+2:n,k+2:n)

              for (I = J; I <= N; I++) {
                A[I][J] -= (A[I][K] / D.toComplex()) * WK.conjugate() +
                    (A[I][K + 1] / D.toComplex()) * WKP1.conjugate();
              }

              // Store L(k) and L(k+1) in cols k and k+1 for row J

              A[J][K] = WK / D.toComplex();
              A[J][K + 1] = WKP1 / D.toComplex();
              // (*) Make sure that diagonal element of pivot is real
              A[J][J] = A[J][J].real.toComplex();
            }
          }

          // Copy subdiagonal elements of D(K) to E(K) and
          // ZERO out subdiagonal entry of A

          E[K] = A[K + 1][K];
          E[K + 1] = Complex.zero;
          A[K + 1][K] = Complex.zero;
        }

        // End column K is nonsingular
      }

      // Store details of the interchanges in IPIV

      if (KSTEP == 1) {
        IPIV[K] = KP;
      } else {
        IPIV[K] = -P;
        IPIV[K + 1] = -KP;
      }

      // Increase K and return to the start of the main loop

      K += KSTEP;
    }
  }
}
