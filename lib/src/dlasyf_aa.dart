import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';

void dlasyf_aa(
  final String UPLO,
  final int J1,
  final int M,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<double> H_,
  final int LDH,
  final Array<double> WORK_,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final H = H_.having(ld: LDH);
  final WORK = WORK_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int J, K, K1, I1, I2, MJ;
  double PIV, ALPHA;

  J = 1;

  // K1 is the first column of the panel to be factorized
  // i.e.,  K1 is 2 for the first block column, and 1 for the rest of the blocks

  K1 = (2 - J1) + 1;

  if (lsame(UPLO, 'U')) {
    // .....................................................
    // Factorize A as U**T*D*U using the upper triangle of A
    // .....................................................

    while (J <= min(M, NB)) {
      // K is the column to be factorized
      //  when being called from DSYTRF_AA,
      //  > for the first block column, J1 is 1, hence J1+J-1 is J,
      //  > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,

      K = J1 + J - 1;
      if (J == M) {
        // Only need to compute T(J, J)

        MJ = 1;
      } else {
        MJ = M - J + 1;
      }

      // H(J:M, J) := A(J, J:M) - H(J:M, 1:(J-1)) * L(J1:(J-1), J),
      //  where H(J:M, J) has been initialized to be A(J, J:M)

      if (K > 2) {
        // K is the column to be factorized
        //  > for the first block column, K is J, skipping the first two
        //    columns
        //  > for the rest of the columns, K is J+1, skipping only the
        //    first column

        dgemv('No transpose', MJ, J - K1, -ONE, H(J, K1), LDH,
            A(1, J).asArray(), 1, ONE, H(J, J).asArray(), 1);
      }

      // Copy H(i:M, i) into WORK

      dcopy(MJ, H(J, J).asArray(), 1, WORK(1), 1);

      if (J > K1) {
        // Compute WORK := WORK - L(J-1, J:M) * T(J-1,J),
        //  where A(J-1, J) stores T(J-1, J) and A(J-2, J:M) stores U(J-1, J:M)

        ALPHA = -A[K - 1][J];
        daxpy(MJ, ALPHA, A(K - 2, J).asArray(), LDA, WORK(1), 1);
      }

      // Set A(J, J) = T(J, J)

      A[K][J] = WORK[1];

      if (J < M) {
        // Compute WORK(2:M) = T(J, J) L(J, (J+1):M)
        //  where A(J, J) stores T(J, J) and A(J-1, (J+1):M) stores U(J, (J+1):M)

        if (K > 1) {
          ALPHA = -A[K][J];
          daxpy(M - J, ALPHA, A(K - 1, J + 1).asArray(), LDA, WORK(2), 1);
        }

        // Find max(|WORK(2:M)|)

        I2 = idamax(M - J, WORK(2), 1) + 1;
        PIV = WORK[I2];

        // Apply symmetric pivot

        if ((I2 != 2) && (PIV != 0)) {
          // Swap WORK(I1) and WORK(I2)

          I1 = 2;
          WORK[I2] = WORK[I1];
          WORK[I1] = PIV;

          // Swap A(I1, I1+1:M) with A(I1+1:M, I2)

          I1 += J - 1;
          I2 += J - 1;
          dswap(I2 - I1 - 1, A(J1 + I1 - 1, I1 + 1).asArray(), LDA,
              A(J1 + I1, I2).asArray(), 1);

          // Swap A(I1, I2+1:M) with A(I2, I2+1:M)

          if (I2 < M) {
            dswap(M - I2, A(J1 + I1 - 1, I2 + 1).asArray(), LDA,
                A(J1 + I2 - 1, I2 + 1).asArray(), LDA);
          }

          // Swap A(I1, I1) with A(I2,I2)

          PIV = A[I1 + J1 - 1][I1];
          A[J1 + I1 - 1][I1] = A[J1 + I2 - 1][I2];
          A[J1 + I2 - 1][I2] = PIV;

          // Swap H(I1, 1:J1) with H(I2, 1:J1)

          dswap(I1 - 1, H(I1, 1).asArray(), LDH, H(I2, 1).asArray(), LDH);
          IPIV[I1] = I2;

          if (I1 > (K1 - 1)) {
            // Swap L(1:I1-1, I1) with L(1:I1-1, I2),
            //  skipping the first column

            dswap(I1 - K1 + 1, A(1, I1).asArray(), 1, A(1, I2).asArray(), 1);
          }
        } else {
          IPIV[J + 1] = J + 1;
        }

        // Set A(J, J+1) = T(J, J+1)

        A[K][J + 1] = WORK[2];

        if (J < NB) {
          // Copy A(J+1:M, J+1) into H(J:M, J),

          dcopy(M - J, A(K + 1, J + 1).asArray(), LDA,
              H(J + 1, J + 1).asArray(), 1);
        }

        // Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1),
        //  where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1)

        if (J < (M - 1)) {
          if (A[K][J + 1] != ZERO) {
            ALPHA = ONE / A[K][J + 1];
            dcopy(M - J - 1, WORK(3), 1, A(K, J + 2).asArray(), LDA);
            dscal(M - J - 1, ALPHA, A(K, J + 2).asArray(), LDA);
          } else {
            dlaset('Full', 1, M - J - 1, ZERO, ZERO, A(K, J + 2), LDA);
          }
        }
      }
      J++;
    }
  } else {
    // .....................................................
    // Factorize A as L*D*L**T using the lower triangle of A
    // .....................................................

    while (J <= min(M, NB)) {
      // K is the column to be factorized
      //  when being called from DSYTRF_AA,
      //  > for the first block column, J1 is 1, hence J1+J-1 is J,
      //  > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,

      K = J1 + J - 1;
      if (J == M) {
        // Only need to compute T(J, J)

        MJ = 1;
      } else {
        MJ = M - J + 1;
      }

      // H(J:M, J) := A(J:M, J) - H(J:M, 1:(J-1)) * L(J, J1:(J-1))^T,
      //  where H(J:M, J) has been initialized to be A(J:M, J)

      if (K > 2) {
        // K is the column to be factorized
        //  > for the first block column, K is J, skipping the first two
        //    columns
        //  > for the rest of the columns, K is J+1, skipping only the
        //    first column

        dgemv('No transpose', MJ, J - K1, -ONE, H(J, K1), LDH,
            A(J, 1).asArray(), LDA, ONE, H(J, J).asArray(), 1);
      }

      // Copy H(J:M, J) into WORK

      dcopy(MJ, H(J, J).asArray(), 1, WORK(1), 1);

      if (J > K1) {
        // Compute WORK := WORK - L(J:M, J-1) * T(J-1,J),
        //  where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1)

        ALPHA = -A[J][K - 1];
        daxpy(MJ, ALPHA, A(J, K - 2).asArray(), 1, WORK(1), 1);
      }

      // Set A(J, J) = T(J, J)

      A[J][K] = WORK[1];

      if (J < M) {
        // Compute WORK(2:M) = T(J, J) L((J+1):M, J)
        //  where A(J, J) = T(J, J) and A((J+1):M, J-1) = L((J+1):M, J)

        if (K > 1) {
          ALPHA = -A[J][K];
          daxpy(M - J, ALPHA, A(J + 1, K - 1).asArray(), 1, WORK(2), 1);
        }

        // Find max(|WORK(2:M)|)

        I2 = idamax(M - J, WORK(2), 1) + 1;
        PIV = WORK[I2];

        // Apply symmetric pivot

        if ((I2 != 2) && (PIV != 0)) {
          // Swap WORK(I1) and WORK(I2)

          I1 = 2;
          WORK[I2] = WORK[I1];
          WORK[I1] = PIV;

          // Swap A(I1+1:M, I1) with A(I2, I1+1:M)

          I1 += J - 1;
          I2 += J - 1;
          dswap(I2 - I1 - 1, A(I1 + 1, J1 + I1 - 1).asArray(), 1,
              A(I2, J1 + I1).asArray(), LDA);

          // Swap A(I2+1:M, I1) with A(I2+1:M, I2)

          if (I2 < M) {
            dswap(M - I2, A(I2 + 1, J1 + I1 - 1).asArray(), 1,
                A(I2 + 1, J1 + I2 - 1).asArray(), 1);
          }

          // Swap A(I1, I1) with A(I2, I2)

          PIV = A[I1][J1 + I1 - 1];
          A[I1][J1 + I1 - 1] = A[I2][J1 + I2 - 1];
          A[I2][J1 + I2 - 1] = PIV;

          // Swap H(I1, I1:J1) with H(I2, I2:J1)

          dswap(I1 - 1, H(I1, 1).asArray(), LDH, H(I2, 1).asArray(), LDH);
          IPIV[I1] = I2;

          if (I1 > (K1 - 1)) {
            // Swap L(1:I1-1, I1) with L(1:I1-1, I2),
            //  skipping the first column

            dswap(
                I1 - K1 + 1, A(I1, 1).asArray(), LDA, A(I2, 1).asArray(), LDA);
          }
        } else {
          IPIV[J + 1] = J + 1;
        }

        // Set A(J+1, J) = T(J+1, J)

        A[J + 1][K] = WORK[2];

        if (J < NB) {
          // Copy A(J+1:M, J+1) into H(J+1:M, J),

          dcopy(M - J, A(J + 1, K + 1).asArray(), 1, H(J + 1, J + 1).asArray(),
              1);
        }

        // Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1),
        //  where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1)

        if (J < (M - 1)) {
          if (A[J + 1][K] != ZERO) {
            ALPHA = ONE / A[J + 1][K];
            dcopy(M - J - 1, WORK(3), 1, A(J + 2, K).asArray(), 1);
            dscal(M - J - 1, ALPHA, A(J + 2, K).asArray(), 1);
          } else {
            dlaset('Full', M - J - 1, 1, ZERO, ZERO, A(J + 2, K), LDA);
          }
        }
      }
      J++;
    }
  }
}
