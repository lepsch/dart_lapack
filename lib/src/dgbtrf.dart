import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbtf2.dart';
import 'package:lapack/src/dlaswp.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgbtrf(
  final int M,
  final int N,
  final int KL,
  final int KU,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final IPIV = IPIV_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NBMAX = 64, LDWORK = NBMAX + 1;
  int I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, JU, K2, KM, KV, NB, NW;
  double TEMP;
  final WORK13 = Matrix<double>(LDWORK, NBMAX),
      WORK31 = Matrix<double>(LDWORK, NBMAX);

  // KV is the number of superdiagonals in the factor U, allowing for
  // fill-in

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
    xerbla('DGBTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Determine the block size for this environment

  NB = ilaenv(1, 'DGBTRF', ' ', M, N, KL, KU);

  // The block size must not exceed the limit set by the size of the
  // local arrays WORK13 and WORK31.

  NB = min(NB, NBMAX);

  if (NB <= 1 || NB > KL) {
    // Use unblocked code

    dgbtf2(M, N, KL, KU, AB, LDAB, IPIV, INFO);
  } else {
    // Use blocked code

    // Zero the superdiagonal elements of the work array WORK13

    for (J = 1; J <= NB; J++) {
      for (I = 1; I <= J - 1; I++) {
        WORK13[I][J] = ZERO;
      }
    }

    // Zero the subdiagonal elements of the work array WORK31

    for (J = 1; J <= NB; J++) {
      for (I = J + 1; I <= NB; I++) {
        WORK31[I][J] = ZERO;
      }
    }

    // Gaussian elimination with partial pivoting

    // Set fill-in elements in columns KU+2 to KV to zero

    for (J = KU + 2; J <= min(KV, N); J++) {
      for (I = KV - J + 2; I <= KL; I++) {
        AB[I][J] = ZERO;
      }
    }

    // JU is the index of the last column affected by the current
    // stage of the factorization

    JU = 1;

    for (J = 1; NB < 0 ? J >= min(M, N) : J <= min(M, N); J += NB) {
      JB = min(NB, min(M, N) - J + 1);

      // The active part of the matrix is partitioned
      //
      //    A11   A12   A13
      //    A21   A22   A23
      //    A31   A32   A33
      //
      // Here A11, A21 and A31 denote the current block of JB columns
      // which is about to be factorized. The number of rows in the
      // partitioning are JB, I2, I3 respectively, and the numbers
      // of columns are JB, J2, J3. The superdiagonal elements of A13
      // and the subdiagonal elements of A31 lie outside the band.

      I2 = min(KL - JB, M - J - JB + 1);
      I3 = min(JB, M - J - KL + 1);

      // J2 and J3 are computed after JU has been updated.

      // Factorize the current block of JB columns

      for (JJ = J; JJ <= J + JB - 1; JJ++) {
        // Set fill-in elements in column JJ+KV to zero

        if (JJ + KV <= N) {
          for (I = 1; I <= KL; I++) {
            AB[I][JJ + KV] = ZERO;
          }
        }

        // Find pivot and test for singularity. KM is the number of
        // subdiagonal elements in the current column.

        KM = min(KL, M - JJ);
        JP = idamax(KM + 1, AB(KV + 1, JJ).asArray(), 1);
        IPIV[JJ] = JP + JJ - J;
        if (AB[KV + JP][JJ] != ZERO) {
          JU = max(JU, min(JJ + KU + JP - 1, N));
          if (JP != 1) {
            // Apply interchange to columns J to J+JB-1

            if (JP + JJ - 1 < J + KL) {
              dswap(JB, AB(KV + 1 + JJ - J, J).asArray(), LDAB - 1,
                  AB(KV + JP + JJ - J, J).asArray(), LDAB - 1);
            } else {
              // The interchange affects columns J to JJ-1 of A31
              // which are stored in the work array WORK31

              dswap(JJ - J, AB(KV + 1 + JJ - J, J).asArray(), LDAB - 1,
                  WORK31(JP + JJ - J - KL, 1).asArray(), LDWORK);
              dswap(J + JB - JJ, AB(KV + 1, JJ).asArray(), LDAB - 1,
                  AB(KV + JP, JJ).asArray(), LDAB - 1);
            }
          }

          // Compute multipliers

          dscal(KM, ONE / AB[KV + 1][JJ], AB(KV + 2, JJ).asArray(), 1);

          // Update trailing submatrix within the band and within
          // the current block. JM is the index of the last column
          // which needs to be updated.

          JM = min(JU, J + JB - 1);
          if (JM > JJ) {
            dger(
                KM,
                JM - JJ,
                -ONE,
                AB(KV + 2, JJ).asArray(),
                1,
                AB(KV, JJ + 1).asArray(),
                LDAB - 1,
                AB(KV + 1, JJ + 1),
                LDAB - 1);
          }
        } else {
          // If pivot is zero, set INFO.value to the index of the pivot
          // unless a zero pivot has already been found.

          if (INFO.value == 0) INFO.value = JJ;
        }

        // Copy current column of A31 into the work array WORK31

        NW = min(JJ - J + 1, I3);
        if (NW > 0) {
          dcopy(NW, AB(KV + KL + 1 - JJ + J, JJ).asArray(), 1,
              WORK31(1, JJ - J + 1).asArray(), 1);
        }
      }
      if (J + JB <= N) {
        // Apply the row interchanges to the other blocks.

        J2 = min(JU - J + 1, KV) - JB;
        J3 = max(0, JU - J - KV + 1);

        // Use DLASWP to apply the row interchanges to A12, A22, and
        // A32.

        dlaswp(J2, AB(KV + 1 - JB, J + JB), LDAB - 1, 1, JB, IPIV(J), 1);

        // Adjust the pivot indices.

        for (I = J; I <= J + JB - 1; I++) {
          IPIV[I] += J - 1;
        }

        // Apply the row interchanges to A13, A23, and A33
        // columnwise.

        K2 = J - 1 + JB + J2;
        for (I = 1; I <= J3; I++) {
          JJ = K2 + I;
          for (II = J + I - 1; II <= J + JB - 1; II++) {
            IP = IPIV[II];
            if (IP != II) {
              TEMP = AB[KV + 1 + II - JJ][JJ];
              AB[KV + 1 + II - JJ][JJ] = AB[KV + 1 + IP - JJ][JJ];
              AB[KV + 1 + IP - JJ][JJ] = TEMP;
            }
          }
        }

        // Update the relevant part of the trailing submatrix

        if (J2 > 0) {
          // Update A12

          dtrsm('Left', 'Lower', 'No transpose', 'Unit', JB, J2, ONE,
              AB(KV + 1, J), LDAB - 1, AB(KV + 1 - JB, J + JB), LDAB - 1);

          if (I2 > 0) {
            // Update A22

            dgemm(
                'No transpose',
                'No transpose',
                I2,
                J2,
                JB,
                -ONE,
                AB(KV + 1 + JB, J),
                LDAB - 1,
                AB(KV + 1 - JB, J + JB),
                LDAB - 1,
                ONE,
                AB(KV + 1, J + JB),
                LDAB - 1);
          }

          if (I3 > 0) {
            // Update A32

            dgemm(
                'No transpose',
                'No transpose',
                I3,
                J2,
                JB,
                -ONE,
                WORK31,
                LDWORK,
                AB(KV + 1 - JB, J + JB),
                LDAB - 1,
                ONE,
                AB(KV + KL + 1 - JB, J + JB),
                LDAB - 1);
          }
        }

        if (J3 > 0) {
          // Copy the lower triangle of A13 into the work array
          // WORK13

          for (JJ = 1; JJ <= J3; JJ++) {
            for (II = JJ; II <= JB; II++) {
              WORK13[II][JJ] = AB[II - JJ + 1][JJ + J + KV - 1];
            }
          }

          // Update A13 in the work array

          dtrsm('Left', 'Lower', 'No transpose', 'Unit', JB, J3, ONE,
              AB(KV + 1, J), LDAB - 1, WORK13, LDWORK);

          if (I2 > 0) {
            // Update A23

            dgemm(
                'No transpose',
                'No transpose',
                I2,
                J3,
                JB,
                -ONE,
                AB(KV + 1 + JB, J),
                LDAB - 1,
                WORK13,
                LDWORK,
                ONE,
                AB(1 + JB, J + KV),
                LDAB - 1);
          }

          if (I3 > 0) {
            // Update A33

            dgemm('No transpose', 'No transpose', I3, J3, JB, -ONE, WORK31,
                LDWORK, WORK13, LDWORK, ONE, AB(1 + KL, J + KV), LDAB - 1);
          }

          // Copy the lower triangle of A13 back into place

          for (JJ = 1; JJ <= J3; JJ++) {
            for (II = JJ; II <= JB; II++) {
              AB[II - JJ + 1][JJ + J + KV - 1] = WORK13[II][JJ];
            }
          }
        }
      } else {
        // Adjust the pivot indices.

        for (I = J; I <= J + JB - 1; I++) {
          IPIV[I] += J - 1;
        }
      }

      // Partially undo the interchanges in the current block to
      // restore the upper triangular form of A31 and copy the upper
      // triangle of A31 back into place

      for (JJ = J + JB - 1; JJ >= J; JJ--) {
        JP = IPIV[JJ] - JJ + 1;
        if (JP != 1) {
          // Apply interchange to columns J to JJ-1

          if (JP + JJ - 1 < J + KL) {
            // The interchange does not affect A31

            dswap(JJ - J, AB(KV + 1 + JJ - J, J).asArray(), LDAB - 1,
                AB(KV + JP + JJ - J, J).asArray(), LDAB - 1);
          } else {
            // The interchange does affect A31

            dswap(JJ - J, AB(KV + 1 + JJ - J, J).asArray(), LDAB - 1,
                WORK31(JP + JJ - J - KL, 1).asArray(), LDWORK);
          }
        }

        // Copy the current column of A31 back into place

        NW = min(I3, JJ - J + 1);
        if (NW > 0) {
          dcopy(NW, WORK31(1, JJ - J + 1).asArray(), 1,
              AB(KV + KL + 1 - JJ + J, JJ).asArray(), 1);
        }
      }
    }
  }
}
