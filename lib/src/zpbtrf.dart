import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zpbtf2.dart';
import 'package:lapack/src/zpotf2.dart';

void zpbtrf(
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.dim(LDAB);
  const ONE = 1.0;
  const NBMAX = 32, LDWORK = NBMAX + 1;
  int I, I2, I3, IB, J, JJ, NB;
  final WORK = Matrix<Complex>(LDWORK, NBMAX);
  final II = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  if ((!lsame(UPLO, 'U')) && (!lsame(UPLO, 'L'))) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KD < 0) {
    INFO.value = -3;
  } else if (LDAB < KD + 1) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('ZPBTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Determine the block size for this environment

  NB = ilaenv(1, 'ZPBTRF', UPLO, N, KD, -1, -1);

  // The block size must not exceed the semi-bandwidth KD, and must not
  // exceed the limit set by the size of the local array WORK.

  NB = min(NB, NBMAX);

  if (NB <= 1 || NB > KD) {
    // Use unblocked code

    zpbtf2(UPLO, N, KD, AB, LDAB, INFO);
  } else {
    // Use blocked code

    if (lsame(UPLO, 'U')) {
      // Compute the Cholesky factorization of a Hermitian band
      // matrix, given the upper triangle of the matrix in band
      // storage.

      // Zero the upper triangle of the work array.

      for (J = 1; J <= NB; J++) {
        // 20
        for (I = 1; I <= J - 1; I++) {
          // 10
          WORK[I][J] = Complex.zero;
        } // 10
      } // 20

      // Process the band matrix one diagonal block at a time.

      for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) {
        // 70
        IB = min(NB, N - I + 1);

        // Factorize the diagonal block

        zpotf2(UPLO, IB, AB(KD + 1, I), LDAB - 1, II);
        if (II.value != 0) {
          INFO.value = I + II.value - 1;
          return;
        }
        if (I + IB <= N) {
          // Update the relevant part of the trailing submatrix.
          // If A11 denotes the diagonal block which has just been
          // factorized, then we need to update the remaining
          // blocks in the diagram:

          // A11   A12   A13
          //       A22   A23
          //             A33

          // The numbers of rows and columns in the partitioning
          // are IB, I2, I3 respectively. The blocks A12, A22 and
          // A23 are empty if IB = KD. The upper triangle of A13
          // lies outside the band.

          I2 = min(KD - IB, N - I - IB + 1);
          I3 = min(IB, N - I - KD + 1);

          if (I2 > 0) {
            // Update A12

            ztrsm(
                'Left',
                'Upper',
                'Conjugate transpose',
                'Non-unit',
                IB,
                I2,
                Complex.one,
                AB(KD + 1, I),
                LDAB - 1,
                AB(KD + 1 - IB, I + IB),
                LDAB - 1);

            // Update A22

            zherk(
                'Upper',
                'Conjugate transpose',
                I2,
                IB,
                -ONE,
                AB(KD + 1 - IB, I + IB),
                LDAB - 1,
                ONE,
                AB(KD + 1, I + IB),
                LDAB - 1);
          }

          if (I3 > 0) {
            // Copy the lower triangle of A13 into the work array.

            for (JJ = 1; JJ <= I3; JJ++) {
              // 40
              for (II.value = JJ; II.value <= IB; II.value++) {
                // 30
                WORK[II.value][JJ] = AB[II.value - JJ + 1][JJ + I + KD - 1];
              } // 30
            } // 40

            // Update A13 (in the work array).

            ztrsm('Left', 'Upper', 'Conjugate transpose', 'Non-unit', IB, I3,
                Complex.one, AB(KD + 1, I), LDAB - 1, WORK, LDWORK);

            // Update A23

            if (I2 > 0) {
              zgemm(
                  'Conjugate transpose',
                  'No transpose',
                  I2,
                  I3,
                  IB,
                  -Complex.one,
                  AB(KD + 1 - IB, I + IB),
                  LDAB - 1,
                  WORK,
                  LDWORK,
                  Complex.one,
                  AB(1 + IB, I + KD),
                  LDAB - 1);
            }

            // Update A33

            zherk('Upper', 'Conjugate transpose', I3, IB, -ONE, WORK, LDWORK,
                ONE, AB(KD + 1, I + KD), LDAB - 1);

            // Copy the lower triangle of A13 back into place.

            for (JJ = 1; JJ <= I3; JJ++) {
              // 60
              for (II.value = JJ; II.value <= IB; II.value++) {
                // 50
                AB[II.value - JJ + 1][JJ + I + KD - 1] = WORK[II.value][JJ];
              } // 50
            } // 60
          }
        }
      } // 70
    } else {
      // Compute the Cholesky factorization of a Hermitian band
      // matrix, given the lower triangle of the matrix in band
      // storage.

      // Zero the lower triangle of the work array.

      for (J = 1; J <= NB; J++) {
        // 90
        for (I = J + 1; I <= NB; I++) {
          // 80
          WORK[I][J] = Complex.zero;
        } // 80
      } // 90

      // Process the band matrix one diagonal block at a time.

      for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) {
        // 140
        IB = min(NB, N - I + 1);

        // Factorize the diagonal block

        zpotf2(UPLO, IB, AB(1, I), LDAB - 1, II);
        if (II.value != 0) {
          INFO.value = I + II.value - 1;
          return;
        }
        if (I + IB <= N) {
          // Update the relevant part of the trailing submatrix.
          // If A11 denotes the diagonal block which has just been
          // factorized, then we need to update the remaining
          // blocks in the diagram:

          // A11
          // A21   A22
          // A31   A32   A33

          // The numbers of rows and columns in the partitioning
          // are IB, I2, I3 respectively. The blocks A21, A22 and
          // A32 are empty if IB = KD. The lower triangle of A31
          // lies outside the band.

          I2 = min(KD - IB, N - I - IB + 1);
          I3 = min(IB, N - I - KD + 1);

          if (I2 > 0) {
            // Update A21

            ztrsm('Right', 'Lower', 'Conjugate transpose', 'Non-unit', I2, IB,
                Complex.one, AB(1, I), LDAB - 1, AB(1 + IB, I), LDAB - 1);

            // Update A22

            zherk('Lower', 'No transpose', I2, IB, -ONE, AB(1 + IB, I),
                LDAB - 1, ONE, AB(1, I + IB), LDAB - 1);
          }

          if (I3 > 0) {
            // Copy the upper triangle of A31 into the work array.

            for (JJ = 1; JJ <= IB; JJ++) {
              // 110
              for (II.value = 1; II.value <= min(JJ, I3); II.value++) {
                // 100
                WORK[II.value][JJ] = AB[KD + 1 - JJ + II.value][JJ + I - 1];
              } // 100
            } // 110

            // Update A31 (in the work array).

            ztrsm('Right', 'Lower', 'Conjugate transpose', 'Non-unit', I3, IB,
                Complex.one, AB(1, I), LDAB - 1, WORK, LDWORK);

            // Update A32

            if (I2 > 0) {
              zgemm(
                  'No transpose',
                  'Conjugate transpose',
                  I3,
                  I2,
                  IB,
                  -Complex.one,
                  WORK,
                  LDWORK,
                  AB(1 + IB, I),
                  LDAB - 1,
                  Complex.one,
                  AB(1 + KD - IB, I + IB),
                  LDAB - 1);
            }

            // Update A33

            zherk('Lower', 'No transpose', I3, IB, -ONE, WORK, LDWORK, ONE,
                AB(1, I + KD), LDAB - 1);

            // Copy the upper triangle of A31 back into place.

            for (JJ = 1; JJ <= IB; JJ++) {
              // 130
              for (II.value = 1; II.value <= min(JJ, I3); II.value++) {
                // 120
                AB[KD + 1 - JJ + II.value][JJ + I - 1] = WORK[II.value][JJ];
              } // 120
            } // 130
          }
        }
      } // 140
    }
  }
}
