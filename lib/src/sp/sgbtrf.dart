      void sgbtrf(M, N, KL, KU, AB, LDAB, IPIV, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDAB, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double               AB( LDAB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NBMAX, LDWORK;
      const              NBMAX = 64, LDWORK = NBMAX+1 ;
      // ..
      // .. Local Scalars ..
      int                I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, JU, K2, KM, KV, NB, NW;
      double               TEMP;
      // ..
      // .. Local Arrays ..
      double               WORK13( LDWORK, NBMAX ), WORK31( LDWORK, NBMAX );
      // ..
      // .. External Functions ..
      //- int                ILAENV, ISAMAX;
      // EXTERNAL ILAENV, ISAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGBTF2, SGEMM, SGER, SLASWP, SSCAL, SSWAP, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // KV is the number of superdiagonals in the factor U, allowing for
      // fill-in

      KV = KU + KL;

      // Test the input parameters.

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KL < 0 ) {
         INFO = -3;
      } else if ( KU < 0 ) {
         INFO = -4;
      } else if ( LDAB < KL+KV+1 ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('SGBTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      // Determine the block size for this environment

      NB = ilaenv( 1, 'SGBTRF', ' ', M, N, KL, KU );

      // The block size must not exceed the limit set by the size of the
      // local arrays WORK13 and WORK31.

      NB = min( NB, NBMAX );

      if ( NB <= 1 || NB > KL ) {

         // Use unblocked code

         sgbtf2(M, N, KL, KU, AB, LDAB, IPIV, INFO );
      } else {

         // Use blocked code

         // Zero the superdiagonal elements of the work array WORK13

         for (J = 1; J <= NB; J++) { // 20
            for (I = 1; I <= J - 1; I++) { // 10
               WORK13[I][J] = ZERO;
            } // 10
         } // 20

         // Zero the subdiagonal elements of the work array WORK31

         for (J = 1; J <= NB; J++) { // 40
            for (I = J + 1; I <= NB; I++) { // 30
               WORK31[I][J] = ZERO;
            } // 30
         } // 40

         // Gaussian elimination with partial pivoting

         // Set fill-in elements in columns KU+2 to KV to zero

         for (J = KU + 2; J <= min( KV, N ); J++) { // 60
            for (I = KV - J + 2; I <= KL; I++) { // 50
               AB[I][J] = ZERO;
            } // 50
         } // 60

         // JU is the index of the last column affected by the current
         // stage of the factorization

         JU = 1;

         for (J = 1; NB < 0 ? J >= min( M, N ) : J <= min( M, N ); J += NB) { // 180
            JB = min( NB, min( M, N )-J+1 );

            // The active part of the matrix is partitioned

               // A11   A12   A13
               // A21   A22   A23
               // A31   A32   A33

            // Here A11, A21 and A31 denote the current block of JB columns
            // which is about to be factorized. The number of rows in the
            // partitioning are JB, I2, I3 respectively, and the numbers
            // of columns are JB, J2, J3. The superdiagonal elements of A13
            // and the subdiagonal elements of A31 lie outside the band.

            I2 = min( KL-JB, M-J-JB+1 );
            I3 = min( JB, M-J-KL+1 );

            // J2 and J3 are computed after JU has been updated.

            // Factorize the current block of JB columns

            for (JJ = J; JJ <= J + JB - 1; JJ++) { // 80

               // Set fill-in elements in column JJ+KV to zero

               if ( JJ+KV <= N ) {
                  for (I = 1; I <= KL; I++) { // 70
                     AB[I, JJ+KV] = ZERO;
                  } // 70
               }

               // Find pivot and test for singularity. KM is the number of
               // subdiagonal elements in the current column.

               KM = min( KL, M-JJ );
               JP = ISAMAX( KM+1, AB( KV+1, JJ ), 1 );
               IPIV[JJ] = JP + JJ - J;
               if ( AB( KV+JP, JJ ) != ZERO ) {
                  JU = max( JU, min( JJ+KU+JP-1, N ) );
                  if ( JP != 1 ) {

                     // Apply interchange to columns J to J+JB-1

                     if ( JP+JJ-1 < J+KL ) {

                        sswap(JB, AB( KV+1+JJ-J, J ), LDAB-1, AB( KV+JP+JJ-J, J ), LDAB-1 );
                     } else {

                        // The interchange affects columns J to JJ-1 of A31
                        // which are stored in the work array WORK31

                        sswap(JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, WORK31( JP+JJ-J-KL, 1 ), LDWORK );
                        sswap(J+JB-JJ, AB( KV+1, JJ ), LDAB-1, AB( KV+JP, JJ ), LDAB-1 );
                     }
                  }

                  // Compute multipliers

                  sscal(KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), 1 );

                  // Update trailing submatrix within the band and within
                  // the current block. JM is the index of the last column
                  // which needs to be updated.

                  JM = min( JU, J+JB-1 );
                  if (JM > JJ) sger( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, AB( KV, JJ+1 ), LDAB-1, AB( KV+1, JJ+1 ), LDAB-1 );
               } else {

                  // If pivot is zero, set INFO to the index of the pivot
                  // unless a zero pivot has already been found.

                  if (INFO == 0) INFO = JJ;
               }

               // Copy current column of A31 into the work array WORK31

               NW = min( JJ-J+1, I3 );
               if (NW > 0) scopy( NW, AB( KV+KL+1-JJ+J, JJ ), 1, WORK31( 1, JJ-J+1 ), 1 );
            } // 80
            if ( J+JB <= N ) {

               // Apply the row interchanges to the other blocks.

               J2 = min( JU-J+1, KV ) - JB;
               J3 = max( 0, JU-J-KV+1 );

               // Use SLASWP to apply the row interchanges to A12, A22, and
               // A32.

               slaswp(J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, IPIV( J ), 1 );

               // Adjust the pivot indices.

               for (I = J; I <= J + JB - 1; I++) { // 90
                  IPIV[I] = IPIV( I ) + J - 1;
               } // 90

               // Apply the row interchanges to A13, A23, and A33
               // columnwise.

               K2 = J - 1 + JB + J2;
               for (I = 1; I <= J3; I++) { // 110
                  JJ = K2 + I;
                  for (II = J + I - 1; II <= J + JB - 1; II++) { // 100
                     IP = IPIV( II );
                     if ( IP != II ) {
                        TEMP = AB( KV+1+II-JJ, JJ );
                        AB[KV+1+II-JJ, JJ] = AB( KV+1+IP-JJ, JJ );
                        AB[KV+1+IP-JJ, JJ] = TEMP;
                     }
                  } // 100
               } // 110

               // Update the relevant part of the trailing submatrix

               if ( J2 > 0 ) {

                  // Update A12

                  strsm('Left', 'Lower', 'No transpose', 'Unit', JB, J2, ONE, AB( KV+1, J ), LDAB-1, AB( KV+1-JB, J+JB ), LDAB-1 );

                  if ( I2 > 0 ) {

                     // Update A22

                     sgemm('No transpose', 'No transpose', I2, J2, JB, -ONE, AB( KV+1+JB, J ), LDAB-1, AB( KV+1-JB, J+JB ), LDAB-1, ONE, AB( KV+1, J+JB ), LDAB-1 );
                  }

                  if ( I3 > 0 ) {

                     // Update A32

                     sgemm('No transpose', 'No transpose', I3, J2, JB, -ONE, WORK31, LDWORK, AB( KV+1-JB, J+JB ), LDAB-1, ONE, AB( KV+KL+1-JB, J+JB ), LDAB-1 );
                  }
               }

               if ( J3 > 0 ) {

                  // Copy the lower triangle of A13 into the work array
                  // WORK13

                  for (JJ = 1; JJ <= J3; JJ++) { // 130
                     for (II = JJ; II <= JB; II++) { // 120
                        WORK13[II][JJ] = AB( II-JJ+1, JJ+J+KV-1 );
                     } // 120
                  } // 130

                  // Update A13 in the work array

                  strsm('Left', 'Lower', 'No transpose', 'Unit', JB, J3, ONE, AB( KV+1, J ), LDAB-1, WORK13, LDWORK );

                  if ( I2 > 0 ) {

                     // Update A23

                     sgemm('No transpose', 'No transpose', I2, J3, JB, -ONE, AB( KV+1+JB, J ), LDAB-1, WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), LDAB-1 );
                  }

                  if ( I3 > 0 ) {

                     // Update A33

                     sgemm('No transpose', 'No transpose', I3, J3, JB, -ONE, WORK31, LDWORK, WORK13, LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 );
                  }

                  // Copy the lower triangle of A13 back into place

                  for (JJ = 1; JJ <= J3; JJ++) { // 150
                     for (II = JJ; II <= JB; II++) { // 140
                        AB[II-JJ+1, JJ+J+KV-1] = WORK13( II, JJ );
                     } // 140
                  } // 150
               }
            } else {

               // Adjust the pivot indices.

               for (I = J; I <= J + JB - 1; I++) { // 160
                  IPIV[I] = IPIV( I ) + J - 1;
               } // 160
            }

            // Partially undo the interchanges in the current block to
            // restore the upper triangular form of A31 and copy the upper
            // triangle of A31 back into place

            for (JJ = J + JB - 1; JJ >= J; JJ--) { // 170
               JP = IPIV( JJ ) - JJ + 1;
               if ( JP != 1 ) {

                  // Apply interchange to columns J to JJ-1

                  if ( JP+JJ-1 < J+KL ) {

                     // The interchange does not affect A31

                     sswap(JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, AB( KV+JP+JJ-J, J ), LDAB-1 );
                  } else {

                     // The interchange does affect A31

                     sswap(JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, WORK31( JP+JJ-J-KL, 1 ), LDWORK );
                  }
               }

               // Copy the current column of A31 back into place

               NW = min( I3, JJ-J+1 );
               if (NW > 0) scopy( NW, WORK31( 1, JJ-J+1 ), 1, AB( KV+KL+1-JJ+J, JJ ), 1 );
            } // 170
         } // 180
      }

      return;
      }
