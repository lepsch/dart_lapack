      SUBROUTINE CGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDAB, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            AB( LDAB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO
      const              ONE = ( 1.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ) ;
      int                NBMAX, LDWORK;
      const              NBMAX = 64, LDWORK = NBMAX+1 ;
      // ..
      // .. Local Scalars ..
      int                I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, JU, K2, KM, KV, NB, NW;
      COMPLEX            TEMP
      // ..
      // .. Local Arrays ..
      COMPLEX            WORK13( LDWORK, NBMAX ), WORK31( LDWORK, NBMAX )
      // ..
      // .. External Functions ..
      int                ICAMAX, ILAENV;
      // EXTERNAL ICAMAX, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGBTF2, CGEMM, CGERU, CLASWP, CSCAL, CSWAP, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // KV is the number of superdiagonals in the factor U, allowing for
      // fill-in

      KV = KU + KL

      // Test the input parameters.

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KL.LT.0 ) {
         INFO = -3
      } else if ( KU.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.KL+KV+1 ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CGBTRF', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Determine the block size for this environment

      NB = ILAENV( 1, 'CGBTRF', ' ', M, N, KL, KU )

      // The block size must not exceed the limit set by the size of the
      // local arrays WORK13 and WORK31.

      NB = MIN( NB, NBMAX )

      if ( NB.LE.1 .OR. NB.GT.KL ) {

         // Use unblocked code

         CALL CGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      } else {

         // Use blocked code

         // Zero the superdiagonal elements of the work array WORK13

         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE

         // Zero the subdiagonal elements of the work array WORK31

         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE

         // Gaussian elimination with partial pivoting

         // Set fill-in elements in columns KU+2 to KV to zero

         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE

         // JU is the index of the last column affected by the current
         // stage of the factorization

         JU = 1

         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )

            // The active part of the matrix is partitioned

               // A11   A12   A13
               // A21   A22   A23
               // A31   A32   A33

            // Here A11, A21 and A31 denote the current block of JB columns
            // which is about to be factorized. The number of rows in the
            // partitioning are JB, I2, I3 respectively, and the numbers
            // of columns are JB, J2, J3. The superdiagonal elements of A13
            // and the subdiagonal elements of A31 lie outside the band.

            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )

            // J2 and J3 are computed after JU has been updated.

            // Factorize the current block of JB columns

            DO 80 JJ = J, J + JB - 1

               // Set fill-in elements in column JJ+KV to zero

               if ( JJ+KV.LE.N ) {
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               }

               // Find pivot and test for singularity. KM is the number of
               // subdiagonal elements in the current column.

               KM = MIN( KL, M-JJ )
               JP = ICAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               if ( AB( KV+JP, JJ ).NE.ZERO ) {
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  if ( JP.NE.1 ) {

                     // Apply interchange to columns J to J+JB-1

                     if ( JP+JJ-1.LT.J+KL ) {

                        CALL CSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, AB( KV+JP+JJ-J, J ), LDAB-1 )
                     } else {

                        // The interchange affects columns J to JJ-1 of A31
                        // which are stored in the work array WORK31

                        CALL CSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, WORK31( JP+JJ-J-KL, 1 ), LDWORK )                         CALL CSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, AB( KV+JP, JJ ), LDAB-1 )
                     }
                  }

                  // Compute multipliers

                  CALL CSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), 1 )

                  // Update trailing submatrix within the band and within
                  // the current block. JM is the index of the last column
                  // which needs to be updated.

                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ ) CALL CGERU( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, AB( KV, JJ+1 ), LDAB-1, AB( KV+1, JJ+1 ), LDAB-1 )
               } else {

                  // If pivot is zero, set INFO to the index of the pivot
                  // unless a zero pivot has already been found.

                  IF( INFO.EQ.0 ) INFO = JJ
               }

               // Copy current column of A31 into the work array WORK31

               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 ) CALL CCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            if ( J+JB.LE.N ) {

               // Apply the row interchanges to the other blocks.

               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )

               // Use CLASWP to apply the row interchanges to A12, A22, and
               // A32.

               CALL CLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, IPIV( J ), 1 )

               // Adjust the pivot indices.

               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE

               // Apply the row interchanges to A13, A23, and A33
               // columnwise.

               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     if ( IP.NE.II ) {
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     }
  100             CONTINUE
  110          CONTINUE

               // Update the relevant part of the trailing submatrix

               if ( J2.GT.0 ) {

                  // Update A12

                  CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, J2, ONE, AB( KV+1, J ), LDAB-1, AB( KV+1-JB, J+JB ), LDAB-1 )

                  if ( I2.GT.0 ) {

                     // Update A22

                     CALL CGEMM( 'No transpose', 'No transpose', I2, J2, JB, -ONE, AB( KV+1+JB, J ), LDAB-1, AB( KV+1-JB, J+JB ), LDAB-1, ONE, AB( KV+1, J+JB ), LDAB-1 )
                  }

                  if ( I3.GT.0 ) {

                     // Update A32

                     CALL CGEMM( 'No transpose', 'No transpose', I3, J2, JB, -ONE, WORK31, LDWORK, AB( KV+1-JB, J+JB ), LDAB-1, ONE, AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  }
               }

               if ( J3.GT.0 ) {

                  // Copy the lower triangle of A13 into the work array
                  // WORK13

                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE

                  // Update A13 in the work array

                  CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, J3, ONE, AB( KV+1, J ), LDAB-1, WORK13, LDWORK )

                  if ( I2.GT.0 ) {

                     // Update A23

                     CALL CGEMM( 'No transpose', 'No transpose', I2, J3, JB, -ONE, AB( KV+1+JB, J ), LDAB-1, WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), LDAB-1 )
                  }

                  if ( I3.GT.0 ) {

                     // Update A33

                     CALL CGEMM( 'No transpose', 'No transpose', I3, J3, JB, -ONE, WORK31, LDWORK, WORK13, LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  }

                  // Copy the lower triangle of A13 back into place

                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               }
            } else {

               // Adjust the pivot indices.

               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            }

            // Partially undo the interchanges in the current block to
            // restore the upper triangular form of A31 and copy the upper
            // triangle of A31 back into place

            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               if ( JP.NE.1 ) {

                  // Apply interchange to columns J to JJ-1

                  if ( JP+JJ-1.LT.J+KL ) {

                     // The interchange does not affect A31

                     CALL CSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, AB( KV+JP+JJ-J, J ), LDAB-1 )
                  } else {

                     // The interchange does affect A31

                     CALL CSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  }
               }

               // Copy the current column of A31 back into place

               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 ) CALL CCOPY( NW, WORK31( 1, JJ-J+1 ), 1, AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      }

      RETURN

      // End of CGBTRF

      }
