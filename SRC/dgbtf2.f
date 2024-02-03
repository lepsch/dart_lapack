      SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDAB, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             AB( LDAB, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, JP, JU, KM, KV;
      // ..
      // .. External Functions ..
      int                IDAMAX;
      // EXTERNAL IDAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGER, DSCAL, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // KV is the number of superdiagonals in the factor U, allowing for
      // fill-in.

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
         CALL XERBLA( 'DGBTF2', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Gaussian elimination with partial pivoting

      // Set fill-in elements in columns KU+2 to KV to zero.

      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE

      // JU is the index of the last column affected by the current stage
      // of the factorization.

      JU = 1

      DO 40 J = 1, MIN( M, N )

         // Set fill-in elements in column J+KV to zero.

         if ( J+KV.LE.N ) {
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         }

         // Find pivot and test for singularity. KM is the number of
         // subdiagonal elements in the current column.

         KM = MIN( KL, M-J )
         JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         if ( AB( KV+JP, J ).NE.ZERO ) {
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )

            // Apply interchange to columns J to JU.

            IF( JP.NE.1 ) CALL DSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1, AB( KV+1, J ), LDAB-1 )

            if ( KM.GT.0 ) {

               // Compute multipliers.

               CALL DSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )

               // Update trailing submatrix within the band.

               IF( JU.GT.J ) CALL DGER( KM, JU-J, -ONE, AB( KV+2, J ), 1, AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), LDAB-1 )
            }
         } else {

            // If pivot is zero, set INFO to the index of the pivot
            // unless a zero pivot has already been found.

            IF( INFO.EQ.0 ) INFO = J
         }
   40 CONTINUE
      RETURN

      // End of DGBTF2

      }
