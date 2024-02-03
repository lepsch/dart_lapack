      SUBROUTINE SLAORHR_COL_GETRFNP( M, N, A, LDA, D, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), D( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                IINFO, J, JB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLAORHR_COL_GETRFNP2, STRSM, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAORHR_COL_GETRFNP', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( MIN( M, N ).EQ.0 ) RETURN

      // Determine the block size for this environment.


      NB = ILAENV( 1, 'SLAORHR_COL_GETRFNP', ' ', M, N, -1, -1 )

      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN

         // Use unblocked code.

         CALL SLAORHR_COL_GETRFNP2( M, N, A, LDA, D, INFO )
      ELSE

         // Use blocked code.

         DO J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )

            // Factor diagonal and subdiagonal blocks.

            CALL SLAORHR_COL_GETRFNP2( M-J+1, JB, A( J, J ), LDA, D( J ), IINFO )

            IF( J+JB.LE.N ) THEN

               // Compute block row of U.

               CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA )
               IF( J+JB.LE.M ) THEN

                  // Update trailing submatrix.

                  CALL SGEMM( 'No transpose', 'No transpose', M-J-JB+1, N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA )
               END IF
            END IF
         END DO
      END IF
      RETURN

      // End of SLAORHR_COL_GETRFNP

      END
