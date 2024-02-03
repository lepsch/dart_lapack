      SUBROUTINE DLAORHR_COL_GETRFNP( M, N, A, LDA, D, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), D( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                IINFO, J, JB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLAORHR_COL_GETRFNP2, DTRSM, XERBLA
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
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DLAORHR_COL_GETRFNP', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( MIN( M, N ).EQ.0 ) RETURN

      // Determine the block size for this environment.


      NB = ILAENV( 1, 'DLAORHR_COL_GETRFNP', ' ', M, N, -1, -1 )

      if ( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) {

         // Use unblocked code.

         CALL DLAORHR_COL_GETRFNP2( M, N, A, LDA, D, INFO )
      } else {

         // Use blocked code.

         DO J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )

            // Factor diagonal and subdiagonal blocks.

            CALL DLAORHR_COL_GETRFNP2( M-J+1, JB, A( J, J ), LDA, D( J ), IINFO )

            if ( J+JB.LE.N ) {

               // Compute block row of U.

               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA )
               if ( J+JB.LE.M ) {

                  // Update trailing submatrix.

                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA )
               }
            }
         END DO
      }
      RETURN

      // End of DLAORHR_COL_GETRFNP

      }
