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
      const              ONE = 1.0E+0 ;
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
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('SLAORHR_COL_GETRFNP', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( MIN( M, N ) == 0 ) RETURN

      // Determine the block size for this environment.


      NB = ILAENV( 1, 'SLAORHR_COL_GETRFNP', ' ', M, N, -1, -1 )

      if ( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) {

         // Use unblocked code.

         slaorhr_col_getrfnp2(M, N, A, LDA, D, INFO );
      } else {

         // Use blocked code.

         DO J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )

            // Factor diagonal and subdiagonal blocks.

            slaorhr_col_getrfnp2(M-J+1, JB, A( J, J ), LDA, D( J ), IINFO );

            if ( J+JB.LE.N ) {

               // Compute block row of U.

               strsm('Left', 'Lower', 'No transpose', 'Unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA );
               if ( J+JB.LE.M ) {

                  // Update trailing submatrix.

                  sgemm('No transpose', 'No transpose', M-J-JB+1, N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA );
               }
            }
         }
      }
      RETURN

      // End of SLAORHR_COL_GETRFNP

      }
