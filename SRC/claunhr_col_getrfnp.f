      SUBROUTINE CLAUNHR_COL_GETRFNP( M, N, A, LDA, D, INFO );
      // IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), D( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                IINFO, J, JB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLAUNHR_COL_GETRFNP2, CTRSM, XERBLA
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

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CLAUNHR_COL_GETRFNP', -INFO );
         RETURN;
      }

      // Quick return if possible

      IF( MIN( M, N ) == 0 ) RETURN;

      // Determine the block size for this environment.


      NB = ILAENV( 1, 'CLAUNHR_COL_GETRFNP', ' ', M, N, -1, -1 );

      if ( NB <= 1 || NB >= MIN( M, N ) ) {

         // Use unblocked code.

         claunhr_col_getrfnp2(M, N, A, LDA, D, INFO );
      } else {

         // Use blocked code.

         DO J = 1, MIN( M, N ), NB;
            JB = MIN( MIN( M, N )-J+1, NB );

            // Factor diagonal and subdiagonal blocks.

            claunhr_col_getrfnp2(M-J+1, JB, A( J, J ), LDA, D( J ), IINFO );

            if ( J+JB <= N ) {

               // Compute block row of U.

               ctrsm('Left', 'Lower', 'No transpose', 'Unit', JB, N-J-JB+1, CONE, A( J, J ), LDA, A( J, J+JB ), LDA );
               if ( J+JB <= M ) {

                  // Update trailing submatrix.

                  cgemm('No transpose', 'No transpose', M-J-JB+1, N-J-JB+1, JB, -CONE, A( J+JB, J ), LDA, A( J, J+JB ), LDA, CONE, A( J+JB, J+JB ), LDA );
               }
            }
         }
      }
      RETURN;

      // End of CLAUNHR_COL_GETRFNP

      }
