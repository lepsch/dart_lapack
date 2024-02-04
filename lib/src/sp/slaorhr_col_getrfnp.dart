      void slaorhr_col_getrfnp(M, N, A, LDA, D, INFO ) {
      // IMPLICIT NONE

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * ), D( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                IINFO, J, JB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLAORHR_COL_GETRFNP2, STRSM, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
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
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('SLAORHR_COL_GETRFNP', -INFO );
         return;
      }

      // Quick return if possible

      if( min( M, N ) == 0 ) return;

      // Determine the block size for this environment.


      NB = ILAENV( 1, 'SLAORHR_COL_GETRFNP', ' ', M, N, -1, -1 );

      if ( NB <= 1 || NB >= min( M, N ) ) {

         // Use unblocked code.

         slaorhr_col_getrfnp2(M, N, A, LDA, D, INFO );
      } else {

         // Use blocked code.

         for (J = 1; NB < 0 ? J >= min( M, N ) : J <= min( M, N ); J += NB) {
            JB = min( min( M, N )-J+1, NB );

            // Factor diagonal and subdiagonal blocks.

            slaorhr_col_getrfnp2(M-J+1, JB, A( J, J ), LDA, D( J ), IINFO );

            if ( J+JB <= N ) {

               // Compute block row of U.

               strsm('Left', 'Lower', 'No transpose', 'Unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA );
               if ( J+JB <= M ) {

                  // Update trailing submatrix.

                  sgemm('No transpose', 'No transpose', M-J-JB+1, N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA );
               }
            }
         }
      }
      return;
      }
