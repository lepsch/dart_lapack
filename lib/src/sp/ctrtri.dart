      void ctrtri(UPLO, DIAG, N, A, LDA, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      Complex            A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex            ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J, JB, NB, NN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTRMM, CTRSM, CTRTI2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      NOUNIT = lsame( DIAG, 'N' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('CTRTRI', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Check for singularity if non-unit.

      if ( NOUNIT ) {
         for (INFO = 1; INFO <= N; INFO++) { // 10
            if( A( INFO, INFO ) == ZERO ) return;
         } // 10
         INFO = 0;
      }

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'CTRTRI', UPLO // DIAG, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code

         ctrti2(UPLO, DIAG, N, A, LDA, INFO );
      } else {

         // Use blocked code

         if ( UPPER ) {

            // Compute inverse of upper triangular matrix

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 20
               JB = min( NB, N-J+1 );

               // Compute rows 1:j-1 of current block column

               ctrmm('Left', 'Upper', 'No transpose', DIAG, J-1, JB, ONE, A, LDA, A( 1, J ), LDA );
               ctrsm('Right', 'Upper', 'No transpose', DIAG, J-1, JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA );

               // Compute inverse of current diagonal block

               ctrti2('Upper', DIAG, JB, A( J, J ), LDA, INFO );
            } // 20
         } else {

            // Compute inverse of lower triangular matrix

            NN = ( ( N-1 ) / NB )*NB + 1;
            for (J = NN; -NB < 0 ? J >= 1 : J <= 1; J += -NB) { // 30
               JB = min( NB, N-J+1 );
               if ( J+JB <= N ) {

                  // Compute rows j+jb:n of current block column

                  ctrmm('Left', 'Lower', 'No transpose', DIAG, N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA, A( J+JB, J ), LDA );
                  ctrsm('Right', 'Lower', 'No transpose', DIAG, N-J-JB+1, JB, -ONE, A( J, J ), LDA, A( J+JB, J ), LDA );
               }

               // Compute inverse of current diagonal block

               ctrti2('Lower', DIAG, JB, A( J, J ), LDA, INFO );
            } // 30
         }
      }

      return;
      }
