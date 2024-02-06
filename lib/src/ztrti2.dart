      void ztrti2(UPLO, DIAG, N, A, LDA, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, UPLO;
      int                INFO, LDA, N;
      Complex         A( LDA, * );
      // ..

      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               NOUNIT, UPPER;
      int                J;
      Complex         AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZSCAL, ZTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

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
         xerbla('ZTRTI2', -INFO );
         return;
      }

      if ( UPPER ) {

         // Compute inverse of upper triangular matrix.

         for (J = 1; J <= N; J++) { // 10
            if ( NOUNIT ) {
               A[J][J] = ONE / A( J, J );
               AJJ = -A( J, J );
            } else {
               AJJ = -ONE;
            }

            // Compute elements 1:j-1 of j-th column.

            ztrmv('Upper', 'No transpose', DIAG, J-1, A, LDA, A( 1, J ), 1 );
            zscal(J-1, AJJ, A( 1, J ), 1 );
         } // 10
      } else {

         // Compute inverse of lower triangular matrix.

         for (J = N; J >= 1; J--) { // 20
            if ( NOUNIT ) {
               A[J][J] = ONE / A( J, J );
               AJJ = -A( J, J );
            } else {
               AJJ = -ONE;
            }
            if ( J < N ) {

               // Compute elements j+1:n of j-th column.

               ztrmv('Lower', 'No transpose', DIAG, N-J, A( J+1, J+1 ), LDA, A( J+1, J ), 1 );
               zscal(N-J, AJJ, A( J+1, J ), 1 );
            }
         } // 20
      }

      return;
      }
