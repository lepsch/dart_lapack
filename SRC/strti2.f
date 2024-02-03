      void strti2(UPLO, DIAG, N, A, LDA, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J;
      REAL               AJJ;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, STRMV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      NOUNIT = LSAME( DIAG, 'N' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOUNIT && !LSAME( DIAG, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('STRTI2', -INFO );
         return;
      }

      if ( UPPER ) {

         // Compute inverse of upper triangular matrix.

         for (J = 1; J <= N; J++) { // 10
            if ( NOUNIT ) {
               A( J, J ) = ONE / A( J, J );
               AJJ = -A( J, J );
            } else {
               AJJ = -ONE;
            }

            // Compute elements 1:j-1 of j-th column.

            strmv('Upper', 'No transpose', DIAG, J-1, A, LDA, A( 1, J ), 1 );
            sscal(J-1, AJJ, A( 1, J ), 1 );
         } // 10
      } else {

         // Compute inverse of lower triangular matrix.

         DO 20 J = N, 1, -1;
            if ( NOUNIT ) {
               A( J, J ) = ONE / A( J, J );
               AJJ = -A( J, J );
            } else {
               AJJ = -ONE;
            }
            if ( J < N ) {

               // Compute elements j+1:n of j-th column.

               strmv('Lower', 'No transpose', DIAG, N-J, A( J+1, J+1 ), LDA, A( J+1, J ), 1 );
               sscal(N-J, AJJ, A( J+1, J ), 1 );
            }
         } // 20
      }

      return;
      }
