      void dlauu2(UPLO, N, A, LDA, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      double             AII;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- double             DDOT;
      // EXTERNAL LSAME, DDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DLAUU2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Compute the product U * U**T.

         for (I = 1; I <= N; I++) { // 10
            AII = A( I, I );
            if ( I < N ) {
               A( I, I ) = DDOT( N-I+1, A( I, I ), LDA, A( I, I ), LDA );
               dgemv('No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, AII, A( 1, I ), 1 );
            } else {
               dscal(I, AII, A( 1, I ), 1 );
            }
         } // 10

      } else {

         // Compute the product L**T * L.

         for (I = 1; I <= N; I++) { // 20
            AII = A( I, I );
            if ( I < N ) {
               A( I, I ) = DDOT( N-I+1, A( I, I ), 1, A( I, I ), 1 );
               dgemv('Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, AII, A( I, 1 ), LDA );
            } else {
               dscal(I, AII, A( I, 1 ), LDA );
            }
         } // 20
      }

      return;
      }
