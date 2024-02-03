      SUBROUTINE ZLAUU2( UPLO, N, A, LDA, INFO );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      double             AII;
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX*16         ZDOTC;
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZGEMV, ZLACGV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('ZLAUU2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( UPPER ) {

         // Compute the product U * U**H.

         for (I = 1; I <= N; I++) { // 10
            AII = DBLE( A( I, I ) );
            if ( I < N ) {
               A( I, I ) = AII*AII + DBLE( ZDOTC( N-I, A( I, I+1 ), LDA, A( I, I+1 ), LDA ) );
               zlacgv(N-I, A( I, I+1 ), LDA );
               zgemv('No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, DCMPLX( AII ), A( 1, I ), 1 );
               zlacgv(N-I, A( I, I+1 ), LDA );
            } else {
               zdscal(I, AII, A( 1, I ), 1 );
            }
         } // 10

      } else {

         // Compute the product L**H * L.

         for (I = 1; I <= N; I++) { // 20
            AII = DBLE( A( I, I ) );
            if ( I < N ) {
               A( I, I ) = AII*AII + DBLE( ZDOTC( N-I, A( I+1, I ), 1, A( I+1, I ), 1 ) );
               zlacgv(I-1, A( I, 1 ), LDA );
               zgemv('Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, DCMPLX( AII ), A( I, 1 ), LDA );
               zlacgv(I-1, A( I, 1 ), LDA );
            } else {
               zdscal(I, AII, A( I, 1 ), LDA );
            }
         } // 20
      }

      return;

      // End of ZLAUU2

      }
