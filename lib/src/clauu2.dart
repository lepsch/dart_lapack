      void clauu2(UPLO, N, A, LDA, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      REAL               AII;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- COMPLEX            CDOTC;
      // EXTERNAL LSAME, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLACGV, CSSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
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
         xerbla('CLAUU2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Compute the product U * U**H.

         for (I = 1; I <= N; I++) { // 10
            AII = REAL( A( I, I ) );
            if ( I < N ) {
               A( I, I ) = AII*AII + REAL( CDOTC( N-I, A( I, I+1 ), LDA, A( I, I+1 ), LDA ) );
               clacgv(N-I, A( I, I+1 ), LDA );
               cgemv('No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, CMPLX( AII ), A( 1, I ), 1 );
               clacgv(N-I, A( I, I+1 ), LDA );
            } else {
               csscal(I, AII, A( 1, I ), 1 );
            }
         } // 10

      } else {

         // Compute the product L**H * L.

         for (I = 1; I <= N; I++) { // 20
            AII = REAL( A( I, I ) );
            if ( I < N ) {
               A( I, I ) = AII*AII + REAL( CDOTC( N-I, A( I+1, I ), 1, A( I+1, I ), 1 ) );
               clacgv(I-1, A( I, 1 ), LDA );
               cgemv('Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, CMPLX( AII ), A( I, 1 ), LDA );
               clacgv(I-1, A( I, 1 ), LDA );
            } else {
               csscal(I, AII, A( I, 1 ), LDA );
            }
         } // 20
      }

      return;
      }
