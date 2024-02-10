      void clauu2(final int UPLO, final int N, final Matrix<double> A, final int LDA, final Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      Complex            A( LDA, * );
      // ..

      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               UPPER;
      int                I;
      double               AII;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- COMPLEX            CDOTC;
      // EXTERNAL lsame, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLACGV, CSSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
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
            AII = double( A( I, I ) );
            if ( I < N ) {
               A[I][I] = AII*AII + double( CDOTC( N-I, A( I, I+1 ), LDA, A( I, I+1 ), LDA ) );
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
            AII = double( A( I, I ) );
            if ( I < N ) {
               A[I][I] = AII*AII + double( CDOTC( N-I, A( I+1, I ), 1, A( I+1, I ), 1 ) );
               clacgv(I-1, A( I, 1 ), LDA );
               cgemv('Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, CMPLX( AII ), A( I, 1 ), LDA );
               clacgv(I-1, A( I, 1 ), LDA );
            } else {
               csscal(I, AII, A( I, 1 ), LDA );
            }
         } // 20
      }

      }
