      void zlauum(UPLO, N, A, LDA, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IB, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMM, ZHERK, ZLAUU2, ZTRMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

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
         xerbla('ZLAUUM', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'ZLAUUM', UPLO, N, -1, -1, -1 );

      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code

         zlauu2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code

         if ( UPPER ) {

            // Compute the product U * U**H.

            for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) { // 10
               IB = min( NB, N-I+1 );
               ztrmm('Right', 'Upper', 'Conjugate transpose', 'Non-unit', I-1, IB, CONE, A( I, I ), LDA, A( 1, I ), LDA );
               zlauu2('Upper', IB, A( I, I ), LDA, INFO );
               if ( I+IB <= N ) {
                  zgemm('No transpose', 'Conjugate transpose', I-1, IB, N-I-IB+1, CONE, A( 1, I+IB ), LDA, A( I, I+IB ), LDA, CONE, A( 1, I ), LDA );
                  zherk('Upper', 'No transpose', IB, N-I-IB+1, ONE, A( I, I+IB ), LDA, ONE, A( I, I ), LDA );
               }
            } // 10
         } else {

            // Compute the product L**H * L.

            for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) { // 20
               IB = min( NB, N-I+1 );
               ztrmm('Left', 'Lower', 'Conjugate transpose', 'Non-unit', IB, I-1, CONE, A( I, I ), LDA, A( I, 1 ), LDA );
               zlauu2('Lower', IB, A( I, I ), LDA, INFO );
               if ( I+IB <= N ) {
                  zgemm('Conjugate transpose', 'No transpose', IB, I-1, N-I-IB+1, CONE, A( I+IB, I ), LDA, A( I+IB, 1 ), LDA, CONE, A( I, 1 ), LDA );
                  zherk('Lower', 'Conjugate transpose', IB, N-I-IB+1, ONE, A( I+IB, I ), LDA, ONE, A( I, I ), LDA );
               }
            } // 20
         }
      }

      return;
      }
