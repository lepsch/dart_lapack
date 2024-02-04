      void zpotf2(UPLO, N, A, LDA, INFO ) {

// -- LAPACK computational routine --
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
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J;
      double             AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame, DISNAN;
      //- Complex         ZDOTC;
      // EXTERNAL lsame, ZDOTC, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZGEMV, ZLACGV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, SQRT
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
         xerbla('ZPOTF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**H *U.

         for (J = 1; J <= N; J++) { // 10

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = (A( J, J )).toDouble() - DBLE( ZDOTC( J-1, A( 1, J ), 1, A( 1, J ), 1 ) );
            if ( AJJ <= ZERO || DISNAN( AJJ ) ) {
               A[J, J] = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            A[J, J] = AJJ;

            // Compute elements J+1:N of row J.

            if ( J < N ) {
               zlacgv(J-1, A( 1, J ), 1 );
               zgemv('Transpose', J-1, N-J, -CONE, A( 1, J+1 ), LDA, A( 1, J ), 1, CONE, A( J, J+1 ), LDA );
               zlacgv(J-1, A( 1, J ), 1 );
               zdscal(N-J, ONE / AJJ, A( J, J+1 ), LDA );
            }
         } // 10
      } else {

         // Compute the Cholesky factorization A = L*L**H.

         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = (A( J, J )).toDouble() - DBLE( ZDOTC( J-1, A( J, 1 ), LDA, A( J, 1 ), LDA ) );
            if ( AJJ <= ZERO || DISNAN( AJJ ) ) {
               A[J, J] = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            A[J, J] = AJJ;

            // Compute elements J+1:N of column J.

            if ( J < N ) {
               zlacgv(J-1, A( J, 1 ), LDA );
               zgemv('No transpose', N-J, J-1, -CONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, CONE, A( J+1, J ), 1 );
               zlacgv(J-1, A( J, 1 ), LDA );
               zdscal(N-J, ONE / AJJ, A( J+1, J ), 1 );
            }
         } // 20
      }
      GO TO 40;

      } // 30
      INFO = J;

      } // 40
      return;
      }