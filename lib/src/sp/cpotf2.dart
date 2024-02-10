      void cpotf2(UPLO, N, A, LDA, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      Complex            A( LDA, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               UPPER;
      int                J;
      double               AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      //- COMPLEX            CDOTC;
      // EXTERNAL lsame, CDOTC, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLACGV, CSSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL, SQRT

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
         xerbla('CPOTF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**H *U.

         for (J = 1; J <= N; J++) { // 10

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = REAL( double( A( J, J ) ) - CDOTC( J-1, A( 1, J ), 1, A( 1, J ), 1 ) );
            if ( AJJ <= ZERO || SISNAN( AJJ ) ) {
               A[J][J] = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            A[J][J] = AJJ;

            // Compute elements J+1:N of row J.

            if ( J < N ) {
               clacgv(J-1, A( 1, J ), 1 );
               cgemv('Transpose', J-1, N-J, -CONE, A( 1, J+1 ), LDA, A( 1, J ), 1, CONE, A( J, J+1 ), LDA );
               clacgv(J-1, A( 1, J ), 1 );
               csscal(N-J, ONE / AJJ, A( J, J+1 ), LDA );
            }
         } // 10
      } else {

         // Compute the Cholesky factorization A = L*L**H.

         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = REAL( double( A( J, J ) ) - CDOTC( J-1, A( J, 1 ), LDA, A( J, 1 ), LDA ) );
            if ( AJJ <= ZERO || SISNAN( AJJ ) ) {
               A[J][J] = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            A[J][J] = AJJ;

            // Compute elements J+1:N of column J.

            if ( J < N ) {
               clacgv(J-1, A( J, 1 ), LDA );
               cgemv('No transpose', N-J, J-1, -CONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, CONE, A( J+1, J ), 1 );
               clacgv(J-1, A( J, 1 ), LDA );
               csscal(N-J, ONE / AJJ, A( J+1, J ), 1 );
            }
         } // 20
      }
      GO TO 40;

      } // 30
      INFO = J;

      } // 40
      }
