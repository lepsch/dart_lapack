      void dpotf2(UPLO, N, A, LDA, INFO ) {

// -- LAPACK computational routine --
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
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J;
      double             AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame, DISNAN;
      //- double             DDOT;
      // EXTERNAL lsame, DDOT, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
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
         xerbla('DPOTF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**T *U.

         for (J = 1; J <= N; J++) { // 10

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 );
            if ( AJJ <= ZERO || DISNAN( AJJ ) ) {
               A[J, J] = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            A[J, J] = AJJ;

            // Compute elements J+1:N of row J.

            if ( J < N ) {
               dgemv('Transpose', J-1, N-J, -ONE, A( 1, J+1 ), LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA );
               dscal(N-J, ONE / AJJ, A( J, J+1 ), LDA );
            }
         } // 10
      } else {

         // Compute the Cholesky factorization A = L*L**T.

         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ), LDA );
            if ( AJJ <= ZERO || DISNAN( AJJ ) ) {
               A[J, J] = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            A[J, J] = AJJ;

            // Compute elements J+1:N of column J.

            if ( J < N ) {
               dgemv('No transpose', N-J, J-1, -ONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 );
               dscal(N-J, ONE / AJJ, A( J+1, J ), 1 );
            }
         } // 20
      }
      GO TO 40;

      } // 30
      INFO = J;

      } // 40
      return;
      }
