      void spotf2(final int UPLO, final int N, final Matrix<double> A, final int LDA, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      double               A( LDA, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      int                J;
      double               AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      //- REAL               SDOT;
      // EXTERNAL lsame, SDOT, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

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
         xerbla('SPOTF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**T *U.

         for (J = 1; J <= N; J++) { // 10

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = A( J, J ) - SDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 );
            if ( AJJ <= ZERO || SISNAN( AJJ ) ) {
               A[J][J] = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            A[J][J] = AJJ;

            // Compute elements J+1:N of row J.

            if ( J < N ) {
               sgemv('Transpose', J-1, N-J, -ONE, A( 1, J+1 ), LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA );
               sscal(N-J, ONE / AJJ, A( J, J+1 ), LDA );
            }
         } // 10
      } else {

         // Compute the Cholesky factorization A = L*L**T.

         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = A( J, J ) - SDOT( J-1, A( J, 1 ), LDA, A( J, 1 ), LDA );
            if ( AJJ <= ZERO || SISNAN( AJJ ) ) {
               A[J][J] = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            A[J][J] = AJJ;

            // Compute elements J+1:N of column J.

            if ( J < N ) {
               sgemv('No transpose', N-J, J-1, -ONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 );
               sscal(N-J, ONE / AJJ, A( J+1, J ), 1 );
            }
         } // 20
      }
      GO TO 40;

      } // 30
      INFO = J;

      } // 40
      }
