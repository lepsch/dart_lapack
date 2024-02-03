      void zpotrf(UPLO, N, A, LDA, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      Complex            A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      Complex         CONE;
      const              ONE = 1.0, CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, JB, NB;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMM, ZHERK, ZPOTRF2, ZTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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
         xerbla('ZPOTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'ZPOTRF', UPLO, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code.

         zpotrf2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code.

         if ( UPPER ) {

            // Compute the Cholesky factorization A = U'*U.

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 10

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = min( NB, N-J+1 );

               zpotrf2('Upper', JB, A( J, J ), LDA, INFO );
                if (INFO != 0) GO TO 30;

               if ( J+JB <= N ) {

                  // Updating the trailing submatrix.

                  ztrsm('Left', 'Upper', 'Conjugate Transpose', 'Non-unit', JB, N-J-JB+1, CONE, A( J, J ), LDA, A( J, J+JB ), LDA );
                  zherk('Upper', 'Conjugate transpose', N-J-JB+1, JB, -ONE, A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA );
               }
            } // 10

         } else {

            // Compute the Cholesky factorization A = L*L'.

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 20

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = min( NB, N-J+1 );

               zpotrf2('Lower', JB, A( J, J ), LDA, INFO );
                if (INFO != 0) GO TO 30;

               if ( J+JB <= N ) {

                 // Updating the trailing submatrix.

                 ztrsm('Right', 'Lower', 'Conjugate Transpose', 'Non-unit', N-J-JB+1, JB, CONE, A( J, J ), LDA, A( J+JB, J ), LDA );
                  zherk('Lower', 'No Transpose', N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, ONE, A( J+JB, J+JB ), LDA );
               }
            } // 20
         }
      }
      GO TO 40;

      } // 30
      INFO = INFO + J - 1;

      } // 40
      return;
      }
