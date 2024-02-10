      void cpotrf(final int UPLO, final int N, final Matrix<double> A, final int LDA, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      Complex            A( LDA, * );
      // ..

      double               ONE;
      Complex            CONE;
      const              ONE = 1.0, CONE = ( 1.0, 0.0 ) ;
      bool               UPPER;
      int                J, JB, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHERK, CPOTRF2, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

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
         xerbla('CPOTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Determine the block size for this environment.

      NB = ilaenv( 1, 'CPOTRF', UPLO, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code.

         cpotrf2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code.

         if ( UPPER ) {

            // Compute the Cholesky factorization A = U**H *U.

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 10

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = min( NB, N-J+1 );
               cherk('Upper', 'Conjugate transpose', JB, J-1, -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA );
               cpotrf2('Upper', JB, A( J, J ), LDA, INFO );
               if (INFO != 0) GO TO 30;
               if ( J+JB <= N ) {

                  // Compute the current block row.

                  cgemm('Conjugate transpose', 'No transpose', JB, N-J-JB+1, J-1, -CONE, A( 1, J ), LDA, A( 1, J+JB ), LDA, CONE, A( J, J+JB ), LDA );
                  ctrsm('Left', 'Upper', 'Conjugate transpose', 'Non-unit', JB, N-J-JB+1, CONE, A( J, J ), LDA, A( J, J+JB ), LDA );
               }
            } // 10

         } else {

            // Compute the Cholesky factorization A = L*L**H.

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 20

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = min( NB, N-J+1 );
               cherk('Lower', 'No transpose', JB, J-1, -ONE, A( J, 1 ), LDA, ONE, A( J, J ), LDA );
               cpotrf2('Lower', JB, A( J, J ), LDA, INFO );
               if (INFO != 0) GO TO 30;
               if ( J+JB <= N ) {

                  // Compute the current block column.

                  cgemm('No transpose', 'Conjugate transpose', N-J-JB+1, JB, J-1, -CONE, A( J+JB, 1 ), LDA, A( J, 1 ), LDA, CONE, A( J+JB, J ), LDA );
                  ctrsm('Right', 'Lower', 'Conjugate transpose', 'Non-unit', N-J-JB+1, JB, CONE, A( J, J ), LDA, A( J+JB, J ), LDA );
               }
            } // 20
         }
      }
      GO TO 40;

      } // 30
      INFO = INFO + J - 1;

      } // 40
      }
