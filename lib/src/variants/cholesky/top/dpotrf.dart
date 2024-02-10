      void dpotrf(UPLO, N, final Matrix<double> A, final int LDA, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      double             A( LDA, * );
      // ..

      double             ONE;
      const              ONE = 1.0 ;
      bool               UPPER;
      int                J, JB, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DPOTRF2, DSYRK, DTRSM, XERBLA
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
         xerbla('DPOTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Determine the block size for this environment.

      NB = ilaenv( 1, 'DPOTRF', UPLO, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code.

         dpotrf2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code.

         if ( UPPER ) {

            // Compute the Cholesky factorization A = U'*U.

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 10

               JB = min( NB, N-J+1 );

               // Compute the current block.

               dtrsm('Left', 'Upper', 'Transpose', 'Non-unit', J-1, JB, ONE, A( 1, 1 ), LDA, A( 1, J ), LDA );
                dsyrk('Upper', 'Transpose', JB, J-1, -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA );

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               dpotrf2('Upper', JB, A( J, J ), LDA, INFO );
               if (INFO != 0) GO TO 30;

            } // 10

         } else {

            // Compute the Cholesky factorization A = L*L'.

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 20

               JB = min( NB, N-J+1 );

               // Compute the current block.

               dtrsm('Right', 'Lower', 'Transpose', 'Non-unit', JB, J-1, ONE, A( 1, 1 ), LDA, A( J, 1 ), LDA );
                dsyrk('Lower', 'No Transpose', JB, J-1, -ONE, A( J, 1 ), LDA, ONE, A( J, J ), LDA );


               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               dpotrf2('Lower', JB, A( J, J ), LDA, INFO );
               if (INFO != 0) GO TO 30;

            } // 20
         }
      }
      GO TO 40;

      } // 30
      INFO = INFO + J - 1;

      } // 40
      }
