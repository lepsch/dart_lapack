      void spotrf(UPLO, N, A, LDA, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, JB, NB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SPOTRF2, SSYRK, STRSM, XERBLA
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
         xerbla('SPOTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'SPOTRF', UPLO, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code.

         spotrf2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code.

         if ( UPPER ) {

            // Compute the Cholesky factorization A = U'*U.

            DO 10 J = 1, N, NB;

               JB = min( NB, N-J+1 );

               // Compute the current block.

               strsm('Left', 'Upper', 'Transpose', 'Non-unit', J-1, JB, ONE, A( 1, 1 ), LDA, A( 1, J ), LDA );
                ssyrk('Upper', 'Transpose', JB, J-1, -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA );

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               spotrf2('Upper', JB, A( J, J ), LDA, INFO );
               if (INFO != 0) GO TO 30;

            } // 10

         } else {

            // Compute the Cholesky factorization A = L*L'.

            DO 20 J = 1, N, NB;

               JB = min( NB, N-J+1 );

               // Compute the current block.

               strsm('Right', 'Lower', 'Transpose', 'Non-unit', JB, J-1, ONE, A( 1, 1 ), LDA, A( J, 1 ), LDA );
                ssyrk('Lower', 'No Transpose', JB, J-1, -ONE, A( J, 1 ), LDA, ONE, A( J, J ), LDA );

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               spotrf2('Lower', JB, A( J, J ), LDA, INFO );
               if (INFO != 0) GO TO 30;

            } // 20
         }
      }
      GO TO 40;

      } // 30
      INFO = INFO + J - 1;

      } // 40
      return;

      // End of SPOTRF

      }
