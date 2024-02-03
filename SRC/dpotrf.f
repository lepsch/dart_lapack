      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
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
      // EXTERNAL DGEMM, DPOTRF2, DSYRK, DTRSM, XERBLA
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
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DPOTRF', -INFO );
         RETURN;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code.

         dpotrf2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code.

         if ( UPPER ) {

            // Compute the Cholesky factorization A = U**T*U.

            DO 10 J = 1, N, NB;

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = MIN( NB, N-J+1 );
               dsyrk('Upper', 'Transpose', JB, J-1, -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA );
               dpotrf2('Upper', JB, A( J, J ), LDA, INFO );
               if (INFO != 0) GO TO 30;
               if ( J+JB <= N ) {

                  // Compute the current block row.

                  dgemm('Transpose', 'No transpose', JB, N-J-JB+1, J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ), LDA, ONE, A( J, J+JB ), LDA );
                  dtrsm('Left', 'Upper', 'Transpose', 'Non-unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA );
               }
            } // 10

         } else {

            // Compute the Cholesky factorization A = L*L**T.

            DO 20 J = 1, N, NB;

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = MIN( NB, N-J+1 );
               dsyrk('Lower', 'No transpose', JB, J-1, -ONE, A( J, 1 ), LDA, ONE, A( J, J ), LDA );
               dpotrf2('Lower', JB, A( J, J ), LDA, INFO );
               if (INFO != 0) GO TO 30;
               if ( J+JB <= N ) {

                  // Compute the current block column.

                  dgemm('No transpose', 'Transpose', N-J-JB+1, JB, J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ), LDA, ONE, A( J+JB, J ), LDA );
                  dtrsm('Right', 'Lower', 'Transpose', 'Non-unit', N-J-JB+1, JB, ONE, A( J, J ), LDA, A( J+JB, J ), LDA );
               }
            } // 20
         }
      }
      GO TO 40;

      } // 30
      INFO = INFO + J - 1;

      } // 40
      RETURN;

      // End of DPOTRF

      }
