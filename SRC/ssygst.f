      SUBROUTINE SSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITYPE, LDA, LDB, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, HALF
      const              ONE = 1.0, HALF = 0.5 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                K, KB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSYGS2, SSYMM, SSYR2K, STRMM, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1
      } else if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('SSYGST', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'SSYGST', UPLO, N, -1, -1, -1 )

      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code

         ssygs2(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      } else {

         // Use blocked code

         if ( ITYPE == 1 ) {
            if ( UPPER ) {

               // Compute inv(U**T)*A*inv(U)

               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the upper triangle of A(k:n,k:n)

                  ssygs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
                  if ( K+KB <= N ) {
                     strsm('Left', UPLO, 'Transpose', 'Non-unit', KB, N-K-KB+1, ONE, B( K, K ), LDB, A( K, K+KB ), LDA );
                     ssymm('Left', UPLO, KB, N-K-KB+1, -HALF, A( K, K ), LDA, B( K, K+KB ), LDB, ONE, A( K, K+KB ), LDA );
                     ssyr2k(UPLO, 'Transpose', N-K-KB+1, KB, -ONE, A( K, K+KB ), LDA, B( K, K+KB ), LDB, ONE, A( K+KB, K+KB ), LDA );
                     ssymm('Left', UPLO, KB, N-K-KB+1, -HALF, A( K, K ), LDA, B( K, K+KB ), LDB, ONE, A( K, K+KB ), LDA );
                     strsm('Right', UPLO, 'No transpose', 'Non-unit', KB, N-K-KB+1, ONE, B( K+KB, K+KB ), LDB, A( K, K+KB ), LDA );
                  }
               } // 10
            } else {

               // Compute inv(L)*A*inv(L**T)

               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the lower triangle of A(k:n,k:n)

                  ssygs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
                  if ( K+KB <= N ) {
                     strsm('Right', UPLO, 'Transpose', 'Non-unit', N-K-KB+1, KB, ONE, B( K, K ), LDB, A( K+KB, K ), LDA );
                     ssymm('Right', UPLO, N-K-KB+1, KB, -HALF, A( K, K ), LDA, B( K+KB, K ), LDB, ONE, A( K+KB, K ), LDA );
                     ssyr2k(UPLO, 'No transpose', N-K-KB+1, KB, -ONE, A( K+KB, K ), LDA, B( K+KB, K ), LDB, ONE, A( K+KB, K+KB ), LDA );
                     ssymm('Right', UPLO, N-K-KB+1, KB, -HALF, A( K, K ), LDA, B( K+KB, K ), LDB, ONE, A( K+KB, K ), LDA );
                     strsm('Left', UPLO, 'No transpose', 'Non-unit', N-K-KB+1, KB, ONE, B( K+KB, K+KB ), LDB, A( K+KB, K ), LDA );
                  }
               } // 20
            }
         } else {
            if ( UPPER ) {

               // Compute U*A*U**T

               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the upper triangle of A(1:k+kb-1,1:k+kb-1)

                  strmm('Left', UPLO, 'No transpose', 'Non-unit', K-1, KB, ONE, B, LDB, A( 1, K ), LDA );
                  ssymm('Right', UPLO, K-1, KB, HALF, A( K, K ), LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA );
                  ssyr2k(UPLO, 'No transpose', K-1, KB, ONE, A( 1, K ), LDA, B( 1, K ), LDB, ONE, A, LDA );
                  ssymm('Right', UPLO, K-1, KB, HALF, A( K, K ), LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA );
                  strmm('Right', UPLO, 'Transpose', 'Non-unit', K-1, KB, ONE, B( K, K ), LDB, A( 1, K ), LDA );
                  ssygs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
               } // 30
            } else {

               // Compute L**T*A*L

               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the lower triangle of A(1:k+kb-1,1:k+kb-1)

                  strmm('Right', UPLO, 'No transpose', 'Non-unit', KB, K-1, ONE, B, LDB, A( K, 1 ), LDA );
                  ssymm('Left', UPLO, KB, K-1, HALF, A( K, K ), LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA );
                  ssyr2k(UPLO, 'Transpose', K-1, KB, ONE, A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A, LDA );
                  ssymm('Left', UPLO, KB, K-1, HALF, A( K, K ), LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA );
                  strmm('Left', UPLO, 'Transpose', 'Non-unit', KB, K-1, ONE, B( K, K ), LDB, A( K, 1 ), LDA );
                  ssygs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
               } // 40
            }
         }
      }
      RETURN

      // End of SSYGST

      }
