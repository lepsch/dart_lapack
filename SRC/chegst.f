      SUBROUTINE CHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITYPE, LDA, LDB, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      COMPLEX            CONE, HALF
      const              CONE = ( 1.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                K, KB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHEGS2, CHEMM, CHER2K, CTRMM, CTRSM, XERBLA
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
      } else if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('CHEGST', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'CHEGST', UPLO, N, -1, -1, -1 )

      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code

         chegs2(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      } else {

         // Use blocked code

         if ( ITYPE == 1 ) {
            if ( UPPER ) {

               // Compute inv(U**H)*A*inv(U)

               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the upper triangle of A(k:n,k:n)

                  chegs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
                  if ( K+KB <= N ) {
                     ctrsm('Left', UPLO, 'Conjugate transpose', 'Non-unit', KB, N-K-KB+1, CONE, B( K, K ), LDB, A( K, K+KB ), LDA );
                     chemm('Left', UPLO, KB, N-K-KB+1, -HALF, A( K, K ), LDA, B( K, K+KB ), LDB, CONE, A( K, K+KB ), LDA );
                     cher2k(UPLO, 'Conjugate transpose', N-K-KB+1, KB, -CONE, A( K, K+KB ), LDA, B( K, K+KB ), LDB, ONE, A( K+KB, K+KB ), LDA );
                     chemm('Left', UPLO, KB, N-K-KB+1, -HALF, A( K, K ), LDA, B( K, K+KB ), LDB, CONE, A( K, K+KB ), LDA );
                     ctrsm('Right', UPLO, 'No transpose', 'Non-unit', KB, N-K-KB+1, CONE, B( K+KB, K+KB ), LDB, A( K, K+KB ), LDA );
                  }
               } // 10
            } else {

               // Compute inv(L)*A*inv(L**H)

               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the lower triangle of A(k:n,k:n)

                  chegs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
                  if ( K+KB <= N ) {
                     ctrsm('Right', UPLO, 'Conjugate transpose', 'Non-unit', N-K-KB+1, KB, CONE, B( K, K ), LDB, A( K+KB, K ), LDA );
                     chemm('Right', UPLO, N-K-KB+1, KB, -HALF, A( K, K ), LDA, B( K+KB, K ), LDB, CONE, A( K+KB, K ), LDA );
                     cher2k(UPLO, 'No transpose', N-K-KB+1, KB, -CONE, A( K+KB, K ), LDA, B( K+KB, K ), LDB, ONE, A( K+KB, K+KB ), LDA );
                     chemm('Right', UPLO, N-K-KB+1, KB, -HALF, A( K, K ), LDA, B( K+KB, K ), LDB, CONE, A( K+KB, K ), LDA );
                     ctrsm('Left', UPLO, 'No transpose', 'Non-unit', N-K-KB+1, KB, CONE, B( K+KB, K+KB ), LDB, A( K+KB, K ), LDA );
                  }
               } // 20
            }
         } else {
            if ( UPPER ) {

               // Compute U*A*U**H

               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the upper triangle of A(1:k+kb-1,1:k+kb-1)

                  ctrmm('Left', UPLO, 'No transpose', 'Non-unit', K-1, KB, CONE, B, LDB, A( 1, K ), LDA );
                  chemm('Right', UPLO, K-1, KB, HALF, A( K, K ), LDA, B( 1, K ), LDB, CONE, A( 1, K ), LDA );
                  cher2k(UPLO, 'No transpose', K-1, KB, CONE, A( 1, K ), LDA, B( 1, K ), LDB, ONE, A, LDA );
                  chemm('Right', UPLO, K-1, KB, HALF, A( K, K ), LDA, B( 1, K ), LDB, CONE, A( 1, K ), LDA );
                  ctrmm('Right', UPLO, 'Conjugate transpose', 'Non-unit', K-1, KB, CONE, B( K, K ), LDB, A( 1, K ), LDA );
                  chegs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
               } // 30
            } else {

               // Compute L**H*A*L

               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the lower triangle of A(1:k+kb-1,1:k+kb-1)

                  ctrmm('Right', UPLO, 'No transpose', 'Non-unit', KB, K-1, CONE, B, LDB, A( K, 1 ), LDA );
                  chemm('Left', UPLO, KB, K-1, HALF, A( K, K ), LDA, B( K, 1 ), LDB, CONE, A( K, 1 ), LDA );
                  cher2k(UPLO, 'Conjugate transpose', K-1, KB, CONE, A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A, LDA );
                  chemm('Left', UPLO, KB, K-1, HALF, A( K, K ), LDA, B( K, 1 ), LDB, CONE, A( K, 1 ), LDA );
                  ctrmm('Left', UPLO, 'Conjugate transpose', 'Non-unit', KB, K-1, CONE, B( K, K ), LDB, A( K, 1 ), LDA );
                  chegs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
               } // 40
            }
         }
      }
      RETURN

      // End of CHEGST

      }
