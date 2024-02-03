      SUBROUTINE ZHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITYPE, LDA, LDB, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
      COMPLEX*16         CONE, HALF
      const              CONE = ( 1.0D+0, 0.0D+0 ), HALF = ( 0.5D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                K, KB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZHEGS2, ZHEMM, ZHER2K, ZTRMM, ZTRSM
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
      if ( ITYPE.LT.1 .OR. ITYPE.GT.3 ) {
         INFO = -1
      } else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO.NE.0 ) {
         xerbla('ZHEGST', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'ZHEGST', UPLO, N, -1, -1, -1 )

      if ( NB.LE.1 .OR. NB.GE.N ) {

         // Use unblocked code

         zhegs2(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      } else {

         // Use blocked code

         if ( ITYPE == 1 ) {
            if ( UPPER ) {

               // Compute inv(U**H)*A*inv(U)

               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the upper triangle of A(k:n,k:n)

                  zhegs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
                  if ( K+KB.LE.N ) {
                     ztrsm('Left', UPLO, 'Conjugate transpose', 'Non-unit', KB, N-K-KB+1, CONE, B( K, K ), LDB, A( K, K+KB ), LDA );
                     zhemm('Left', UPLO, KB, N-K-KB+1, -HALF, A( K, K ), LDA, B( K, K+KB ), LDB, CONE, A( K, K+KB ), LDA );
                     zher2k(UPLO, 'Conjugate transpose', N-K-KB+1, KB, -CONE, A( K, K+KB ), LDA, B( K, K+KB ), LDB, ONE, A( K+KB, K+KB ), LDA );
                     zhemm('Left', UPLO, KB, N-K-KB+1, -HALF, A( K, K ), LDA, B( K, K+KB ), LDB, CONE, A( K, K+KB ), LDA );
                     ztrsm('Right', UPLO, 'No transpose', 'Non-unit', KB, N-K-KB+1, CONE, B( K+KB, K+KB ), LDB, A( K, K+KB ), LDA );
                  }
               } // 10
            } else {

               // Compute inv(L)*A*inv(L**H)

               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the lower triangle of A(k:n,k:n)

                  zhegs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
                  if ( K+KB.LE.N ) {
                     ztrsm('Right', UPLO, 'Conjugate transpose', 'Non-unit', N-K-KB+1, KB, CONE, B( K, K ), LDB, A( K+KB, K ), LDA );
                     zhemm('Right', UPLO, N-K-KB+1, KB, -HALF, A( K, K ), LDA, B( K+KB, K ), LDB, CONE, A( K+KB, K ), LDA );
                     zher2k(UPLO, 'No transpose', N-K-KB+1, KB, -CONE, A( K+KB, K ), LDA, B( K+KB, K ), LDB, ONE, A( K+KB, K+KB ), LDA );
                     zhemm('Right', UPLO, N-K-KB+1, KB, -HALF, A( K, K ), LDA, B( K+KB, K ), LDB, CONE, A( K+KB, K ), LDA );
                     ztrsm('Left', UPLO, 'No transpose', 'Non-unit', N-K-KB+1, KB, CONE, B( K+KB, K+KB ), LDB, A( K+KB, K ), LDA );
                  }
               } // 20
            }
         } else {
            if ( UPPER ) {

               // Compute U*A*U**H

               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the upper triangle of A(1:k+kb-1,1:k+kb-1)

                  ztrmm('Left', UPLO, 'No transpose', 'Non-unit', K-1, KB, CONE, B, LDB, A( 1, K ), LDA );
                  zhemm('Right', UPLO, K-1, KB, HALF, A( K, K ), LDA, B( 1, K ), LDB, CONE, A( 1, K ), LDA );
                  zher2k(UPLO, 'No transpose', K-1, KB, CONE, A( 1, K ), LDA, B( 1, K ), LDB, ONE, A, LDA );
                  zhemm('Right', UPLO, K-1, KB, HALF, A( K, K ), LDA, B( 1, K ), LDB, CONE, A( 1, K ), LDA );
                  ztrmm('Right', UPLO, 'Conjugate transpose', 'Non-unit', K-1, KB, CONE, B( K, K ), LDB, A( 1, K ), LDA );
                  zhegs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
               } // 30
            } else {

               // Compute L**H*A*L

               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )

                  // Update the lower triangle of A(1:k+kb-1,1:k+kb-1)

                  ztrmm('Right', UPLO, 'No transpose', 'Non-unit', KB, K-1, CONE, B, LDB, A( K, 1 ), LDA );
                  zhemm('Left', UPLO, KB, K-1, HALF, A( K, K ), LDA, B( K, 1 ), LDB, CONE, A( K, 1 ), LDA );
                  zher2k(UPLO, 'Conjugate transpose', K-1, KB, CONE, A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A, LDA );
                  zhemm('Left', UPLO, KB, K-1, HALF, A( K, K ), LDA, B( K, 1 ), LDB, CONE, A( K, 1 ), LDA );
                  ztrmm('Left', UPLO, 'Conjugate transpose', 'Non-unit', KB, K-1, CONE, B( K, K ), LDB, A( K, 1 ), LDA );
                  zhegs2(ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO );
               } // 40
            }
         }
      }
      RETURN

      // End of ZHEGST

      }
