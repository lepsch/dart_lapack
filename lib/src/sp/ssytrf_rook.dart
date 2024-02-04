      void ssytrf_rook(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double               A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASYF_ROOK, SSYTF2_ROOK, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LWORK < 1 && !LQUERY ) {
         INFO = -7;
      }

      if ( INFO == 0 ) {

         // Determine the block size

         NB = ILAENV( 1, 'SSYTRF_ROOK', UPLO, N, -1, -1, -1 );
         LWKOPT = max( 1, N*NB );
         WORK[1] = SROUNDUP_LWORK( LWKOPT );
      }

      if ( INFO != 0 ) {
         xerbla('SSYTRF_ROOK', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      NBMIN = 2;
      LDWORK = N;
      if ( NB > 1 && NB < N ) {
         IWS = LDWORK*NB;
         if ( LWORK < IWS ) {
            NB = max( LWORK / LDWORK, 1 );
            NBMIN = max( 2, ILAENV( 2, 'SSYTRF_ROOK', UPLO, N, -1, -1, -1 ) );
         }
      } else {
         IWS = 1;
      }
      if (NB < NBMIN) NB = N;

      if ( UPPER ) {

         // Factorize A as U*D*U**T using the upper triangle of A

         // K is the main loop index, decreasing from N to 1 in steps of
         // KB, where KB is the number of columns factorized by SLASYF_ROOK;
         // KB is either NB or NB-1, or K for the last block

         K = N;
         } // 10

         // If K < 1, exit from loop

         if (K < 1) GO TO 40;

         if ( K > NB ) {

            // Factorize columns k-kb+1:k of A and use blocked code to
            // update columns 1:k-kb

            slasyf_rook(UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK, IINFO );
         } else {

            // Use unblocked code to factorize columns 1:k of A

            ssytf2_rook(UPLO, K, A, LDA, IPIV, IINFO );
            KB = K;
         }

         // Set INFO on the first occurrence of a zero pivot

         if (INFO == 0 && IINFO > 0) INFO = IINFO;

         // No need to adjust IPIV

         // Decrease K and return to the start of the main loop

         K = K - KB;
         GO TO 10;

      } else {

         // Factorize A as L*D*L**T using the lower triangle of A

         // K is the main loop index, increasing from 1 to N in steps of
         // KB, where KB is the number of columns factorized by SLASYF_ROOK;
         // KB is either NB or NB-1, or N-K+1 for the last block

         K = 1;
         } // 20

         // If K > N, exit from loop

         if (K > N) GO TO 40;

         if ( K <= N-NB ) {

            // Factorize columns k:k+kb-1 of A and use blocked code to
            // update columns k+kb:n

            slasyf_rook(UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ), WORK, LDWORK, IINFO );
         } else {

            // Use unblocked code to factorize columns k:n of A

            ssytf2_rook(UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO );
            KB = N - K + 1;
         }

         // Set INFO on the first occurrence of a zero pivot

         if (INFO == 0 && IINFO > 0) INFO = IINFO + K - 1;

         // Adjust IPIV

         for (J = K; J <= K + KB - 1; J++) { // 30
            if ( IPIV( J ) > 0 ) {
               IPIV[J] = IPIV( J ) + K - 1;
            } else {
               IPIV[J] = IPIV( J ) - K + 1;
            }
         } // 30

         // Increase K and return to the start of the main loop

         K = K + KB;
         GO TO 20;

      }

      } // 40

      WORK[1] = SROUNDUP_LWORK( LWKOPT );
      return;
      }