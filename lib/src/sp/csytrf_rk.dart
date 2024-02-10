      void csytrf_rk(final int UPLO, final int N, final Matrix<double> A, final int LDA, final int E, final Array<int> IPIV, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      int                IPIV( * );
      Complex            A( LDA, * ), E( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                I, IINFO, IP, IWS, K, KB, LDWORK, LWKOPT, NB, NBMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASYF_RK, CSYTF2_RK, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

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
         INFO = -8;
      }

      if ( INFO == 0 ) {

         // Determine the block size

         NB = ilaenv( 1, 'CSYTRF_RK', UPLO, N, -1, -1, -1 );
         LWKOPT = max( 1, N*NB );
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('CSYTRF_RK', -INFO );
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
            NBMIN = max( 2, ilaenv( 2, 'CSYTRF_RK', UPLO, N, -1, -1, -1 ) );
         }
      } else {
         IWS = 1;
      }
      if (NB < NBMIN) NB = N;

      if ( UPPER ) {

         // Factorize A as U*D*U**T using the upper triangle of A

         // K is the main loop index, decreasing from N to 1 in steps of
         // KB, where KB is the number of columns factorized by CLASYF_RK;
         // KB is either NB or NB-1, or K for the last block

         K = N;
         } // 10

         // If K < 1, exit from loop

         if (K < 1) GO TO 15;

         if ( K > NB ) {

            // Factorize columns k-kb+1:k of A and use blocked code to
            // update columns 1:k-kb

            clasyf_rk(UPLO, K, NB, KB, A, LDA, E, IPIV, WORK, LDWORK, IINFO );
         } else {

            // Use unblocked code to factorize columns 1:k of A

            csytf2_rk(UPLO, K, A, LDA, E, IPIV, IINFO );
            KB = K;
         }

         // Set INFO on the first occurrence of a zero pivot

         if (INFO == 0 && IINFO > 0) INFO = IINFO;

         // No need to adjust IPIV


         // Apply permutations to the leading panel 1:k-1

         // Read IPIV from the last block factored, i.e.
         // indices  k-kb+1:k and apply row permutations to the
         // last k+1 colunms k+1:N after that block
         // (We can do the simple loop over IPIV with decrement -1,
         // since the ABS value of IPIV( I ) represents the row index
         // of the interchange with row i in both 1x1 and 2x2 pivot cases)

         if ( K < N ) {
            for (I = K; I >= ( K - KB + 1 ); I--) {
               IP = ( IPIV( I ) ).abs();
               if ( IP != I ) {
                  cswap(N-K, A( I, K+1 ), LDA, A( IP, K+1 ), LDA );
               }
            }
         }

         // Decrease K and return to the start of the main loop

         K = K - KB;
         GO TO 10;

         // This label is the exit from main loop over K decreasing
         // from N to 1 in steps of KB

         } // 15

      } else {

         // Factorize A as L*D*L**T using the lower triangle of A

         // K is the main loop index, increasing from 1 to N in steps of
         // KB, where KB is the number of columns factorized by CLASYF_RK;
         // KB is either NB or NB-1, or N-K+1 for the last block

         K = 1;
         } // 20

         // If K > N, exit from loop

         if (K > N) GO TO 35;

         if ( K <= N-NB ) {

            // Factorize columns k:k+kb-1 of A and use blocked code to
            // update columns k+kb:n

            clasyf_rk(UPLO, N-K+1, NB, KB, A( K, K ), LDA, E( K ), IPIV( K ), WORK, LDWORK, IINFO );


         } else {

            // Use unblocked code to factorize columns k:n of A

            csytf2_rk(UPLO, N-K+1, A( K, K ), LDA, E( K ), IPIV( K ), IINFO );
            KB = N - K + 1;

         }

         // Set INFO on the first occurrence of a zero pivot

         if (INFO == 0 && IINFO > 0) INFO = IINFO + K - 1;

         // Adjust IPIV

         for (I = K; I <= K + KB - 1; I++) {
            if ( IPIV( I ) > 0 ) {
               IPIV[I] = IPIV( I ) + K - 1;
            } else {
               IPIV[I] = IPIV( I ) - K + 1;
            }
         }

         // Apply permutations to the leading panel 1:k-1

         // Read IPIV from the last block factored, i.e.
         // indices  k:k+kb-1 and apply row permutations to the
         // first k-1 colunms 1:k-1 before that block
         // (We can do the simple loop over IPIV with increment 1,
         // since the ABS value of IPIV( I ) represents the row index
         // of the interchange with row i in both 1x1 and 2x2 pivot cases)

         if ( K > 1 ) {
            for (I = K; I <= ( K + KB - 1 ); I++) {
               IP = ( IPIV( I ) ).abs();
               if ( IP != I ) {
                  cswap(K-1, A( I, 1 ), LDA, A( IP, 1 ), LDA );
               }
            }
         }

         // Increase K and return to the start of the main loop

         K = K + KB;
         GO TO 20;

         // This label is the exit from main loop over K increasing
         // from 1 to N in steps of KB

         } // 35

      // End Lower

      }

      WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }
