      void sgerqf(final int M, final int N, final Matrix<double> A, final int LDA, final int TAU, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LWORK, M, N;
      double               A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, K, KI, KK, LDWORK, LWKOPT, MU, NB, NBMIN, NU, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGERQ2, SLARFB, SLARFT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SROUNDUP_LWORK

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }

      if ( INFO == 0 ) {
         K = min( M, N );
         if ( K == 0 ) {
            LWKOPT = 1;
         } else {
            NB = ilaenv( 1, 'SGERQF', ' ', M, N, -1, -1 );
            LWKOPT = M*NB;
         }
         WORK[1] = SROUNDUP_LWORK(LWKOPT);

         if ( !LQUERY ) {
            if( LWORK <= 0 || ( N > 0 && LWORK < max( 1, M ) ) ) INFO = -7;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SGERQF', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( K == 0 ) {
         return;
      }

      NBMIN = 2;
      NX = 1;
      IWS = M;
      if ( NB > 1 && NB < K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = max( 0, ilaenv( 3, 'SGERQF', ' ', M, N, -1, -1 ) );
         if ( NX < K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = M;
            IWS = LDWORK*NB;
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK;
               NBMIN = max( 2, ilaenv( 2, 'SGERQF', ' ', M, N, -1, -1 ) );
            }
         }
      }

      if ( NB >= NBMIN && NB < K && NX < K ) {

         // Use blocked code initially.
         // The last kk rows are handled by the block method.

         KI = ( ( K-NX-1 ) / NB )*NB;
         KK = min( K, KI+NB );

         for (I = K - KK + KI + 1; -NB < 0 ? I >= K - KK + 1 : I <= K - KK + 1; I += -NB) { // 10
            IB = min( K-I+1, NB );

            // Compute the RQ factorization of the current block
            // A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1)

            sgerq2(IB, N-K+I+IB-1, A( M-K+I, 1 ), LDA, TAU( I ), WORK, IINFO );
            if ( M-K+I > 1 ) {

               // Form the triangular factor of the block reflector
               // H = H(i+ib-1) . . . H(i+1) H(i)

               slarft('Backward', 'Rowwise', N-K+I+IB-1, IB, A( M-K+I, 1 ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right

               slarfb('Right', 'No transpose', 'Backward', 'Rowwise', M-K+I-1, N-K+I+IB-1, IB, A( M-K+I, 1 ), LDA, WORK, LDWORK, A, LDA, WORK( IB+1 ), LDWORK );
            }
         } // 10
         MU = M - K + I + NB - 1;
         NU = N - K + I + NB - 1;
      } else {
         MU = M;
         NU = N;
      }

      // Use unblocked code to factor the last or only block

      if (MU > 0 && NU > 0) sgerq2( MU, NU, A, LDA, TAU, WORK, IINFO );

      WORK[1] = SROUNDUP_LWORK(IWS);
      }
