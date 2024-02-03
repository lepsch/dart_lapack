      void sorgrq(M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, II, IINFO, IWS, J, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARFB, SLARFT, SORGR2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < M ) {
         INFO = -2;
      } else if ( K < 0 || K > M ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      }

      if ( INFO == 0 ) {
         if ( M <= 0 ) {
            LWKOPT = 1;
         } else {
            NB = ILAENV( 1, 'SORGRQ', ' ', M, N, K, -1 );
            LWKOPT = M*NB;
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT);

         if ( LWORK < max( 1, M ) && !LQUERY ) {
            INFO = -8;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SORGRQ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M <= 0 ) {
         return;
      }

      NBMIN = 2;
      NX = 0;
      IWS = M;
      if ( NB > 1 && NB < K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = max( 0, ILAENV( 3, 'SORGRQ', ' ', M, N, K, -1 ) );
         if ( NX < K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = M;
            IWS = LDWORK*NB;
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK;
               NBMIN = max( 2, ILAENV( 2, 'SORGRQ', ' ', M, N, K, -1 ) );
            }
         }
      }

      if ( NB >= NBMIN && NB < K && NX < K ) {

         // Use blocked code after the first block.
         // The last kk rows are handled by the block method.

         KK = min( K, ( ( K-NX+NB-1 ) / NB )*NB );

         // Set A(1:m-kk,n-kk+1:n) to zero.

         for (J = N - KK + 1; J <= N; J++) { // 20
            for (I = 1; I <= M - KK; I++) { // 10
               A( I, J ) = ZERO;
            } // 10
         } // 20
      } else {
         KK = 0;
      }

      // Use unblocked code for the first or only block.

      sorgr2(M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO );

      if ( KK > 0 ) {

         // Use blocked code

         DO 50 I = K - KK + 1, K, NB;
            IB = min( NB, K-I+1 );
            II = M - K + I;
            if ( II > 1 ) {

               // Form the triangular factor of the block reflector
               // H = H(i+ib-1) . . . H(i+1) H(i)

               slarft('Backward', 'Rowwise', N-K+I+IB-1, IB, A( II, 1 ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H**T to A(1:m-k+i-1,1:n-k+i+ib-1) from the right

               slarfb('Right', 'Transpose', 'Backward', 'Rowwise', II-1, N-K+I+IB-1, IB, A( II, 1 ), LDA, WORK, LDWORK, A, LDA, WORK( IB+1 ), LDWORK );
            }

            // Apply H**T to columns 1:n-k+i+ib-1 of current block

            sorgr2(IB, N-K+I+IB-1, IB, A( II, 1 ), LDA, TAU( I ), WORK, IINFO );

            // Set columns n-k+i+ib:n of current block to zero

            for (L = N - K + I + IB; L <= N; L++) { // 40
               for (J = II; J <= II + IB - 1; J++) { // 30
                  A( J, L ) = ZERO;
               } // 30
            } // 40
         } // 50
      }

      WORK( 1 ) = SROUNDUP_LWORK(IWS);
      return;
      }
