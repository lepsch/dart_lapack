      void ctzrzf(M, N, A, LDA, TAU, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LWORK, M, N;
      Complex            A( LDA, * ), TAU( * ), WORK( * );
      // ..

      Complex            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      bool               LQUERY;
      int                I, IB, IWS, KI, KK, LDWORK, LWKMIN, LWKOPT, M1, MU, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CLARZB, CLARZT, CLATRZ
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
      } else if ( N < M ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }

      if ( INFO == 0 ) {
         if ( M == 0 || M == N ) {
            LWKOPT = 1;
            LWKMIN = 1;
         } else {

            // Determine the block size.

            NB = ilaenv( 1, 'CGERQF', ' ', M, N, -1, -1 );
            LWKOPT = M*NB;
            LWKMIN = max( 1, M );
         }
         WORK[1] = SROUNDUP_LWORK(LWKOPT);

         if ( LWORK < LWKMIN && !LQUERY ) {
            INFO = -7;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CTZRZF', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == 0 ) {
         return;
      } else if ( M == N ) {
         for (I = 1; I <= N; I++) { // 10
            TAU[I] = ZERO;
         } // 10
         return;
      }

      NBMIN = 2;
      NX = 1;
      IWS = M;
      if ( NB > 1 && NB < M ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = max( 0, ilaenv( 3, 'CGERQF', ' ', M, N, -1, -1 ) );
         if ( NX < M ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = M;
            IWS = LDWORK*NB;
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK;
               NBMIN = max( 2, ilaenv( 2, 'CGERQF', ' ', M, N, -1, -1 ) );
            }
         }
      }

      if ( NB >= NBMIN && NB < M && NX < M ) {

         // Use blocked code initially.
         // The last kk rows are handled by the block method.

         M1 = min( M+1, N );
         KI = ( ( M-NX-1 ) / NB )*NB;
         KK = min( M, KI+NB );

         for (I = M - KK + KI + 1; -NB < 0 ? I >= M - KK + 1 : I <= M - KK + 1; I += -NB) { // 20
            IB = min( M-I+1, NB );

            // Compute the TZ factorization of the current block
            // A(i:i+ib-1,i:n)

            clatrz(IB, N-I+1, N-M, A( I, I ), LDA, TAU( I ), WORK );
            if ( I > 1 ) {

               // Form the triangular factor of the block reflector
               // H = H(i+ib-1) . . . H(i+1) H(i)

               clarzt('Backward', 'Rowwise', N-M, IB, A( I, M1 ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H to A(1:i-1,i:n) from the right

               clarzb('Right', 'No transpose', 'Backward', 'Rowwise', I-1, N-I+1, IB, N-M, A( I, M1 ), LDA, WORK, LDWORK, A( 1, I ), LDA, WORK( IB+1 ), LDWORK );
            }
         } // 20
         MU = I + NB - 1;
      } else {
         MU = M;
      }

      // Use unblocked code to factor the last or only block

      if (MU > 0) clatrz( MU, N, N-M, A, LDA, TAU, WORK );

      WORK[1] = SROUNDUP_LWORK(LWKOPT);

      }
