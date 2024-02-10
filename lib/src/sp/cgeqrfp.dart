      void cgeqrfp(M, N, final Matrix<double> A, final int LDA, TAU, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LWORK, M, N;
      Complex            A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, K, LDWORK, LWKMIN, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQR2P, CLARFB, CLARFT, XERBLA
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
      NB = ilaenv( 1, 'CGEQRF', ' ', M, N, -1, -1 );
      K = min( M, N );
      if ( K == 0 ) {
         LWKMIN = 1;
         LWKOPT = 1;
      } else {
         LWKMIN = N;
         LWKOPT = N*NB;
      }
      WORK[1] = SROUNDUP_LWORK( LWKOPT );

      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('CGEQRFP', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( K == 0 ) {
         WORK[1] = 1;
         return;
      }

      NBMIN = 2;
      NX = 0;
      IWS = LWKMIN;
      if ( NB > 1 && NB < K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = max( 0, ilaenv( 3, 'CGEQRF', ' ', M, N, -1, -1 ) );
         if ( NX < K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = N;
            IWS = LDWORK*NB;
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK;
               NBMIN = max( 2, ilaenv( 2, 'CGEQRF', ' ', M, N, -1, -1 ) );
            }
         }
      }

      if ( NB >= NBMIN && NB < K && NX < K ) {

         // Use blocked code initially

         for (I = 1; NB < 0 ? I >= K - NX : I <= K - NX; I += NB) { // 10
            IB = min( K-I+1, NB );

            // Compute the QR factorization of the current block
            // A(i:m,i:i+ib-1)

            cgeqr2p(M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO );
            if ( I+IB <= N ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               clarft('Forward', 'Columnwise', M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H**H to A(i:m,i+ib:n) from the left

               clarfb('Left', 'Conjugate transpose', 'Forward', 'Columnwise', M-I+1, N-I-IB+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), LDA, WORK( IB+1 ), LDWORK );
            }
         } // 10
      } else {
         I = 1;
      }

      // Use unblocked code to factor the last or only block.

      if (I <= K) cgeqr2p( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO );

      WORK[1] = SROUNDUP_LWORK( IWS );
      }
