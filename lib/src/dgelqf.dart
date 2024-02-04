      void dgelqf(M, N, A, LDA, TAU, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGELQ2, DLARFB, DLARFT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      K = min( M, N );
      NB = ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 );
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      } else if ( !LQUERY ) {
         if( LWORK <= 0 || ( N > 0 && LWORK < max( 1, M ) ) ) INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('DGELQF', -INFO );
         return;
      } else if ( LQUERY ) {
         if ( K == 0 ) {
            LWKOPT = 1;
         } else {
            LWKOPT = M*NB;
         }
         WORK[1] = LWKOPT;
         return;
      }

      // Quick return if possible

      if ( K == 0 ) {
         WORK[1] = 1;
         return;
      }

      NBMIN = 2;
      NX = 0;
      IWS = M;
      if ( NB > 1 && NB < K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = max( 0, ILAENV( 3, 'DGELQF', ' ', M, N, -1, -1 ) );
         if ( NX < K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = M;
            IWS = LDWORK*NB;
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK;
               NBMIN = max( 2, ILAENV( 2, 'DGELQF', ' ', M, N, -1, -1 ) );
            }
         }
      }

      if ( NB >= NBMIN && NB < K && NX < K ) {

         // Use blocked code initially

         for (I = 1; NB < 0 ? I >= K - NX : I <= K - NX; I += NB) { // 10
            IB = min( K-I+1, NB );

            // Compute the LQ factorization of the current block
            // A(i:i+ib-1,i:n)

            dgelq2(IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO );
            if ( I+IB <= M ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               dlarft('Forward', 'Rowwise', N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H to A(i+ib:m,i:n) from the right

               dlarfb('Right', 'No transpose', 'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I+IB, I ), LDA, WORK( IB+1 ), LDWORK );
            }
         } // 10
      } else {
         I = 1;
      }

      // Use unblocked code to factor the last or only block.

      if (I <= K) dgelq2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO );

      WORK[1] = IWS;
      return;
      }