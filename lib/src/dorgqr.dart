      void dorgqr(M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFB, DLARFT, DORG2R, XERBLA
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
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 );
      LWKOPT = max( 1, N )*NB;
      WORK[1] = LWKOPT;
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 || N > M ) {
         INFO = -2;
      } else if ( K < 0 || K > N ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LWORK < max( 1, N ) && !LQUERY ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('DORGQR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N <= 0 ) {
         WORK[1] = 1;
         return;
      }

      NBMIN = 2;
      NX = 0;
      IWS = N;
      if ( NB > 1 && NB < K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = max( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) );
         if ( NX < K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = N;
            IWS = LDWORK*NB;
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK;
               NBMIN = max( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) );
            }
         }
      }

      if ( NB >= NBMIN && NB < K && NX < K ) {

         // Use blocked code after the last block.
         // The first kk columns are handled by the block method.

         KI = ( ( K-NX-1 ) / NB )*NB;
         KK = min( K, KI+NB );

         // Set A(1:kk,kk+1:n) to zero.

         for (J = KK + 1; J <= N; J++) { // 20
            for (I = 1; I <= KK; I++) { // 10
               A[I, J] = ZERO;
            } // 10
         } // 20
      } else {
         KK = 0;
      }

      // Use unblocked code for the last or only block.

      if (KK < N) dorg2r( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, TAU( KK+1 ), WORK, IINFO );

      if ( KK > 0 ) {

         // Use blocked code

         for (I = KI + 1; -NB < 0 ? I >= 1 : I <= 1; I += -NB) { // 50
            IB = min( NB, K-I+1 );
            if ( I+IB <= N ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               dlarft('Forward', 'Columnwise', M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H to A(i:m,i+ib:n) from the left

               dlarfb('Left', 'No transpose', 'Forward', 'Columnwise', M-I+1, N-I-IB+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), LDA, WORK( IB+1 ), LDWORK );
            }

            // Apply H to rows i:m of current block

            dorg2r(M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO );

            // Set rows 1:i-1 of current block to zero

            for (J = I; J <= I + IB - 1; J++) { // 40
               for (L = 1; L <= I - 1; L++) { // 30
                  A[L, J] = ZERO;
               } // 30
            } // 40
         } // 50
      }

      WORK[1] = IWS;
      return;
      }
