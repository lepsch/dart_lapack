      void zunglq(final int M, final int N, final int K, final Matrix<double> A, final int LDA, final int TAU, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, K, LDA, LWORK, M, N;
      Complex         A( LDA, * ), TAU( * ), WORK( * );
      // ..

      Complex         ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      bool               LQUERY;
      int                I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARFB, ZLARFT, ZUNGL2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV

      // Test the input arguments

      INFO = 0;
      NB = ilaenv( 1, 'ZUNGLQ', ' ', M, N, K, -1 );
      LWKOPT = max( 1, M )*NB;
      WORK[1] = LWKOPT;
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < M ) {
         INFO = -2;
      } else if ( K < 0 || K > M ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LWORK < max( 1, M ) && !LQUERY ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('ZUNGLQ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M <= 0 ) {
         WORK[1] = 1;
         return;
      }

      NBMIN = 2;
      NX = 0;
      IWS = M;
      if ( NB > 1 && NB < K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = max( 0, ilaenv( 3, 'ZUNGLQ', ' ', M, N, K, -1 ) );
         if ( NX < K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = M;
            IWS = LDWORK*NB;
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK;
               NBMIN = max( 2, ilaenv( 2, 'ZUNGLQ', ' ', M, N, K, -1 ) );
            }
         }
      }

      if ( NB >= NBMIN && NB < K && NX < K ) {

         // Use blocked code after the last block.
         // The first kk rows are handled by the block method.

         KI = ( ( K-NX-1 ) / NB )*NB;
         KK = min( K, KI+NB );

         // Set A(kk+1:m,1:kk) to zero.

         for (J = 1; J <= KK; J++) { // 20
            for (I = KK + 1; I <= M; I++) { // 10
               A[I][J] = ZERO;
            } // 10
         } // 20
      } else {
         KK = 0;
      }

      // Use unblocked code for the last or only block.

      if (KK < M) zungl2( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, TAU( KK+1 ), WORK, IINFO );

      if ( KK > 0 ) {

         // Use blocked code

         for (I = KI + 1; -NB < 0 ? I >= 1 : I <= 1; I += -NB) { // 50
            IB = min( NB, K-I+1 );
            if ( I+IB <= M ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               zlarft('Forward', 'Rowwise', N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H**H to A(i+ib:m,i:n) from the right

               zlarfb('Right', 'Conjugate transpose', 'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I+IB, I ), LDA, WORK( IB+1 ), LDWORK );
            }

            // Apply H**H to columns i:n of current block

            zungl2(IB, N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO );

            // Set columns 1:i-1 of current block to zero

            for (J = 1; J <= I - 1; J++) { // 40
               for (L = I; L <= I + IB - 1; L++) { // 30
                  A[L][J] = ZERO;
               } // 30
            } // 40
         } // 50
      }

      WORK[1] = IWS;
      }
