      SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFB, DLARFT, DORGL2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      NB = ILAENV( 1, 'DORGLQ', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, M )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK == -1 )
      if ( M < 0 ) {
         INFO = -1
      } else if ( N < M ) {
         INFO = -2
      } else if ( K < 0 || K > M ) {
         INFO = -3
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -5
      } else if ( LWORK < MAX( 1, M ) && .NOT.LQUERY ) {
         INFO = -8
      }
      if ( INFO != 0 ) {
         xerbla('DORGLQ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M.LE.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      NBMIN = 2
      NX = 0
      IWS = M
      if ( NB > 1 && NB < K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
         if ( NX < K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = M
            IWS = LDWORK*NB
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGLQ', ' ', M, N, K, -1 ) )
            }
         }
      }

      if ( NB.GE.NBMIN && NB < K && NX < K ) {

         // Use blocked code after the last block.
         // The first kk rows are handled by the block method.

         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )

         // Set A(kk+1:m,1:kk) to zero.

         for (J = 1; J <= KK; J++) { // 20
            for (I = KK + 1; I <= M; I++) { // 10
               A( I, J ) = ZERO
            } // 10
         } // 20
      } else {
         KK = 0
      }

      // Use unblocked code for the last or only block.

      if (KK < M) CALL DORGL2( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, TAU( KK+1 ), WORK, IINFO );

      if ( KK > 0 ) {

         // Use blocked code

         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            if ( I+IB.LE.M ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               dlarft('Forward', 'Rowwise', N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H**T to A(i+ib:m,i:n) from the right

               dlarfb('Right', 'Transpose', 'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I+IB, I ), LDA, WORK( IB+1 ), LDWORK );
            }

            // Apply H**T to columns i:n of current block

            dorgl2(IB, N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO );

            // Set columns 1:i-1 of current block to zero

            for (J = 1; J <= I - 1; J++) { // 40
               for (L = I; L <= I + IB - 1; L++) { // 30
                  A( L, J ) = ZERO
               } // 30
            } // 40
         } // 50
      }

      WORK( 1 ) = IWS
      RETURN

      // End of DORGLQ

      }
