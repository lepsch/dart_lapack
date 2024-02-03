      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

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
      int                I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFB, DLARFT, DORG2L, XERBLA
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
      LQUERY = ( LWORK.EQ.-1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 .OR. N.GT.M ) {
         INFO = -2
      } else if ( K.LT.0 .OR. K.GT.N ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      }

      if ( INFO.EQ.0 ) {
         if ( N.EQ.0 ) {
            LWKOPT = 1
         } else {
            NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
            LWKOPT = N*NB
         }
         WORK( 1 ) = LWKOPT

         if ( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) {
            INFO = -8
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('DORGQL', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.LE.0 ) {
         RETURN
      }

      NBMIN = 2
      NX = 0
      IWS = N
      if ( NB.GT.1 .AND. NB.LT.K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
         if ( NX.LT.K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = N
            IWS = LDWORK*NB
            if ( LWORK.LT.IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
            }
         }
      }

      if ( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) {

         // Use blocked code after the first block.
         // The last kk columns are handled by the block method.

         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )

         // Set A(m-kk+1:m,1:n-kk) to zero.

         for (J = 1; J <= N - KK; J++) { // 20
            for (I = M - KK + 1; I <= M; I++) { // 10
               A( I, J ) = ZERO
            } // 10
         } // 20
      } else {
         KK = 0
      }

      // Use unblocked code for the first or only block.

      dorg2l(M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO );

      if ( KK.GT.0 ) {

         // Use blocked code

         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            if ( N-K+I.GT.1 ) {

               // Form the triangular factor of the block reflector
               // H = H(i+ib-1) . . . H(i+1) H(i)

               dlarft('Backward', 'Columnwise', M-K+I+IB-1, IB, A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left

               dlarfb('Left', 'No transpose', 'Backward', 'Columnwise', M-K+I+IB-1, N-K+I-1, IB, A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA, WORK( IB+1 ), LDWORK );
            }

            // Apply H to rows 1:m-k+i+ib-1 of current block

            dorg2l(M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA, TAU( I ), WORK, IINFO );

            // Set rows m-k+i+ib:m of current block to zero

            for (J = N - K + I; J <= N - K + I + IB - 1; J++) { // 40
               for (L = M - K + I + IB; L <= M; L++) { // 30
                  A( L, J ) = ZERO
               } // 30
            } // 40
         } // 50
      }

      WORK( 1 ) = IWS
      RETURN

      // End of DORGQL

      }
