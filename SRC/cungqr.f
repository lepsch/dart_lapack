      SUBROUTINE CUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARFB, CLARFT, CUNG2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      NB = ILAENV( 1, 'CUNGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      LQUERY = ( LWORK.EQ.-1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 .OR. N.GT.M ) {
         INFO = -2
      } else if ( K.LT.0 .OR. K.GT.N ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('CUNGQR', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.LE.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      NBMIN = 2
      NX = 0
      IWS = N
      if ( NB.GT.1 .AND. NB.LT.K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'CUNGQR', ' ', M, N, K, -1 ) )
         if ( NX.LT.K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = N
            IWS = LDWORK*NB
            if ( LWORK.LT.IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'CUNGQR', ' ', M, N, K, -1 ) )
            }
         }
      }

      if ( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) {

         // Use blocked code after the last block.
         // The first kk columns are handled by the block method.

         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )

         // Set A(1:kk,kk+1:n) to zero.

         for (J = KK + 1; J <= N; J++) { // 20
            for (I = 1; I <= KK; I++) { // 10
               A( I, J ) = ZERO
            } // 10
         } // 20
      } else {
         KK = 0
      }

      // Use unblocked code for the last or only block.

      IF( KK.LT.N ) CALL CUNG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, TAU( KK+1 ), WORK, IINFO )

      if ( KK.GT.0 ) {

         // Use blocked code

         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            if ( I+IB.LE.N ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               clarft('Forward', 'Columnwise', M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H to A(i:m,i+ib:n) from the left

               clarfb('Left', 'No transpose', 'Forward', 'Columnwise', M-I+1, N-I-IB+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), LDA, WORK( IB+1 ), LDWORK );
            }

            // Apply H to rows i:m of current block

            cung2r(M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO );

            // Set rows 1:i-1 of current block to zero

            for (J = I; J <= I + IB - 1; J++) { // 40
               for (L = 1; L <= I - 1; L++) { // 30
                  A( L, J ) = ZERO
               } // 30
            } // 40
         } // 50
      }

      WORK( 1 ) = SROUNDUP_LWORK(IWS)
      RETURN

      // End of CUNGQR

      }
