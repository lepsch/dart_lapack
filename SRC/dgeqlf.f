      SUBROUTINE DGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, K, KI, KK, LDWORK, LWKOPT, MU, NB, NBMIN, NU, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQL2, DLARFB, DLARFT, XERBLA
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
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }

      if ( INFO.EQ.0 ) {
         K = MIN( M, N )
         if ( K.EQ.0 ) {
            LWKOPT = 1
         } else {
            NB = ILAENV( 1, 'DGEQLF', ' ', M, N, -1, -1 )
            LWKOPT = N*NB
         }
         WORK( 1 ) = LWKOPT

         if ( .NOT.LQUERY ) {
            IF( LWORK.LE.0 .OR. ( M.GT.0 .AND. LWORK.LT.MAX( 1, N ) ) ) INFO = -7
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('DGEQLF', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( K.EQ.0 ) {
         RETURN
      }

      NBMIN = 2
      NX = 1
      IWS = N
      if ( NB.GT.1 .AND. NB.LT.K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'DGEQLF', ' ', M, N, -1, -1 ) )
         if ( NX.LT.K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = N
            IWS = LDWORK*NB
            if ( LWORK.LT.IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGEQLF', ' ', M, N, -1, -1 ) )
            }
         }
      }

      if ( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) {

         // Use blocked code initially.
         // The last kk columns are handled by the block method.

         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )

         DO 10 I = K - KK + KI + 1, K - KK + 1, -NB
            IB = MIN( K-I+1, NB )

            // Compute the QL factorization of the current block
            // A(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)

            dgeql2(M-K+I+IB-1, IB, A( 1, N-K+I ), LDA, TAU( I ), WORK, IINFO );
            if ( N-K+I.GT.1 ) {

               // Form the triangular factor of the block reflector
               // H = H(i+ib-1) . . . H(i+1) H(i)

               dlarft('Backward', 'Columnwise', M-K+I+IB-1, IB, A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H**T to A(1:m-k+i+ib-1,1:n-k+i-1) from the left

               dlarfb('Left', 'Transpose', 'Backward', 'Columnwise', M-K+I+IB-1, N-K+I-1, IB, A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA, WORK( IB+1 ), LDWORK );
            }
   10    CONTINUE
         MU = M - K + I + NB - 1
         NU = N - K + I + NB - 1
      } else {
         MU = M
         NU = N
      }

      // Use unblocked code to factor the last or only block

      IF( MU.GT.0 .AND. NU.GT.0 ) CALL DGEQL2( MU, NU, A, LDA, TAU, WORK, IINFO )

      WORK( 1 ) = IWS
      RETURN

      // End of DGEQLF

      }
