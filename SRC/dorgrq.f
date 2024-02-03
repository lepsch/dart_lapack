      SUBROUTINE DORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

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
      PARAMETER          ( ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, II, IINFO, IWS, J, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFB, DLARFT, DORGR2, XERBLA
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
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF

      IF( INFO.EQ.0 ) THEN
         IF( M.LE.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'DORGRQ', ' ', M, N, K, -1 )
            LWKOPT = M*NB
         END IF
         WORK( 1 ) = LWKOPT

         IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGRQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( M.LE.0 ) THEN
         RETURN
      END IF

      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'DORGRQ', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN

            // Determine if workspace is large enough for blocked code.

            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGRQ', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF

      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN

         // Use blocked code after the first block.
         // The last kk rows are handled by the block method.

         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )

         // Set A(1:m-kk,n-kk+1:n) to zero.

         DO 20 J = N - KK + 1, N
            DO 10 I = 1, M - KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF

      // Use unblocked code for the first or only block.

      CALL DORGR2( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )

      IF( KK.GT.0 ) THEN

         // Use blocked code

         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            II = M - K + I
            IF( II.GT.1 ) THEN

               // Form the triangular factor of the block reflector
               // H = H(i+ib-1) . . . H(i+1) H(i)

               CALL DLARFT( 'Backward', 'Rowwise', N-K+I+IB-1, IB, A( II, 1 ), LDA, TAU( I ), WORK, LDWORK )

               // Apply H**T to A(1:m-k+i-1,1:n-k+i+ib-1) from the right

               CALL DLARFB( 'Right', 'Transpose', 'Backward', 'Rowwise', II-1, N-K+I+IB-1, IB, A( II, 1 ), LDA, WORK, LDWORK, A, LDA, WORK( IB+1 ), LDWORK )
            END IF

            // Apply H**T to columns 1:n-k+i+ib-1 of current block

            CALL DORGR2( IB, N-K+I+IB-1, IB, A( II, 1 ), LDA, TAU( I ), WORK, IINFO )

            // Set columns n-k+i+ib:n of current block to zero

            DO 40 L = N - K + I + IB, N
               DO 30 J = II, II + IB - 1
                  A( J, L ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF

      WORK( 1 ) = IWS
      RETURN

      // End of DORGRQ

      END
