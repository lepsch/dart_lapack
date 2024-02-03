      SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )

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
      int                I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQR2, DLARFB, DLARFT, XERBLA
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

      K = MIN( M, N )
      INFO = 0
      NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( .NOT.LQUERY ) THEN
         IF( LWORK.LE.0 .OR. ( M.GT.0 .AND. LWORK.LT.MAX( 1, N ) ) ) INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         IF( K.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            LWKOPT = N*NB
         END IF
         WORK( 1 ) = LWKOPT
         RETURN
      END IF

      // Quick return if possible

      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF

      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN

            // Determine if workspace is large enough for blocked code.

            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', M, N, -1, -1 ) )
            END IF
         END IF
      END IF

      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN

         // Use blocked code initially

         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )

            // Compute the QR factorization of the current block
            // A(i:m,i:i+ib-1)

            CALL DGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO )
            IF( I+IB.LE.N ) THEN

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, LDWORK )

               // Apply H**T to A(i:m,i+ib:n) from the left

               CALL DLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-I+1, N-I-IB+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF

      // Use unblocked code to factor the last or only block.

      IF( I.LE.K ) CALL DGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO )

      WORK( 1 ) = IWS
      RETURN

      // End of DGEQRF

      END
