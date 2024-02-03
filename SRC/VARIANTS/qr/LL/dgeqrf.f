      SUBROUTINE DGEQRF ( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, J, K, LWKOPT, NB, NBMIN, NX, LBWORK, NT, LLWORK;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQR2, DLARFB, DLARFT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, MAX, MIN, REAL
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DROUNDUP_LWORK;
      // EXTERNAL ILAENV, DROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      INFO = 0
      NBMIN = 2
      NX = 0
      IWS = N
      K = MIN( M, N )
      NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )

      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
         // Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
      END IF
*
      // Get NT, the size of the very last T, which is the left-over from in-between K-NX and K to K, eg.:
*
             // NB=3     2NB=6       K=10
             // |        |           |
       // 1--2--3--4--5--6--7--8--9--10
                   // |     \________/
                // K-NX=5      NT=4
*
      // So here 4 x 4 is the last T stored in the workspace
*
      NT = K-CEILING(REAL(K-NX)/REAL(NB))*NB

*
      // optimal workspace = space for dlarfb + space for normal T's + space for the last T
*
      LLWORK = MAX (MAX((N-M)*K, (N-M)*NB), MAX(K*NB, NB*NB))
      LLWORK = CEILING(REAL(LLWORK)/REAL(NB))

      IF( K.EQ.0 ) THEN

         LBWORK = 0
         LWKOPT = 1
         WORK( 1 ) = LWKOPT

      ELSE IF ( NT.GT.NB ) THEN

          LBWORK = K-NT
*
          // Optimal workspace for dlarfb = MAX(1,N)*NT
*
          LWKOPT = (LBWORK+LLWORK)*NB
          WORK( 1 ) = DROUNDUP_LWORK(LWKOPT+NT*NT)

      ELSE

          LBWORK = CEILING(REAL(K)/REAL(NB))*NB
          LWKOPT = (LBWORK+LLWORK-NB)*NB
          WORK( 1 ) = DROUNDUP_LWORK(LWKOPT)

      END IF

*
      // Test the input arguments
*
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF ( .NOT.LQUERY ) THEN
         IF( LWORK.LE.0 .OR. ( M.GT.0 .AND. LWORK.LT.MAX( 1, N ) ) ) INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( K.EQ.0 ) THEN
         RETURN
      END IF
*
      IF( NB.GT.1 .AND. NB.LT.K ) THEN

         IF( NX.LT.K ) THEN
*
            // Determine if workspace is large enough for blocked code.
*
            IF ( NT.LE.NB ) THEN
                IWS = (LBWORK+LLWORK-NB)*NB
            ELSE
                IWS = (LBWORK+LLWORK)*NB+NT*NT
            END IF

            IF( LWORK.LT.IWS ) THEN
*
               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.
*
               IF ( NT.LE.NB ) THEN
                    NB = LWORK / (LLWORK+(LBWORK-NB))
               ELSE
                    NB = (LWORK-NT*NT)/(LBWORK+LLWORK)
               END IF
                NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', M, N, -1, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
         // Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
            // Update the current column using old T's
*
            DO 20 J = 1, I - NB, NB
*
               // Apply H' to A(J:M,I:I+IB-1) from the left
*
               CALL DLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, IB, NB, A( J, J ), LDA, WORK(J), LBWORK, A( J, I ), LDA, WORK(LBWORK*NB+NT*NT+1), IB)

20          CONTINUE
*
            // Compute the QR factorization of the current block
            // A(I:M,I:I+IB-1)
*
            CALL DGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK(LBWORK*NB+NT*NT+1), IINFO )

            IF( I+IB.LE.N ) THEN
*
               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK(I), LBWORK )
*
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
      // Use unblocked code to factor the last or only block.
*
      IF( I.LE.K ) THEN

         IF ( I .NE. 1 )   THEN

             DO 30 J = 1, I - NB, NB
*
                 // Apply H' to A(J:M,I:K) from the left
*
                 CALL DLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, K-I+1, NB, A( J, J ), LDA, WORK(J), LBWORK, A( J, I ), LDA, WORK(LBWORK*NB+NT*NT+1), K-I+1)
30           CONTINUE
              CALL DGEQR2( M-I+1, K-I+1, A( I, I ), LDA, TAU( I ), WORK(LBWORK*NB+NT*NT+1),IINFO )

         ELSE
*
         // Use unblocked code to factor the last or only block.
*
         CALL DGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,IINFO )

         END IF
      END IF


*
      // Apply update to the column M+1:N when N > M
*
      IF ( M.LT.N .AND. I.NE.1) THEN
*
          // Form the last triangular factor of the block reflector
          // H = H(i) H(i+1) . . . H(i+ib-1)
*
          IF ( NT .LE. NB ) THEN
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, K-I+1, A( I, I ), LDA, TAU( I ), WORK(I), LBWORK )
          ELSE
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, K-I+1, A( I, I ), LDA, TAU( I ), WORK(LBWORK*NB+1), NT )
          END IF

*
          // Apply H' to A(1:M,M+1:N) from the left
*
          DO 40 J = 1, K-NX, NB

               IB = MIN( K-J+1, NB )
                CALL DLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, N-M, IB, A( J, J ), LDA, WORK(J), LBWORK, A( J, M+1 ), LDA, WORK(LBWORK*NB+NT*NT+1), N-M)

40       CONTINUE

         IF ( NT.LE.NB ) THEN
             CALL DLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, N-M, K-J+1, A( J, J ), LDA, WORK(J), LBWORK, A( J, M+1 ), LDA, WORK(LBWORK*NB+NT*NT+1), N-M)
         ELSE
             CALL DLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, N-M, K-J+1, A( J, J ), LDA, WORK(LBWORK*NB+1), NT, A( J, M+1 ), LDA, WORK(LBWORK*NB+NT*NT+1), N-M)
         END IF

      END IF

      WORK( 1 ) = DROUNDUP_LWORK(IWS)
      RETURN
*
      // End of DGEQRF
*
      END
