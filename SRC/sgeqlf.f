      SUBROUTINE SGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, K, KI, KK, LDWORK, LWKOPT, MU, NB, NBMIN, NU, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEQL2, SLARFB, SLARFT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      int                ILAENV
      REAL               SROUNDUP_LWORK
      EXTERNAL           ILAENV, SROUNDUP_LWORK
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
*
      IF( INFO.EQ.0 ) THEN
         K = MIN( M, N )
         IF( K.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'SGEQLF', ' ', M, N, -1, -1 )
            LWKOPT = N*NB
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
         IF( .NOT.LQUERY ) THEN
            IF( LWORK.LE.0 .OR. ( M.GT.0 .AND. LWORK.LT.MAX( 1, N ) ) ) INFO = -7
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQLF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( K.EQ.0 ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      NX = 1
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'SGEQLF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'SGEQLF', ' ', M, N, -1, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially.
*        The last kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
         DO 10 I = K - KK + KI + 1, K - KK + 1, -NB
            IB = MIN( K-I+1, NB )
*
*           Compute the QL factorization of the current block
*           A(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
*
            CALL SGEQL2( M-K+I+IB-1, IB, A( 1, N-K+I ), LDA, TAU( I ), WORK, IINFO )
            IF( N-K+I.GT.1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i+ib-1) . . . H(i+1) H(i)
*
               CALL SLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB, A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H**T to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
*
               CALL SLARFB( 'Left', 'Transpose', 'Backward', 'Columnwise', M-K+I+IB-1, N-K+I-1, IB, A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
         MU = M - K + I + NB - 1
         NU = N - K + I + NB - 1
      ELSE
         MU = M
         NU = N
      END IF
*
*     Use unblocked code to factor the last or only block
*
      IF( MU.GT.0 .AND. NU.GT.0 ) CALL SGEQL2( MU, NU, A, LDA, TAU, WORK, IINFO )
*
      WORK( 1 ) = SROUNDUP_LWORK(IWS)
      RETURN
*
*     End of SGEQLF
*
      END
