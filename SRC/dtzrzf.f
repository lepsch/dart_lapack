      SUBROUTINE DTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      int                I, IB, IWS, KI, KK, LDWORK, LWKMIN, LWKOPT, M1, MU, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DLARZB, DLARZT, DLATRZ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      int                ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( M.EQ.0 .OR. M.EQ.N ) THEN
            LWKOPT = 1
            LWKMIN = 1
         ELSE
*
*           Determine the block size.
*
            NB = ILAENV( 1, 'DGERQF', ' ', M, N, -1, -1 )
            LWKOPT = M*NB
            LWKMIN = MAX( 1, M )
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -7
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTZRZF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 ) THEN
         RETURN
      ELSE IF( M.EQ.N ) THEN
         DO 10 I = 1, N
            TAU( I ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
      NBMIN = 2
      NX = 1
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.M ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DGERQF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.M ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGERQF', ' ', M, N, -1, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.M .AND. NX.LT.M ) THEN
*
*        Use blocked code initially.
*        The last kk rows are handled by the block method.
*
         M1 = MIN( M+1, N )
         KI = ( ( M-NX-1 ) / NB )*NB
         KK = MIN( M, KI+NB )
*
         DO 20 I = M - KK + KI + 1, M - KK + 1, -NB
            IB = MIN( M-I+1, NB )
*
*           Compute the TZ factorization of the current block
*           A(i:i+ib-1,i:n)
*
            CALL DLATRZ( IB, N-I+1, N-M, A( I, I ), LDA, TAU( I ), WORK )
            IF( I.GT.1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i+ib-1) . . . H(i+1) H(i)
*
               CALL DLARZT( 'Backward', 'Rowwise', N-M, IB, A( I, M1 ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(1:i-1,i:n) from the right
*
               CALL DLARZB( 'Right', 'No transpose', 'Backward', 'Rowwise', I-1, N-I+1, IB, N-M, A( I, M1 ), LDA, WORK, LDWORK, A( 1, I ), LDA, WORK( IB+1 ), LDWORK )
            END IF
   20    CONTINUE
         MU = I + NB - 1
      ELSE
         MU = M
      END IF
*
*     Use unblocked code to factor the last or only block
*
      IF( MU.GT.0 ) CALL DLATRZ( MU, N, N-M, A, LDA, TAU, WORK )
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DTZRZF
*
      END
