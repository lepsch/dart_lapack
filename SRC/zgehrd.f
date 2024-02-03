      SUBROUTINE ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1, TSIZE = LDT*NBMAX )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      int                I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB, NBMIN, NH, NX
      COMPLEX*16         EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZGEHD2, ZGEMM, ZLAHR2, ZLARFB, ZTRMM, XERBLA
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
*     Test the input parameters
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      NH = IHI - ILO + 1
      IF( INFO.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         IF( NH.LE.1 ) THEN
            LWKOPT = 1
         ELSE
            NB = MIN( NBMAX, ILAENV( 1, 'ZGEHRD', ' ', N, ILO, IHI, -1 ) )
            LWKOPT = N*NB + TSIZE
         END IF
         WORK( 1 ) = LWKOPT
      ENDIF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEHRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
*
      DO 10 I = 1, ILO - 1
         TAU( I ) = ZERO
   10 CONTINUE
      DO 20 I = MAX( 1, IHI ), N - 1
         TAU( I ) = ZERO
   20 CONTINUE
*
*     Quick return if possible
*
      IF( NH.LE.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Determine the block size
*
      NB = MIN( NBMAX, ILAENV( 1, 'ZGEHRD', ' ', N, ILO, IHI, -1 ) )
      NBMIN = 2
      IF( NB.GT.1 .AND. NB.LT.NH ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code)
*
         NX = MAX( NB, ILAENV( 3, 'ZGEHRD', ' ', N, ILO, IHI, -1 ) )
         IF( NX.LT.NH ) THEN
*
*           Determine if workspace is large enough for blocked code
*
            IF( LWORK.LT.LWKOPT ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code
*
               NBMIN = MAX( 2, ILAENV( 2, 'ZGEHRD', ' ', N, ILO, IHI, -1 ) )
               IF( LWORK.GE.(N*NBMIN + TSIZE) ) THEN
                  NB = (LWORK-TSIZE) / N
               ELSE
                  NB = 1
               END IF
            END IF
         END IF
      END IF
      LDWORK = N
*
      IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
*
*        Use unblocked code below
*
         I = ILO
*
      ELSE
*
*        Use blocked code
*
         IWT = 1 + N*NB
         DO 40 I = ILO, IHI - 1 - NX, NB
            IB = MIN( NB, IHI-I )
*
*           Reduce columns i:i+ib-1 to Hessenberg form, returning the
*           matrices V and T of the block reflector H = I - V*T*V**H
*           which performs the reduction, and also the matrix Y = A*V*T
*
            CALL ZLAHR2( IHI, I, IB, A( 1, I ), LDA, TAU( I ), WORK( IWT ), LDT, WORK, LDWORK )
*
*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
*           right, computing  A := A - Y * V**H. V(i+ib,ib-1) must be set
*           to 1
*
            EI = A( I+IB, I+IB-1 )
            A( I+IB, I+IB-1 ) = ONE
            CALL ZGEMM( 'No transpose', 'Conjugate transpose', IHI, IHI-I-IB+1, IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE, A( 1, I+IB ), LDA )
            A( I+IB, I+IB-1 ) = EI
*
*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
*           right
*
            CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', 'Unit', I, IB-1, ONE, A( I+1, I ), LDA, WORK, LDWORK )
            DO 30 J = 0, IB-2
               CALL ZAXPY( I, -ONE, WORK( LDWORK*J+1 ), 1, A( 1, I+J+1 ), 1 )
   30       CONTINUE
*
*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
*           left
*
            CALL ZLARFB( 'Left', 'Conjugate transpose', 'Forward', 'Columnwise', IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA, WORK( IWT ), LDT, A( I+1, I+IB ), LDA, WORK, LDWORK )
   40    CONTINUE
      END IF
*
*     Use unblocked code to reduce the rest of the matrix
*
      CALL ZGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of ZGEHRD
*
      END
