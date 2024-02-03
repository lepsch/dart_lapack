      SUBROUTINE DORMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      int                INFO, K, L, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1, TSIZE = LDT*NBMAX )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      int                I, I1, I2, I3, IB, IC, IINFO, IWT, JA, JC, LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARZB, DLARZT, DORMR3, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = MAX( 1, N )
      ELSE
         NQ = N
         NW = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( L.LT.0 .OR. ( LEFT .AND. ( L.GT.M ) ) .OR. ( .NOT.LEFT .AND. ( L.GT.N ) ) ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.NW .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         IF( M.EQ.0 .OR. N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = MIN( NBMAX, ILAENV( 1, 'DORMRQ', SIDE // TRANS, M, N, K, -1 ) )
            LWKOPT = NW*NB + TSIZE
         END IF
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMRZ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IF( LWORK.LT.LWKOPT ) THEN
            NB = (LWORK-TSIZE) / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMRQ', SIDE // TRANS, M, N, K, -1 ) )
         END IF
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL DORMR3( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, WORK, IINFO )
      ELSE
*
*        Use blocked code
*
         IWT = 1 + NW*NB
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
            JA = M - L + 1
         ELSE
            MI = M
            IC = 1
            JA = N - L + 1
         END IF
*
         IF( NOTRAN ) THEN
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i+ib-1) . . . H(i+1) H(i)
*
            CALL DLARZT( 'Backward', 'Rowwise', L, IB, A( I, JA ), LDA, TAU( I ), WORK( IWT ), LDT )
*
            IF( LEFT ) THEN
*
*              H or H**T is applied to C(i:m,1:n)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H or H**T is applied to C(1:m,i:n)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           Apply H or H**T
*
            CALL DLARZB( SIDE, TRANST, 'Backward', 'Rowwise', MI, NI, IB, L, A( I, JA ), LDA, WORK( IWT ), LDT, C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
*
      END IF
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DORMRZ
*
      END
