      SUBROUTINE CUNMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             SIDE, TRANS, VECT;
      int                INFO, K, LDA, LDC, LWORK, M, N;
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      bool               APPLYQ, LEFT, LQUERY, NOTRAN;
      String             TRANST;
      int                I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW;
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      EXTERNAL           ILAENV, LSAME, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           CUNMLQ, CUNMQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      APPLYQ = LSAME( VECT, 'Q' )
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q or P and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = MAX( 1, N )
      ELSE
         NQ = N
         NW = MAX( 1, M )
      END IF
      IF( .NOT.APPLYQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( K.LT.0 ) THEN
         INFO = -6
      ELSE IF( ( APPLYQ .AND. LDA.LT.MAX( 1, NQ ) ) .OR. ( .NOT.APPLYQ .AND. LDA.LT.MAX( 1, MIN( NQ, K ) ) ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.NW .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( M.GT.0 .AND. N.GT.0 ) THEN
            IF( APPLYQ ) THEN
               IF( LEFT ) THEN
                  NB = ILAENV( 1, 'CUNMQR', SIDE // TRANS, M-1, N, M-1, -1 )
               ELSE
                  NB = ILAENV( 1, 'CUNMQR', SIDE // TRANS, M, N-1, N-1, -1 )
               END IF
            ELSE
               IF( LEFT ) THEN
                  NB = ILAENV( 1, 'CUNMLQ', SIDE // TRANS, M-1, N, M-1, -1 )
               ELSE
                  NB = ILAENV( 1, 'CUNMLQ', SIDE // TRANS, M, N-1, N-1, -1 )
               END IF
            END IF
            LWKOPT = NW*NB
         ELSE
            LWKOPT = 1
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNMBR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      IF( APPLYQ ) THEN
*
*        Apply Q
*
         IF( NQ.GE.K ) THEN
*
*           Q was determined by a call to CGEBRD with nq >= k
*
            CALL CUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN
*
*           Q was determined by a call to CGEBRD with nq < k
*
            IF( LEFT ) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL CUNMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         END IF
      ELSE
*
*        Apply P
*
         IF( NOTRAN ) THEN
            TRANST = 'C'
         ELSE
            TRANST = 'N'
         END IF
         IF( NQ.GT.K ) THEN
*
*           P was determined by a call to CGEBRD with nq > k
*
            CALL CUNMLQ( SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN
*
*           P was determined by a call to CGEBRD with nq <= k
*
            IF( LEFT ) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL CUNMLQ( SIDE, TRANST, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN
*
*     End of CUNMBR
*
      END
