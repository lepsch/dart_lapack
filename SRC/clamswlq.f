      SUBROUTINE CLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T,
     $                     LDT, C, LDC, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), WORK( * ), C( LDC, * ),
     $                   T( LDT, * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN, LQUERY
      INTEGER            I, II, KK, LW, CTR, MINMNK, LWMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SROUNDUP_LWORK
      EXTERNAL           LSAME, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           CTPMLQT, CGEMLQT, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY  = ( LWORK.EQ.-1 )
      NOTRAN  = LSAME( TRANS, 'N' )
      TRAN    = LSAME( TRANS, 'C' )
      LEFT    = LSAME( SIDE, 'L' )
      RIGHT   = LSAME( SIDE, 'R' )
      IF( LEFT ) THEN
        LW = N * MB
      ELSE
        LW = M * MB
      END IF
*
      MINMNK = MIN( M, N, K )
      IF( MINMNK.EQ.0 ) THEN
        LWMIN = 1
      ELSE
        LWMIN = MAX( 1, LW )
      END IF
*
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
        INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
        INFO = -2
      ELSE IF( K.LT.0 ) THEN
        INFO = -5
      ELSE IF( M.LT.K ) THEN
        INFO = -3
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( K.LT.MB .OR. MB.LT.1 ) THEN
        INFO = -6
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
        INFO = -9
      ELSE IF( LDT.LT.MAX( 1, MB ) ) THEN
        INFO = -11
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -13
      ELSE IF( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) THEN
        INFO = -15
      END IF
*
      IF( INFO.EQ.0 ) THEN
        WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'CLAMSWLQ', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( MINMNK.EQ.0 ) THEN
        RETURN
      END IF
*
      IF((NB.LE.K).OR.(NB.GE.MAX(M,N,K))) THEN
        CALL CGEMLQT( SIDE, TRANS, M, N, K, MB, A, LDA,
     $        T, LDT, C, LDC, WORK, INFO )
        RETURN
      END IF
*
      IF(LEFT.AND.TRAN) THEN
*
*         Multiply Q to the last block of C
*
          KK = MOD((M-K),(NB-K))
          CTR = (M-K)/(NB-K)
          IF (KK.GT.0) THEN
            II=M-KK+1
            CALL CTPMLQT('L','C',KK , N, K, 0, MB, A(1,II), LDA,
     $        T(1,CTR*K+1), LDT, C(1,1), LDC,
     $        C(II,1), LDC, WORK, INFO )
          ELSE
            II=M+1
          END IF
*
          DO I=II-(NB-K),NB+1,-(NB-K)
*
*         Multiply Q to the current block of C (1:M,I:I+NB)
*
            CTR = CTR - 1
            CALL CTPMLQT('L','C',NB-K , N, K, 0,MB, A(1,I), LDA,
     $          T(1,CTR*K+1),LDT, C(1,1), LDC,
     $          C(I,1), LDC, WORK, INFO )

          END DO
*
*         Multiply Q to the first block of C (1:M,1:NB)
*
          CALL CGEMLQT('L','C',NB , N, K, MB, A(1,1), LDA, T
     $              ,LDT ,C(1,1), LDC, WORK, INFO )
*
      ELSE IF (LEFT.AND.NOTRAN) THEN
*
*         Multiply Q to the first block of C
*
         KK  = MOD((M-K),(NB-K))
         II  = M-KK+1
         CTR = 1
         CALL CGEMLQT('L','N',NB , N, K, MB, A(1,1), LDA, T
     $              ,LDT ,C(1,1), LDC, WORK, INFO )
*
         DO I=NB+1,II-NB+K,(NB-K)
*
*         Multiply Q to the current block of C (I:I+NB,1:N)
*
          CALL CTPMLQT('L','N',NB-K , N, K, 0,MB, A(1,I), LDA,
     $         T(1, CTR *K+1), LDT, C(1,1), LDC,
     $         C(I,1), LDC, WORK, INFO )
          CTR = CTR + 1
*
         END DO
         IF(II.LE.M) THEN
*
*         Multiply Q to the last block of C
*
          CALL CTPMLQT('L','N',KK , N, K, 0, MB, A(1,II), LDA,
     $        T(1, CTR*K+1), LDT, C(1,1), LDC,
     $        C(II,1), LDC, WORK, INFO )
*
         END IF
*
      ELSE IF(RIGHT.AND.NOTRAN) THEN
*
*         Multiply Q to the last block of C
*
          KK = MOD((N-K),(NB-K))
          CTR = (N-K)/(NB-K)
          IF (KK.GT.0) THEN
            II=N-KK+1
            CALL CTPMLQT('R','N',M , KK, K, 0, MB, A(1, II), LDA,
     $        T(1,CTR*K+1), LDT, C(1,1), LDC,
     $        C(1,II), LDC, WORK, INFO )
          ELSE
            II=N+1
          END IF
*
          DO I=II-(NB-K),NB+1,-(NB-K)
*
*         Multiply Q to the current block of C (1:M,I:I+MB)
*
              CTR = CTR - 1
              CALL CTPMLQT('R','N', M, NB-K, K, 0, MB, A(1, I), LDA,
     $            T(1,CTR*K+1), LDT, C(1,1), LDC,
     $            C(1,I), LDC, WORK, INFO )
          END DO
*
*         Multiply Q to the first block of C (1:M,1:MB)
*
          CALL CGEMLQT('R','N',M , NB, K, MB, A(1,1), LDA, T
     $            ,LDT ,C(1,1), LDC, WORK, INFO )
*
      ELSE IF (RIGHT.AND.TRAN) THEN
*
*       Multiply Q to the first block of C
*
         KK = MOD((N-K),(NB-K))
         II=N-KK+1
         CTR = 1
         CALL CGEMLQT('R','C',M , NB, K, MB, A(1,1), LDA, T
     $            ,LDT ,C(1,1), LDC, WORK, INFO )
*
         DO I=NB+1,II-NB+K,(NB-K)
*
*         Multiply Q to the current block of C (1:M,I:I+MB)
*
          CALL CTPMLQT('R','C',M , NB-K, K, 0,MB, A(1,I), LDA,
     $       T(1,CTR*K+1), LDT, C(1,1), LDC,
     $       C(1,I), LDC, WORK, INFO )
          CTR = CTR + 1
*
         END DO
         IF(II.LE.N) THEN
*
*       Multiply Q to the last block of C
*
          CALL CTPMLQT('R','C',M , KK, K, 0,MB, A(1,II), LDA,
     $      T(1,CTR*K+1),LDT, C(1,1), LDC,
     $      C(1,II), LDC, WORK, INFO )
*
         END IF
*
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      RETURN
*
*     End of CLAMSWLQ
*
      END