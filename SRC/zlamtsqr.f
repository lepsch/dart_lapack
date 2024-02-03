      SUBROUTINE ZLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, LDT, C, LDC, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      int                INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), WORK( * ), C( LDC, * ), T( LDT, * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN, LQUERY
      int                I, II, KK, LW, CTR, Q, MINMNK, LWMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMQRT, ZTPMQRT, XERBLA
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
        LW = N * NB
        Q = M
      ELSE
        LW = M * NB
        Q = N
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
      ELSE IF( M.LT.K ) THEN
        INFO = -3
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( K.LT.0 ) THEN
        INFO = -5
      ELSE IF( K.LT.NB .OR. NB.LT.1 ) THEN
        INFO = -7
      ELSE IF( LDA.LT.MAX( 1, Q ) ) THEN
        INFO = -9
      ELSE IF( LDT.LT.MAX( 1, NB ) ) THEN
        INFO = -11
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -13
      ELSE IF( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) THEN
        INFO = -15
      END IF
*
      IF( INFO.EQ.0 )  THEN
        WORK( 1 ) = LWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'ZLAMTSQR', -INFO )
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
*     Determine the block size if it is tall skinny or short and wide
*
      IF((MB.LE.K).OR.(MB.GE.MAX(M,N,K))) THEN
        CALL ZGEMQRT( SIDE, TRANS, M, N, K, NB, A, LDA, T, LDT, C, LDC, WORK, INFO )
        RETURN
      END IF
*
      IF(LEFT.AND.NOTRAN) THEN
*
*         Multiply Q to the last block of C
*
         KK = MOD((M-K),(MB-K))
         CTR = (M-K)/(MB-K)
         IF (KK.GT.0) THEN
           II=M-KK+1
           CALL ZTPMQRT('L','N',KK , N, K, 0, NB, A(II,1), LDA, T(1, CTR * K + 1),LDT , C(1,1), LDC, C(II,1), LDC, WORK, INFO )
         ELSE
           II=M+1
         END IF
*
         DO I=II-(MB-K),MB+1,-(MB-K)
*
*         Multiply Q to the current block of C (I:I+MB,1:N)
*
           CTR = CTR - 1
           CALL ZTPMQRT('L','N',MB-K , N, K, 0,NB, A(I,1), LDA, T(1,CTR * K + 1),LDT, C(1,1), LDC, C(I,1), LDC, WORK, INFO )

         END DO
*
*         Multiply Q to the first block of C (1:MB,1:N)
*
         CALL ZGEMQRT('L','N',MB , N, K, NB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO )
*
      ELSE IF (LEFT.AND.TRAN) THEN
*
*         Multiply Q to the first block of C
*
         KK = MOD((M-K),(MB-K))
         II=M-KK+1
         CTR = 1
         CALL ZGEMQRT('L','C',MB , N, K, NB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO )
*
         DO I=MB+1,II-MB+K,(MB-K)
*
*         Multiply Q to the current block of C (I:I+MB,1:N)
*
          CALL ZTPMQRT('L','C',MB-K , N, K, 0,NB, A(I,1), LDA, T(1,CTR * K + 1),LDT, C(1,1), LDC, C(I,1), LDC, WORK, INFO )
          CTR = CTR + 1
*
         END DO
         IF(II.LE.M) THEN
*
*         Multiply Q to the last block of C
*
          CALL ZTPMQRT('L','C',KK , N, K, 0,NB, A(II,1), LDA, T(1, CTR * K + 1), LDT, C(1,1), LDC, C(II,1), LDC, WORK, INFO )
*
         END IF
*
      ELSE IF(RIGHT.AND.TRAN) THEN
*
*         Multiply Q to the last block of C
*
          KK = MOD((N-K),(MB-K))
          CTR = (N-K)/(MB-K)
          IF (KK.GT.0) THEN
            II=N-KK+1
            CALL ZTPMQRT('R','C',M , KK, K, 0, NB, A(II,1), LDA, T(1,CTR * K + 1), LDT, C(1,1), LDC, C(1,II), LDC, WORK, INFO )
          ELSE
            II=N+1
          END IF
*
          DO I=II-(MB-K),MB+1,-(MB-K)
*
*         Multiply Q to the current block of C (1:M,I:I+MB)
*
            CTR = CTR - 1
            CALL ZTPMQRT('R','C',M , MB-K, K, 0,NB, A(I,1), LDA, T(1, CTR * K + 1), LDT, C(1,1), LDC, C(1,I), LDC, WORK, INFO )

          END DO
*
*         Multiply Q to the first block of C (1:M,1:MB)
*
          CALL ZGEMQRT('R','C',M , MB, K, NB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO )
*
      ELSE IF (RIGHT.AND.NOTRAN) THEN
*
*         Multiply Q to the first block of C
*
         KK = MOD((N-K),(MB-K))
         II=N-KK+1
         CTR = 1
         CALL ZGEMQRT('R','N', M, MB , K, NB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO )
*
         DO I=MB+1,II-MB+K,(MB-K)
*
*         Multiply Q to the current block of C (1:M,I:I+MB)
*
          CALL ZTPMQRT('R','N', M, MB-K, K, 0,NB, A(I,1), LDA, T(1, CTR * K + 1),LDT, C(1,1), LDC, C(1,I), LDC, WORK, INFO )
          CTR = CTR + 1
*
         END DO
         IF(II.LE.N) THEN
*
*         Multiply Q to the last block of C
*
          CALL ZTPMQRT('R','N', M, KK , K, 0,NB, A(II,1), LDA, T(1,CTR * K + 1),LDT, C(1,1), LDC, C(1,II), LDC, WORK, INFO )
*
         END IF
*
      END IF
*
      WORK( 1 ) = LWMIN
      RETURN
*
*     End of ZLAMTSQR
*
      END
