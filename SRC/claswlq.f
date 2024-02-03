      SUBROUTINE CLASWLQ( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
*
*     .. Scalar Arguments ..
      int                INFO, LDA, M, N, MB, NB, LWORK, LDT;
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), WORK( * ), T( LDT, * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      bool               LQUERY;
      int                I, II, KK, CTR, MINMN, LWMIN;
*     ..
*     .. EXTERNAL FUNCTIONS ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SROUNDUP_LWORK
*     ..
*     .. EXTERNAL SUBROUTINES ..
      // EXTERNAL CGELQT, CTPLQT, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      // INTRINSIC MAX, MIN, MOD
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
*
      LQUERY = ( LWORK.EQ.-1 )
*
      MINMN = MIN( M, N )
      IF( MINMN.EQ.0 ) THEN
        LWMIN = 1
      ELSE
        LWMIN = M*MB
      END IF
*
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 .OR. N.LT.M ) THEN
        INFO = -2
      ELSE IF( MB.LT.1 .OR. ( MB.GT.M .AND. M.GT.0 ) ) THEN
        INFO = -3
      ELSE IF( NB.LE.0 ) THEN
        INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -6
      ELSE IF( LDT.LT.MB ) THEN
        INFO = -8
      ELSE IF( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) THEN
        INFO = -10
      END IF
*
      IF( INFO.EQ.0 ) THEN
        WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      END IF
*
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'CLASWLQ', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( MINMN.EQ.0 ) THEN
        RETURN
      END IF
*
*     The LQ Decomposition
*
      IF( (M.GE.N) .OR. (NB.LE.M) .OR. (NB.GE.N) ) THEN
        CALL CGELQT( M, N, MB, A, LDA, T, LDT, WORK, INFO)
        RETURN
      END IF
*
      KK = MOD((N-M),(NB-M))
      II = N-KK+1
*
*     Compute the LQ factorization of the first block A(1:M,1:NB)
*
      CALL CGELQT( M, NB, MB, A(1,1), LDA, T, LDT, WORK, INFO)
      CTR = 1
*
      DO I = NB+1, II-NB+M , (NB-M)
*
*       Compute the QR factorization of the current block A(1:M,I:I+NB-M)
*
        CALL CTPLQT( M, NB-M, 0, MB, A(1,1), LDA, A( 1, I ), LDA, T(1,CTR*M+1), LDT, WORK, INFO )
        CTR = CTR + 1
      END DO
*
*     Compute the QR factorization of the last block A(1:M,II:N)
*
      IF( II.LE.N ) THEN
        CALL CTPLQT( M, KK, 0, MB, A(1,1), LDA, A( 1, II ), LDA, T(1,CTR*M+1), LDT, WORK, INFO )
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      RETURN
*
*     End of CLASWLQ
*
      END
