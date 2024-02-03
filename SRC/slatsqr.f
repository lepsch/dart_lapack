      SUBROUTINE SLATSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
*
*     .. Scalar Arguments ..
      int                INFO, LDA, M, N, MB, NB, LDT, LWORK
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * ), T( LDT, * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      bool               LQUERY;
      int                I, II, KK, CTR, MINMN, LWMIN
*     ..
*     .. EXTERNAL FUNCTIONS ..
      bool               LSAME;
      EXTERNAL           LSAME
      REAL               SROUNDUP_LWORK
      EXTERNAL           SROUNDUP_LWORK
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           SGEQRT, STPQRT, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN, MOD
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
        LWMIN = N*NB
      END IF
*
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 .OR. M.LT.N ) THEN
        INFO = -2
      ELSE IF( MB.LT.1 ) THEN
        INFO = -3
      ELSE IF( NB.LT.1 .OR. ( NB.GT.N .AND. N.GT.0 ) ) THEN
        INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -6
      ELSE IF( LDT.LT.NB ) THEN
        INFO = -8
      ELSE IF( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) THEN
        INFO = -10
      END IF
*
      IF( INFO.EQ.0 )  THEN
        WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'SLATSQR', -INFO )
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
*     The QR Decomposition
*
      IF( (MB.LE.N) .OR. (MB.GE.M) ) THEN
        CALL SGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
        RETURN
      END IF
      KK = MOD((M-N),(MB-N))
      II = M-KK+1
*
*     Compute the QR factorization of the first block A(1:MB,1:N)
*
      CALL SGEQRT( MB, N, NB, A(1,1), LDA, T, LDT, WORK, INFO )
*
      CTR = 1
      DO I = MB+1, II-MB+N, (MB-N)
*
*       Compute the QR factorization of the current block A(I:I+MB-N,1:N)
*
        CALL STPQRT( MB-N, N, 0, NB, A(1,1), LDA, A( I, 1 ), LDA, T(1, CTR * N + 1), LDT, WORK, INFO )
        CTR = CTR + 1
      END DO
*
*     Compute the QR factorization of the last block A(II:M,1:N)
*
      IF( II.LE.M ) THEN
        CALL STPQRT( KK, N, 0, NB, A(1,1), LDA, A( II, 1 ), LDA, T(1, CTR * N + 1), LDT, WORK, INFO )
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      RETURN
*
*     End of SLATSQR
*
      END
