      SUBROUTINE DSYTRD_SY2SB( UPLO, N, KD, A, LDA, AB, LDAB, TAU,  WORK, LWORK, INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDAB, LWORK, N, KD
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), AB( LDAB, * ),  TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   RONE
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( RONE = 1.0D+0, ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                I, J, IINFO, LWMIN, PN, PK, LK, LDT, LDW, LDS2, LDS1, LS2, LS1, LW, LT, TPOS, WPOS, S2POS, S1POS
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DSYR2K, DSYMM, DGEMM, DCOPY, DLARFT, DGELQF, DGEQRF, DLASET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE
      EXTERNAL           LSAME, ILAENV2STAGE
*     ..
*     .. Executable Statements ..
*
*     Determine the minimal workspace size required
*     and test the input parameters
*
      INFO   = 0
      UPPER  = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LE.KD+1 ) THEN
         LWMIN = 1
      ELSE
         LWMIN = ILAENV2STAGE( 4, 'DSYTRD_SY2SB', ' ', N, KD, -1, -1 )
      END IF
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.MAX( 1, KD+1 ) ) THEN
         INFO = -7
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRD_SY2SB', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         WORK( 1 ) = LWMIN
         RETURN
      END IF
*
*     Quick return if possible
*     Copy the upper/lower portion of A into AB
*
      IF( N.LE.KD+1 ) THEN
          IF( UPPER ) THEN
              DO 100 I = 1, N
                  LK = MIN( KD+1, I )
                  CALL DCOPY( LK, A( I-LK+1, I ), 1,  AB( KD+1-LK+1, I ), 1 )
  100         CONTINUE
          ELSE
              DO 110 I = 1, N
                  LK = MIN( KD+1, N-I+1 )
                  CALL DCOPY( LK, A( I, I ), 1, AB( 1, I ), 1 )
  110         CONTINUE
          ENDIF
          WORK( 1 ) = 1
          RETURN
      END IF
*
*     Determine the pointer position for the workspace
*
      LDT    = KD
      LDS1   = KD
      LT     = LDT*KD
      LW     = N*KD
      LS1    = LDS1*KD
      LS2    = LWMIN - LT - LW - LS1
*      LS2 = N*MAX(KD,FACTOPTNB)
      TPOS   = 1
      WPOS   = TPOS  + LT
      S1POS  = WPOS  + LW
      S2POS  = S1POS + LS1
      IF( UPPER ) THEN
          LDW    = KD
          LDS2   = KD
      ELSE
          LDW    = N
          LDS2   = N
      ENDIF
*
*
*     Set the workspace of the triangular matrix T to zero once such a
*     way every time T is generated the upper/lower portion will be always zero
*
      CALL DLASET( "A", LDT, KD, ZERO, ZERO, WORK( TPOS ), LDT )
*
      IF( UPPER ) THEN
          DO 10 I = 1, N - KD, KD
             PN = N-I-KD+1
             PK = MIN( N-I-KD+1, KD )
*
*            Compute the LQ factorization of the current block
*
             CALL DGELQF( KD, PN, A( I, I+KD ), LDA, TAU( I ), WORK( S2POS ), LS2, IINFO )
*
*            Copy the upper portion of A into AB
*
             DO 20 J = I, I+PK-1
                LK = MIN( KD, N-J ) + 1
                CALL DCOPY( LK, A( J, J ), LDA, AB( KD+1, J ), LDAB-1 )
   20        CONTINUE
*
             CALL DLASET( 'Lower', PK, PK, ZERO, ONE,  A( I, I+KD ), LDA )
*
*            Form the matrix T
*
             CALL DLARFT( 'Forward', 'Rowwise', PN, PK, A( I, I+KD ), LDA, TAU( I ), WORK( TPOS ), LDT )
*
*            Compute W:
*
             CALL DGEMM( 'Conjugate', 'No transpose', PK, PN, PK, ONE,  WORK( TPOS ), LDT, A( I, I+KD ), LDA, ZERO, WORK( S2POS ), LDS2 )
*
             CALL DSYMM( 'Right', UPLO, PK, PN, ONE,  A( I+KD, I+KD ), LDA, WORK( S2POS ), LDS2, ZERO, WORK( WPOS ), LDW )
*
             CALL DGEMM( 'No transpose', 'Conjugate', PK, PK, PN, ONE,  WORK( WPOS ), LDW, WORK( S2POS ), LDS2, ZERO, WORK( S1POS ), LDS1 )
*
             CALL DGEMM( 'No transpose', 'No transpose', PK, PN, PK, -HALF, WORK( S1POS ), LDS1, A( I, I+KD ), LDA, ONE,   WORK( WPOS ), LDW )
*
*
*            Update the unreduced submatrix A(i+kd:n,i+kd:n), using
*            an update of the form:  A := A - V'*W - W'*V
*
             CALL DSYR2K( UPLO, 'Conjugate', PN, PK, -ONE, A( I, I+KD ), LDA, WORK( WPOS ), LDW, RONE, A( I+KD, I+KD ), LDA )
   10     CONTINUE
*
*        Copy the upper band to AB which is the band storage matrix
*
         DO 30 J = N-KD+1, N
            LK = MIN(KD, N-J) + 1
            CALL DCOPY( LK, A( J, J ), LDA, AB( KD+1, J ), LDAB-1 )
   30    CONTINUE
*
      ELSE
*
*         Reduce the lower triangle of A to lower band matrix
*
          DO 40 I = 1, N - KD, KD
             PN = N-I-KD+1
             PK = MIN( N-I-KD+1, KD )
*
*            Compute the QR factorization of the current block
*
             CALL DGEQRF( PN, KD, A( I+KD, I ), LDA, TAU( I ), WORK( S2POS ), LS2, IINFO )
*
*            Copy the upper portion of A into AB
*
             DO 50 J = I, I+PK-1
                LK = MIN( KD, N-J ) + 1
                CALL DCOPY( LK, A( J, J ), 1, AB( 1, J ), 1 )
   50        CONTINUE
*
             CALL DLASET( 'Upper', PK, PK, ZERO, ONE,  A( I+KD, I ), LDA )
*
*            Form the matrix T
*
             CALL DLARFT( 'Forward', 'Columnwise', PN, PK, A( I+KD, I ), LDA, TAU( I ), WORK( TPOS ), LDT )
*
*            Compute W:
*
             CALL DGEMM( 'No transpose', 'No transpose', PN, PK, PK, ONE, A( I+KD, I ), LDA, WORK( TPOS ), LDT, ZERO, WORK( S2POS ), LDS2 )
*
             CALL DSYMM( 'Left', UPLO, PN, PK, ONE, A( I+KD, I+KD ), LDA, WORK( S2POS ), LDS2, ZERO, WORK( WPOS ), LDW )
*
             CALL DGEMM( 'Conjugate', 'No transpose', PK, PK, PN, ONE, WORK( S2POS ), LDS2, WORK( WPOS ), LDW, ZERO, WORK( S1POS ), LDS1 )
*
             CALL DGEMM( 'No transpose', 'No transpose', PN, PK, PK, -HALF, A( I+KD, I ), LDA, WORK( S1POS ), LDS1, ONE, WORK( WPOS ), LDW )
*
*
*            Update the unreduced submatrix A(i+kd:n,i+kd:n), using
*            an update of the form:  A := A - V*W' - W*V'
*
             CALL DSYR2K( UPLO, 'No transpose', PN, PK, -ONE, A( I+KD, I ), LDA, WORK( WPOS ), LDW, RONE, A( I+KD, I+KD ), LDA )
*            ==================================================================
*            RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED
*             DO 45 J = I, I+PK-1
*                LK = MIN( KD, N-J ) + 1
*                CALL DCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 )
*   45        CONTINUE
*            ==================================================================
   40     CONTINUE
*
*        Copy the lower band to AB which is the band storage matrix
*
         DO 60 J = N-KD+1, N
            LK = MIN(KD, N-J) + 1
            CALL DCOPY( LK, A( J, J ), 1, AB( 1, J ), 1 )
   60    CONTINUE

      END IF
*
      WORK( 1 ) = LWMIN
      RETURN
*
*     End of DSYTRD_SY2SB
*
      END
