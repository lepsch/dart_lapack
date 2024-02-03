      REAL FUNCTION CLA_SYRCOND_X( UPLO, N, A, LDA, AF, LDAF, IPIV, X, INFO, WORK, RWORK )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                N, LDA, LDAF, INFO
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )
      REAL               RWORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                KASE
      REAL               AINVNM, ANORM, TMP
      int                I, J
      bool               UP, UPPER;
      COMPLEX            ZDUM
*     ..
*     .. Local Arrays ..
      int                ISAVE( 3 )
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACN2, CSYTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      CLA_SYRCOND_X = 0.0E+0
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF ( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLA_SYRCOND_X', -INFO )
         RETURN
      END IF
      UP = .FALSE.
      IF ( LSAME( UPLO, 'U' ) ) UP = .TRUE.
*
*     Compute norm of op(A)*op2(C).
*
      ANORM = 0.0
      IF ( UP ) THEN
         DO I = 1, N
            TMP = 0.0E+0
            DO J = 1, I
               TMP = TMP + CABS1( A( J, I ) * X( J ) )
            END DO
            DO J = I+1, N
               TMP = TMP + CABS1( A( I, J ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0E+0
            DO J = 1, I
               TMP = TMP + CABS1( A( I, J ) * X( J ) )
            END DO
            DO J = I+1, N
               TMP = TMP + CABS1( A( J, I ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         CLA_SYRCOND_X = 1.0E+0
         RETURN
      ELSE IF( ANORM .EQ. 0.0E+0 ) THEN
         RETURN
      END IF
*
*     Estimate the norm of inv(op(A)).
*
      AINVNM = 0.0E+0
*
      KASE = 0
   10 CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.2 ) THEN
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
*
            IF ( UP ) THEN
               CALL CSYTRS( 'U', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ELSE
               CALL CSYTRS( 'L', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ENDIF
*
*           Multiply by inv(X).
*
            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO
         ELSE
*
*           Multiply by inv(X**T).
*
            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO
*
            IF ( UP ) THEN
               CALL CSYTRS( 'U', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ELSE
               CALL CSYTRS( 'L', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            END IF
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM .NE. 0.0E+0 ) CLA_SYRCOND_X = 1.0E+0 / AINVNM
*
      RETURN
*
*     End of CLA_SYRCOND_X
*
      END
