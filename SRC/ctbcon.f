      SUBROUTINE CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            AB( LDAB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, ONENRM, UPPER
      CHARACTER          NORMIN
      INTEGER            IX, KASE, KASE1
      REAL               AINVNM, ANORM, SCALE, SMLNUM, XNORM
      COMPLEX            ZDUM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICAMAX
      REAL               CLANTB, SLAMCH
      EXTERNAL           LSAME, ICAMAX, CLANTB, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACN2, CLATBS, CSRSCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, MAX, REAL
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      NOUNIT = LSAME( DIAG, 'N' )
*
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KD.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTBCON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      END IF
*
      RCOND = ZERO
      SMLNUM = SLAMCH( 'Safe minimum' )*REAL( MAX( N, 1 ) )
*
*     Compute the 1-norm of the triangular matrix A or A**H.
*
      ANORM = CLANTB( NORM, UPLO, DIAG, N, KD, AB, LDAB, RWORK )
*
*     Continue only if ANORM > 0.
*
      IF( ANORM.GT.ZERO ) THEN
*
*        Estimate the 1-norm of the inverse of A.
*
         AINVNM = ZERO
         NORMIN = 'N'
         IF( ONENRM ) THEN
            KASE1 = 1
         ELSE
            KASE1 = 2
         END IF
         KASE = 0
   10    CONTINUE
         CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.KASE1 ) THEN
*
*              Multiply by inv(A).
*
               CALL CLATBS( UPLO, 'No transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, RWORK, INFO )
            ELSE
*
*              Multiply by inv(A**H).
*
               CALL CLATBS( UPLO, 'Conjugate transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, RWORK, INFO )
            END IF
            NORMIN = 'Y'
*
*           Multiply by 1/SCALE if doing so will not cause overflow.
*
            IF( SCALE.NE.ONE ) THEN
               IX = ICAMAX( N, WORK, 1 )
               XNORM = CABS1( WORK( IX ) )
               IF( SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
               CALL CSRSCL( N, SCALE, WORK, 1 )
            END IF
            GO TO 10
         END IF
*
*        Compute the estimate of the reciprocal condition number.
*
         IF( AINVNM.NE.ZERO ) RCOND = ( ONE / ANORM ) / AINVNM
      END IF
*
   20 CONTINUE
      RETURN
*
*     End of CTBCON
*
      END
