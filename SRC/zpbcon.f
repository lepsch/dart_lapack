      SUBROUTINE ZPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int                INFO, KD, LDAB, N
      DOUBLE PRECISION   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         AB( LDAB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      CHARACTER          NORMIN
      int                IX, KASE
      DOUBLE PRECISION   AINVNM, SCALE, SCALEL, SCALEU, SMLNUM
      COMPLEX*16         ZDUM
*     ..
*     .. Local Arrays ..
      int                ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                IZAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IZAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDRSCL, ZLACN2, ZLATBS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPBCON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      SMLNUM = DLAMCH( 'Safe minimum' )
*
*     Estimate the 1-norm of the inverse.
*
      KASE = 0
      NORMIN = 'N'
   10 CONTINUE
      CALL ZLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( UPPER ) THEN
*
*           Multiply by inv(U**H).
*
            CALL ZLATBS( 'Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, RWORK, INFO )
            NORMIN = 'Y'
*
*           Multiply by inv(U).
*
            CALL ZLATBS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, RWORK, INFO )
         ELSE
*
*           Multiply by inv(L).
*
            CALL ZLATBS( 'Lower', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, RWORK, INFO )
            NORMIN = 'Y'
*
*           Multiply by inv(L**H).
*
            CALL ZLATBS( 'Lower', 'Conjugate transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, RWORK, INFO )
         END IF
*
*        Multiply by 1/SCALE if doing so will not cause overflow.
*
         SCALE = SCALEL*SCALEU
         IF( SCALE.NE.ONE ) THEN
            IX = IZAMAX( N, WORK, 1 )
            IF( SCALE.LT.CABS1( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
            CALL ZDRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM
*
   20 CONTINUE
*
      RETURN
*
*     End of ZPBCON
*
      END
