      SUBROUTINE ZGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             NORM;
      int                INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               ONENRM;
      String             NORMIN;
      int                IX, KASE, KASE1
      DOUBLE PRECISION   AINVNM, SCALE, SL, SMLNUM, SU, HUGEVAL
      COMPLEX*16         ZDUM
*     ..
*     .. Local Arrays ..
      int                ISAVE( 3 )
*     ..
*     .. External Functions ..
      bool               LSAME, DISNAN;
      int                IZAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IZAMAX, DLAMCH, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDRSCL, ZLACN2, ZLATRS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      HUGEVAL = DLAMCH( 'Overflow' )
*
*     Test the input parameters.
*
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGECON', -INFO )
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
      ELSE IF( DISNAN( ANORM ) ) THEN
         RCOND = ANORM
         INFO = -5
         RETURN
      ELSE IF( ANORM.GT.HUGEVAL ) THEN
         INFO = -5
         RETURN
      END IF
*
      SMLNUM = DLAMCH( 'Safe minimum' )
*
*     Estimate the norm of inv(A).
*
      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   10 CONTINUE
      CALL ZLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
*
*           Multiply by inv(L).
*
            CALL ZLATRS( 'Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, RWORK, INFO )
*
*           Multiply by inv(U).
*
            CALL ZLATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, RWORK( N+1 ), INFO )
         ELSE
*
*           Multiply by inv(U**H).
*
            CALL ZLATRS( 'Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, RWORK( N+1 ), INFO )
*
*           Multiply by inv(L**H).
*
            CALL ZLATRS( 'Lower', 'Conjugate transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, RWORK, INFO )
         END IF
*
*        Divide X by 1/(SL*SU) if doing so will not cause overflow.
*
         SCALE = SL*SU
         NORMIN = 'Y'
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
      IF( AINVNM.NE.ZERO ) THEN
         RCOND = ( ONE / AINVNM ) / ANORM
      ELSE
         INFO = 1
         RETURN
      END IF
*
*     Check for NaNs and Infs
*
      IF( DISNAN( RCOND ) .OR. RCOND.GT.HUGEVAL ) INFO = 1
*
   20 CONTINUE
      RETURN
*
*     End of ZGECON
*
      END
