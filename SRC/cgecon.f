      SUBROUTINE CGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             NORM;
      int                INFO, LDA, N;
      REAL               ANORM, RCOND
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               ONENRM;
      String             NORMIN;
      int                IX, KASE, KASE1;
      REAL               AINVNM, SCALE, SL, SMLNUM, SU, HUGEVAL
      COMPLEX            ZDUM
*     ..
*     .. Local Arrays ..
      int                ISAVE( 3 );
*     ..
*     .. External Functions ..
      bool               LSAME, SISNAN;
      int                ICAMAX;
      REAL               SLAMCH
      EXTERNAL           LSAME, ICAMAX, SLAMCH, SISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACN2, CLATRS, CSRSCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      HUGEVAL = SLAMCH( 'Overflow' )
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
         CALL XERBLA( 'CGECON', -INFO )
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
      ELSE IF( SISNAN( ANORM ) ) THEN
         RCOND = ANORM
         INFO = -5
         RETURN
      ELSE IF( ANORM.GT.HUGEVAL ) THEN
         INFO = -5
         RETURN
      END IF
*
      SMLNUM = SLAMCH( 'Safe minimum' )
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
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
*
*           Multiply by inv(L).
*
            CALL CLATRS( 'Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, RWORK, INFO )
*
*           Multiply by inv(U).
*
            CALL CLATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, RWORK( N+1 ), INFO )
         ELSE
*
*           Multiply by inv(U**H).
*
            CALL CLATRS( 'Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, RWORK( N+1 ), INFO )
*
*           Multiply by inv(L**H).
*
            CALL CLATRS( 'Lower', 'Conjugate transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, RWORK, INFO )
         END IF
*
*        Divide X by 1/(SL*SU) if doing so will not cause overflow.
*
         SCALE = SL*SU
         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = ICAMAX( N, WORK, 1 )
            IF( SCALE.LT.CABS1( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
            CALL CSRSCL( N, SCALE, WORK, 1 )
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
      IF( SISNAN( RCOND ) .OR. RCOND.GT.HUGEVAL ) INFO = 1
*
   20 CONTINUE
      RETURN
*
*     End of CGECON
*
      END
