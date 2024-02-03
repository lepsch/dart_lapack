      SUBROUTINE CPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      REAL               ANORM, RCOND
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            AB( LDAB, * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             NORMIN;
      int                IX, KASE;
      REAL               AINVNM, SCALE, SCALEL, SCALEU, SMLNUM
      COMPLEX            ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SLAMCH
      // EXTERNAL LSAME, ICAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CLATBS, CSRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
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
         CALL XERBLA( 'CPBCON', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      SMLNUM = SLAMCH( 'Safe minimum' )
*
      // Estimate the 1-norm of the inverse.
*
      KASE = 0
      NORMIN = 'N'
   10 CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( UPPER ) THEN
*
            // Multiply by inv(U**H).
*
            CALL CLATBS( 'Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, RWORK, INFO )
            NORMIN = 'Y'
*
            // Multiply by inv(U).
*
            CALL CLATBS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, RWORK, INFO )
         ELSE
*
            // Multiply by inv(L).
*
            CALL CLATBS( 'Lower', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, RWORK, INFO )
            NORMIN = 'Y'
*
            // Multiply by inv(L**H).
*
            CALL CLATBS( 'Lower', 'Conjugate transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, RWORK, INFO )
         END IF
*
         // Multiply by 1/SCALE if doing so will not cause overflow.
*
         SCALE = SCALEL*SCALEU
         IF( SCALE.NE.ONE ) THEN
            IX = ICAMAX( N, WORK, 1 )
            IF( SCALE.LT.CABS1( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
            CALL CSRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
*
      // Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM
*
   20 CONTINUE
*
      RETURN
*
      // End of CPBCON
*
      END
