      SUBROUTINE SGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, LDA, N;
      REAL               ANORM, RCOND
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               ONENRM;
      String             NORMIN;
      int                IX, KASE, KASE1;
      REAL               AINVNM, SCALE, SL, SMLNUM, SU, HUGEVAL
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      int                ISAMAX;
      REAL               SLAMCH
      // EXTERNAL LSAME, ISAMAX, SLAMCH, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SLATRS, SRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      HUGEVAL = SLAMCH( 'Overflow' )

      // Test the input parameters.

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
         CALL XERBLA( 'SGECON', -INFO )
         RETURN
      END IF

      // Quick return if possible

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

      SMLNUM = SLAMCH( 'Safe minimum' )

      // Estimate the norm of inv(A).

      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   10 CONTINUE
      CALL SLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN

            // Multiply by inv(L).

            CALL SLATRS( 'Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, WORK( 2*N+1 ), INFO )

            // Multiply by inv(U).

            CALL SLATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, WORK( 3*N+1 ), INFO )
         ELSE

            // Multiply by inv(U**T).

            CALL SLATRS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, WORK( 3*N+1 ), INFO )

            // Multiply by inv(L**T).

            CALL SLATRS( 'Lower', 'Transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, WORK( 2*N+1 ), INFO )
         END IF

         // Divide X by 1/(SL*SU) if doing so will not cause overflow.

         SCALE = SL*SU
         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = ISAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
            CALL SRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM.NE.ZERO ) THEN
         RCOND = ( ONE / AINVNM ) / ANORM
      ELSE
         INFO = 1
         RETURN
      END IF

      // Check for NaNs and Infs

      IF( SISNAN( RCOND ) .OR. RCOND.GT.HUGEVAL ) INFO = 1

   20 CONTINUE
      RETURN

      // End of SGECON

      }
