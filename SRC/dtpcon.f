      SUBROUTINE DTPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                INFO, N;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             AP( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, ONENRM, UPPER;
      String             NORMIN;
      int                IX, KASE, KASE1;
      double             AINVNM, ANORM, SCALE, SMLNUM, XNORM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH, DLANTP;
      // EXTERNAL LSAME, IDAMAX, DLAMCH, DLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLATPS, DRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      NOUNIT = LSAME( DIAG, 'N' )

      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTPCON', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      END IF

      RCOND = ZERO
      SMLNUM = DLAMCH( 'Safe minimum' )*DBLE( MAX( 1, N ) )

      // Compute the norm of the triangular matrix A.

      ANORM = DLANTP( NORM, UPLO, DIAG, N, AP, WORK )

      // Continue only if ANORM > 0.

      IF( ANORM.GT.ZERO ) THEN

         // Estimate the norm of the inverse of A.

         AINVNM = ZERO
         NORMIN = 'N'
         IF( ONENRM ) THEN
            KASE1 = 1
         ELSE
            KASE1 = 2
         END IF
         KASE = 0
   10    CONTINUE
         CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.KASE1 ) THEN

               // Multiply by inv(A).

               CALL DLATPS( UPLO, 'No transpose', DIAG, NORMIN, N, AP, WORK, SCALE, WORK( 2*N+1 ), INFO )
            ELSE

               // Multiply by inv(A**T).

               CALL DLATPS( UPLO, 'Transpose', DIAG, NORMIN, N, AP, WORK, SCALE, WORK( 2*N+1 ), INFO )
            END IF
            NORMIN = 'Y'

            // Multiply by 1/SCALE if doing so will not cause overflow.

            IF( SCALE.NE.ONE ) THEN
               IX = IDAMAX( N, WORK, 1 )
               XNORM = ABS( WORK( IX ) )
               IF( SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
               CALL DRSCL( N, SCALE, WORK, 1 )
            END IF
            GO TO 10
         END IF

         // Compute the estimate of the reciprocal condition number.

         IF( AINVNM.NE.ZERO ) RCOND = ( ONE / ANORM ) / AINVNM
      END IF

   20 CONTINUE
      RETURN

      // End of DTPCON

      END
