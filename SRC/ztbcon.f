      SUBROUTINE ZTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                INFO, KD, LDAB, N;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         AB( LDAB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, ONENRM, UPPER;
      String             NORMIN;
      int                IX, KASE, KASE1;
      double             AINVNM, ANORM, SCALE, SMLNUM, XNORM;
      COMPLEX*16         ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH, ZLANTB;
      // EXTERNAL LSAME, IZAMAX, DLAMCH, ZLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDRSCL, ZLACN2, ZLATBS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
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
      ELSE IF( KD.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTBCON', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      END IF

      RCOND = ZERO
      SMLNUM = DLAMCH( 'Safe minimum' )*DBLE( MAX( N, 1 ) )

      // Compute the 1-norm of the triangular matrix A or A**H.

      ANORM = ZLANTB( NORM, UPLO, DIAG, N, KD, AB, LDAB, RWORK )

      // Continue only if ANORM > 0.

      IF( ANORM.GT.ZERO ) THEN

         // Estimate the 1-norm of the inverse of A.

         AINVNM = ZERO
         NORMIN = 'N'
         IF( ONENRM ) THEN
            KASE1 = 1
         } else {
            KASE1 = 2
         END IF
         KASE = 0
   10    CONTINUE
         CALL ZLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.KASE1 ) THEN

               // Multiply by inv(A).

               CALL ZLATBS( UPLO, 'No transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, RWORK, INFO )
            } else {

               // Multiply by inv(A**H).

               CALL ZLATBS( UPLO, 'Conjugate transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, RWORK, INFO )
            END IF
            NORMIN = 'Y'

            // Multiply by 1/SCALE if doing so will not cause overflow.

            IF( SCALE.NE.ONE ) THEN
               IX = IZAMAX( N, WORK, 1 )
               XNORM = CABS1( WORK( IX ) )
               IF( SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
               CALL ZDRSCL( N, SCALE, WORK, 1 )
            END IF
            GO TO 10
         END IF

         // Compute the estimate of the reciprocal condition number.

         IF( AINVNM.NE.ZERO ) RCOND = ( ONE / ANORM ) / AINVNM
      END IF

   20 CONTINUE
      RETURN

      // End of ZTBCON

      }
