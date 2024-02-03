      SUBROUTINE DGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, KL, KU, LDAB, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             AB( LDAB, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               LNOTI, ONENRM;
      String             NORMIN;
      int                IX, J, JP, KASE, KASE1, KD, LM;
      double             AINVNM, SCALE, SMLNUM, T;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DDOT, DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DDOT, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLACN2, DLATBS, DRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
         INFO = -6
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBCON', -INFO )
         RETURN
      END IF

      // Quick return if possible

      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF

      SMLNUM = DLAMCH( 'Safe minimum' )

      // Estimate the norm of inv(A).

      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KD = KL + KU + 1
      LNOTI = KL.GT.0
      KASE = 0
   10 CONTINUE
      CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN

            // Multiply by inv(L).

            IF( LNOTI ) THEN
               DO 20 J = 1, N - 1
                  LM = MIN( KL, N-J )
                  JP = IPIV( J )
                  T = WORK( JP )
                  IF( JP.NE.J ) THEN
                     WORK( JP ) = WORK( J )
                     WORK( J ) = T
                  END IF
                  CALL DAXPY( LM, -T, AB( KD+1, J ), 1, WORK( J+1 ), 1 )
   20          CONTINUE
            END IF

            // Multiply by inv(U).

            CALL DLATBS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N, KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO )
         ELSE

            // Multiply by inv(U**T).

            CALL DLATBS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N, KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO )

            // Multiply by inv(L**T).

            IF( LNOTI ) THEN
               DO 30 J = N - 1, 1, -1
                  LM = MIN( KL, N-J )
                  WORK( J ) = WORK( J ) - DDOT( LM, AB( KD+1, J ), 1, WORK( J+1 ), 1 )
                  JP = IPIV( J )
                  IF( JP.NE.J ) THEN
                     T = WORK( JP )
                     WORK( JP ) = WORK( J )
                     WORK( J ) = T
                  END IF
   30          CONTINUE
            END IF
         END IF

         // Divide X by 1/SCALE if doing so will not cause overflow.

         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 40
            CALL DRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM

   40 CONTINUE
      RETURN

      // End of DGBCON

      END
