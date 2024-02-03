      SUBROUTINE DGET52( LEFT, N, A, LDA, B, LDB, E, LDE, ALPHAR, ALPHAI, BETA, WORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      bool               LEFT;
      int                LDA, LDB, LDE, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), E( LDE, * ), RESULT( 2 ), WORK( * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE, TEN;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TEN = 10.0D0 )
      // ..
      // .. Local Scalars ..
      bool               ILCPLX;
      String             NORMAB, TRANS;
      int                J, JVEC;
      double             ABMAX, ACOEF, ALFMAX, ANORM, BCOEFI, BCOEFR, BETMAX, BNORM, ENORM, ENRMER, ERRNRM, SAFMAX, SAFMIN, SALFI, SALFR, SBETA, SCALE, TEMP1, ULP;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX
      // ..
      // .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 ) RETURN
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
*
      IF( LEFT ) THEN
         TRANS = 'T'
         NORMAB = 'I'
      ELSE
         TRANS = 'N'
         NORMAB = 'O'
      END IF
*
      // Norm of A, B, and E:
*
      ANORM = MAX( DLANGE( NORMAB, N, N, A, LDA, WORK ), SAFMIN )
      BNORM = MAX( DLANGE( NORMAB, N, N, B, LDB, WORK ), SAFMIN )
      ENORM = MAX( DLANGE( 'O', N, N, E, LDE, WORK ), ULP )
      ALFMAX = SAFMAX / MAX( ONE, BNORM )
      BETMAX = SAFMAX / MAX( ONE, ANORM )
*
      // Compute error matrix.
      // Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )
*
      ILCPLX = .FALSE.
      DO 10 JVEC = 1, N
         IF( ILCPLX ) THEN
*
            // 2nd Eigenvalue/-vector of pair -- do nothing
*
            ILCPLX = .FALSE.
         ELSE
            SALFR = ALPHAR( JVEC )
            SALFI = ALPHAI( JVEC )
            SBETA = BETA( JVEC )
            IF( SALFI.EQ.ZERO ) THEN
*
               // Real eigenvalue and -vector
*
               ABMAX = MAX( ABS( SALFR ), ABS( SBETA ) )
               IF( ABS( SALFR ).GT.ALFMAX .OR. ABS( SBETA ).GT. BETMAX .OR. ABMAX.LT.ONE ) THEN
                  SCALE = ONE / MAX( ABMAX, SAFMIN )
                  SALFR = SCALE*SALFR
                  SBETA = SCALE*SBETA
               END IF
               SCALE = ONE / MAX( ABS( SALFR )*BNORM, ABS( SBETA )*ANORM, SAFMIN )
               ACOEF = SCALE*SBETA
               BCOEFR = SCALE*SALFR
               CALL DGEMV( TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1, ZERO, WORK( N*( JVEC-1 )+1 ), 1 )                CALL DGEMV( TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 )
            ELSE
*
               // Complex conjugate pair
*
               ILCPLX = .TRUE.
               IF( JVEC.EQ.N ) THEN
                  RESULT( 1 ) = TEN / ULP
                  RETURN
               END IF
               ABMAX = MAX( ABS( SALFR )+ABS( SALFI ), ABS( SBETA ) )
               IF( ABS( SALFR )+ABS( SALFI ).GT.ALFMAX .OR. ABS( SBETA ).GT.BETMAX .OR. ABMAX.LT.ONE ) THEN
                  SCALE = ONE / MAX( ABMAX, SAFMIN )
                  SALFR = SCALE*SALFR
                  SALFI = SCALE*SALFI
                  SBETA = SCALE*SBETA
               END IF
               SCALE = ONE / MAX( ( ABS( SALFR )+ABS( SALFI ) )*BNORM, ABS( SBETA )*ANORM, SAFMIN )
               ACOEF = SCALE*SBETA
               BCOEFR = SCALE*SALFR
               BCOEFI = SCALE*SALFI
               IF( LEFT ) THEN
                  BCOEFI = -BCOEFI
               END IF
*
               CALL DGEMV( TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1, ZERO, WORK( N*( JVEC-1 )+1 ), 1 )                CALL DGEMV( TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 )                CALL DGEMV( TRANS, N, N, BCOEFI, B, LDA, E( 1, JVEC+1 ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 )
*
               CALL DGEMV( TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC+1 ), 1, ZERO, WORK( N*JVEC+1 ), 1 )                CALL DGEMV( TRANS, N, N, -BCOEFI, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*JVEC+1 ), 1 )                CALL DGEMV( TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC+1 ), 1, ONE, WORK( N*JVEC+1 ), 1 )
            END IF
         END IF
   10 CONTINUE
*
      ERRNRM = DLANGE( 'One', N, N, WORK, N, WORK( N**2+1 ) ) / ENORM
*
      // Compute RESULT(1)
*
      RESULT( 1 ) = ERRNRM / ULP
*
      // Normalization of E:
*
      ENRMER = ZERO
      ILCPLX = .FALSE.
      DO 40 JVEC = 1, N
         IF( ILCPLX ) THEN
            ILCPLX = .FALSE.
         ELSE
            TEMP1 = ZERO
            IF( ALPHAI( JVEC ).EQ.ZERO ) THEN
               DO 20 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) ) )
   20          CONTINUE
               ENRMER = MAX( ENRMER, ABS( TEMP1-ONE ) )
            ELSE
               ILCPLX = .TRUE.
               DO 30 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) )+ ABS( E( J, JVEC+1 ) ) )
   30          CONTINUE
               ENRMER = MAX( ENRMER, ABS( TEMP1-ONE ) )
            END IF
         END IF
   40 CONTINUE
*
      // Compute RESULT(2) : the normalization error in E.
*
      RESULT( 2 ) = ENRMER / ( DBLE( N )*ULP )
*
      RETURN
*
      // End of DGET52
*
      END
