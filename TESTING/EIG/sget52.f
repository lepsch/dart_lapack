      SUBROUTINE SGET52( LEFT, N, A, LDA, B, LDB, E, LDE, ALPHAR,
     $                   ALPHAI, BETA, WORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            LEFT
      INTEGER            LDA, LDB, LDE, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), E( LDE, * ),
     $                   RESULT( 2 ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TEN
      PARAMETER          ( ZERO = 0.0, ONE = 1.0, TEN = 10.0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILCPLX
      CHARACTER          NORMAB, TRANS
      INTEGER            J, JVEC
      REAL               ABMAX, ACOEF, ALFMAX, ANORM, BCOEFI, BCOEFR,
     $                   BETMAX, BNORM, ENORM, ENRMER, ERRNRM, SAFMAX,
     $                   SAFMIN, SALFI, SALFR, SBETA, SCALE, TEMP1, ULP
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, REAL
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFMAX = ONE / SAFMIN
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
*
      IF( LEFT ) THEN
         TRANS = 'T'
         NORMAB = 'I'
      ELSE
         TRANS = 'N'
         NORMAB = 'O'
      END IF
*
*     Norm of A, B, and E:
*
      ANORM = MAX( SLANGE( NORMAB, N, N, A, LDA, WORK ), SAFMIN )
      BNORM = MAX( SLANGE( NORMAB, N, N, B, LDB, WORK ), SAFMIN )
      ENORM = MAX( SLANGE( 'O', N, N, E, LDE, WORK ), ULP )
      ALFMAX = SAFMAX / MAX( ONE, BNORM )
      BETMAX = SAFMAX / MAX( ONE, ANORM )
*
*     Compute error matrix.
*     Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )
*
      ILCPLX = .FALSE.
      DO 10 JVEC = 1, N
         IF( ILCPLX ) THEN
*
*           2nd Eigenvalue/-vector of pair -- do nothing
*
            ILCPLX = .FALSE.
         ELSE
            SALFR = ALPHAR( JVEC )
            SALFI = ALPHAI( JVEC )
            SBETA = BETA( JVEC )
            IF( SALFI.EQ.ZERO ) THEN
*
*              Real eigenvalue and -vector
*
               ABMAX = MAX( ABS( SALFR ), ABS( SBETA ) )
               IF( ABS( SALFR ).GT.ALFMAX .OR. ABS( SBETA ).GT.
     $             BETMAX .OR. ABMAX.LT.ONE ) THEN
                  SCALE = ONE / MAX( ABMAX, SAFMIN )
                  SALFR = SCALE*SALFR
                  SBETA = SCALE*SBETA
               END IF
               SCALE = ONE / MAX( ABS( SALFR )*BNORM,
     $                 ABS( SBETA )*ANORM, SAFMIN )
               ACOEF = SCALE*SBETA
               BCOEFR = SCALE*SALFR
               CALL SGEMV( TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1,
     $                     ZERO, WORK( N*( JVEC-1 )+1 ), 1 )
               CALL SGEMV( TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ),
     $                     1, ONE, WORK( N*( JVEC-1 )+1 ), 1 )
            ELSE
*
*              Complex conjugate pair
*
               ILCPLX = .TRUE.
               IF( JVEC.EQ.N ) THEN
                  RESULT( 1 ) = TEN / ULP
                  RETURN
               END IF
               ABMAX = MAX( ABS( SALFR )+ABS( SALFI ), ABS( SBETA ) )
               IF( ABS( SALFR )+ABS( SALFI ).GT.ALFMAX .OR.
     $             ABS( SBETA ).GT.BETMAX .OR. ABMAX.LT.ONE ) THEN
                  SCALE = ONE / MAX( ABMAX, SAFMIN )
                  SALFR = SCALE*SALFR
                  SALFI = SCALE*SALFI
                  SBETA = SCALE*SBETA
               END IF
               SCALE = ONE / MAX( ( ABS( SALFR )+ABS( SALFI ) )*BNORM,
     $                 ABS( SBETA )*ANORM, SAFMIN )
               ACOEF = SCALE*SBETA
               BCOEFR = SCALE*SALFR
               BCOEFI = SCALE*SALFI
               IF( LEFT ) THEN
                  BCOEFI = -BCOEFI
               END IF
*
               CALL SGEMV( TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1,
     $                     ZERO, WORK( N*( JVEC-1 )+1 ), 1 )
               CALL SGEMV( TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ),
     $                     1, ONE, WORK( N*( JVEC-1 )+1 ), 1 )
               CALL SGEMV( TRANS, N, N, BCOEFI, B, LDA, E( 1, JVEC+1 ),
     $                     1, ONE, WORK( N*( JVEC-1 )+1 ), 1 )
*
               CALL SGEMV( TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC+1 ),
     $                     1, ZERO, WORK( N*JVEC+1 ), 1 )
               CALL SGEMV( TRANS, N, N, -BCOEFI, B, LDA, E( 1, JVEC ),
     $                     1, ONE, WORK( N*JVEC+1 ), 1 )
               CALL SGEMV( TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC+1 ),
     $                     1, ONE, WORK( N*JVEC+1 ), 1 )
            END IF
         END IF
   10 CONTINUE
*
      ERRNRM = SLANGE( 'One', N, N, WORK, N, WORK( N**2+1 ) ) / ENORM
*
*     Compute RESULT(1)
*
      RESULT( 1 ) = ERRNRM / ULP
*
*     Normalization of E:
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
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) )+
     $                    ABS( E( J, JVEC+1 ) ) )
   30          CONTINUE
               ENRMER = MAX( ENRMER, ABS( TEMP1-ONE ) )
            END IF
         END IF
   40 CONTINUE
*
*     Compute RESULT(2) : the normalization error in E.
*
      RESULT( 2 ) = ENRMER / ( REAL( N )*ULP )
*
      RETURN
*
*     End of SGET52
*
      END