      SUBROUTINE SGET22( TRANSA, TRANSE, TRANSW, N, A, LDA, E, LDE, WR, WI, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSA, TRANSE, TRANSW;
      int                LDA, LDE, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), E( LDE, * ), RESULT( 2 ), WI( * ), WORK( * ), WR( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0, ONE = 1.0 )
      // ..
      // .. Local Scalars ..
      String             NORMA, NORME;
      int                IECOL, IEROW, INCE, IPAIR, ITRNSE, J, JCOL, JVEC       REAL               ANORM, ENORM, ENRMAX, ENRMIN, ERRNRM, TEMP1, ULP, UNFL;
      // ..
      // .. Local Arrays ..
      REAL               WMAT( 2, 2 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGE
      // EXTERNAL LSAME, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SGEMM, SLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Initialize RESULT (in case N=0)

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 ) RETURN

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )

      ITRNSE = 0
      INCE = 1
      NORMA = 'O'
      NORME = 'O'

      IF( LSAME( TRANSA, 'T' ) .OR. LSAME( TRANSA, 'C' ) ) THEN
         NORMA = 'I'
      END IF
      IF( LSAME( TRANSE, 'T' ) .OR. LSAME( TRANSE, 'C' ) ) THEN
         NORME = 'I'
         ITRNSE = 1
         INCE = LDE
      END IF

      // Check normalization of E

      ENRMIN = ONE / ULP
      ENRMAX = ZERO
      IF( ITRNSE.EQ.0 ) THEN

         // Eigenvectors are column vectors.

         IPAIR = 0
         DO 30 JVEC = 1, N
            TEMP1 = ZERO
            IF( IPAIR.EQ.0 .AND. JVEC.LT.N .AND. WI( JVEC ).NE.ZERO ) IPAIR = 1
            IF( IPAIR.EQ.1 ) THEN

               // Complex eigenvector

               DO 10 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) )+ ABS( E( J, JVEC+1 ) ) )
   10          CONTINUE
               ENRMIN = MIN( ENRMIN, TEMP1 )
               ENRMAX = MAX( ENRMAX, TEMP1 )
               IPAIR = 2
            ELSE IF( IPAIR.EQ.2 ) THEN
               IPAIR = 0
            ELSE

               // Real eigenvector

               DO 20 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) ) )
   20          CONTINUE
               ENRMIN = MIN( ENRMIN, TEMP1 )
               ENRMAX = MAX( ENRMAX, TEMP1 )
               IPAIR = 0
            END IF
   30    CONTINUE

      ELSE

         // Eigenvectors are row vectors.

         DO 40 JVEC = 1, N
            WORK( JVEC ) = ZERO
   40    CONTINUE

         DO 60 J = 1, N
            IPAIR = 0
            DO 50 JVEC = 1, N
               IF( IPAIR.EQ.0 .AND. JVEC.LT.N .AND. WI( JVEC ).NE.ZERO ) IPAIR = 1
               IF( IPAIR.EQ.1 ) THEN
                  WORK( JVEC ) = MAX( WORK( JVEC ), ABS( E( J, JVEC ) )+ABS( E( J, JVEC+1 ) ) )
                  WORK( JVEC+1 ) = WORK( JVEC )
               ELSE IF( IPAIR.EQ.2 ) THEN
                  IPAIR = 0
               ELSE
                  WORK( JVEC ) = MAX( WORK( JVEC ), ABS( E( J, JVEC ) ) )
                  IPAIR = 0
               END IF
   50       CONTINUE
   60    CONTINUE

         DO 70 JVEC = 1, N
            ENRMIN = MIN( ENRMIN, WORK( JVEC ) )
            ENRMAX = MAX( ENRMAX, WORK( JVEC ) )
   70    CONTINUE
      END IF

      // Norm of A:

      ANORM = MAX( SLANGE( NORMA, N, N, A, LDA, WORK ), UNFL )

      // Norm of E:

      ENORM = MAX( SLANGE( NORME, N, N, E, LDE, WORK ), ULP )

      // Norm of error:

      // Error =  AE - EW

      CALL SLASET( 'Full', N, N, ZERO, ZERO, WORK, N )

      IPAIR = 0
      IEROW = 1
      IECOL = 1

      DO 80 JCOL = 1, N
         IF( ITRNSE.EQ.1 ) THEN
            IEROW = JCOL
         ELSE
            IECOL = JCOL
         END IF

         IF( IPAIR.EQ.0 .AND. WI( JCOL ).NE.ZERO ) IPAIR = 1

         IF( IPAIR.EQ.1 ) THEN
            WMAT( 1, 1 ) = WR( JCOL )
            WMAT( 2, 1 ) = -WI( JCOL )
            WMAT( 1, 2 ) = WI( JCOL )
            WMAT( 2, 2 ) = WR( JCOL )
            CALL SGEMM( TRANSE, TRANSW, N, 2, 2, ONE, E( IEROW, IECOL ), LDE, WMAT, 2, ZERO, WORK( N*( JCOL-1 )+1 ), N )
            IPAIR = 2
         ELSE IF( IPAIR.EQ.2 ) THEN
            IPAIR = 0

         ELSE

            CALL SAXPY( N, WR( JCOL ), E( IEROW, IECOL ), INCE, WORK( N*( JCOL-1 )+1 ), 1 )
            IPAIR = 0
         END IF

   80 CONTINUE

      CALL SGEMM( TRANSA, TRANSE, N, N, N, ONE, A, LDA, E, LDE, -ONE, WORK, N )

      ERRNRM = SLANGE( 'One', N, N, WORK, N, WORK( N*N+1 ) ) / ENORM

      // Compute RESULT(1) (avoiding under/overflow)

      IF( ANORM.GT.ERRNRM ) THEN
         RESULT( 1 ) = ( ERRNRM / ANORM ) / ULP
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ONE / ULP
         ELSE
            RESULT( 1 ) = MIN( ERRNRM / ANORM, ONE ) / ULP
         END IF
      END IF

      // Compute RESULT(2) : the normalization error in E.

      RESULT( 2 ) = MAX( ABS( ENRMAX-ONE ), ABS( ENRMIN-ONE ) ) / ( REAL( N )*ULP )

      RETURN

      // End of SGET22

      END
