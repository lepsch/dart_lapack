      SUBROUTINE CGET22( TRANSA, TRANSE, TRANSW, N, A, LDA, E, LDE, W, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSA, TRANSE, TRANSW;
      int                LDA, LDE, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), E( LDE, * ), W( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      String             NORMA, NORME;
      int                ITRNSE, ITRNSW, J, JCOL, JOFF, JROW, JVEC;
      REAL               ANORM, ENORM, ENRMAX, ENRMIN, ERRNRM, TEMP1, ULP, UNFL
      COMPLEX            WTEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH
      // EXTERNAL LSAME, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Initialize RESULT (in case N=0)

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 ) RETURN

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )

      ITRNSE = 0
      ITRNSW = 0
      NORMA = 'O'
      NORME = 'O'

      IF( LSAME( TRANSA, 'T' ) .OR. LSAME( TRANSA, 'C' ) ) THEN
         NORMA = 'I'
      END IF

      IF( LSAME( TRANSE, 'T' ) ) THEN
         ITRNSE = 1
         NORME = 'I'
      ELSE IF( LSAME( TRANSE, 'C' ) ) THEN
         ITRNSE = 2
         NORME = 'I'
      END IF

      IF( LSAME( TRANSW, 'C' ) ) THEN
         ITRNSW = 1
      END IF

      // Normalization of E:

      ENRMIN = ONE / ULP
      ENRMAX = ZERO
      IF( ITRNSE.EQ.0 ) THEN
         DO 20 JVEC = 1, N
            TEMP1 = ZERO
            DO 10 J = 1, N
               TEMP1 = MAX( TEMP1, ABS( REAL( E( J, JVEC ) ) )+ ABS( AIMAG( E( J, JVEC ) ) ) )
   10       CONTINUE
            ENRMIN = MIN( ENRMIN, TEMP1 )
            ENRMAX = MAX( ENRMAX, TEMP1 )
   20    CONTINUE
      ELSE
         DO 30 JVEC = 1, N
            RWORK( JVEC ) = ZERO
   30    CONTINUE

         DO 50 J = 1, N
            DO 40 JVEC = 1, N
               RWORK( JVEC ) = MAX( RWORK( JVEC ), ABS( REAL( E( JVEC, J ) ) )+ ABS( AIMAG( E( JVEC, J ) ) ) )
   40       CONTINUE
   50    CONTINUE

         DO 60 JVEC = 1, N
            ENRMIN = MIN( ENRMIN, RWORK( JVEC ) )
            ENRMAX = MAX( ENRMAX, RWORK( JVEC ) )
   60    CONTINUE
      END IF

      // Norm of A:

      ANORM = MAX( CLANGE( NORMA, N, N, A, LDA, RWORK ), UNFL )

      // Norm of E:

      ENORM = MAX( CLANGE( NORME, N, N, E, LDE, RWORK ), ULP )

      // Norm of error:

      // Error =  AE - EW

      CALL CLASET( 'Full', N, N, CZERO, CZERO, WORK, N )

      JOFF = 0
      DO 100 JCOL = 1, N
         IF( ITRNSW.EQ.0 ) THEN
            WTEMP = W( JCOL )
         ELSE
            WTEMP = CONJG( W( JCOL ) )
         END IF

         IF( ITRNSE.EQ.0 ) THEN
            DO 70 JROW = 1, N
               WORK( JOFF+JROW ) = E( JROW, JCOL )*WTEMP
   70       CONTINUE
         ELSE IF( ITRNSE.EQ.1 ) THEN
            DO 80 JROW = 1, N
               WORK( JOFF+JROW ) = E( JCOL, JROW )*WTEMP
   80       CONTINUE
         ELSE
            DO 90 JROW = 1, N
               WORK( JOFF+JROW ) = CONJG( E( JCOL, JROW ) )*WTEMP
   90       CONTINUE
         END IF
         JOFF = JOFF + N
  100 CONTINUE

      CALL CGEMM( TRANSA, TRANSE, N, N, N, CONE, A, LDA, E, LDE, -CONE, WORK, N )

      ERRNRM = CLANGE( 'One', N, N, WORK, N, RWORK ) / ENORM

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

      // End of CGET22

      END
