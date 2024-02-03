      SUBROUTINE SSTT22( N, M, KBAND, AD, AE, SD, SE, U, LDU, WORK, LDWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KBAND, LDU, LDWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               AD( * ), AE( * ), RESULT( 2 ), SD( * ), SE( * ), U( LDU, * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      REAL               ANORM, AUKJ, ULP, UNFL, WNORM
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 .OR. M.LE.0 ) RETURN

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )

      // Do Test 1

      // Compute the 1-norm of A.

      IF( N.GT.1 ) THEN
         ANORM = ABS( AD( 1 ) ) + ABS( AE( 1 ) )
         DO 10 J = 2, N - 1
            ANORM = MAX( ANORM, ABS( AD( J ) )+ABS( AE( J ) )+ ABS( AE( J-1 ) ) )
   10    CONTINUE
         ANORM = MAX( ANORM, ABS( AD( N ) )+ABS( AE( N-1 ) ) )
      ELSE
         ANORM = ABS( AD( 1 ) )
      END IF
      ANORM = MAX( ANORM, UNFL )

      // Norm of U'AU - S

      DO 40 I = 1, M
         DO 30 J = 1, M
            WORK( I, J ) = ZERO
            DO 20 K = 1, N
               AUKJ = AD( K )*U( K, J )
               IF( K.NE.N ) AUKJ = AUKJ + AE( K )*U( K+1, J )                IF( K.NE.1 ) AUKJ = AUKJ + AE( K-1 )*U( K-1, J )
               WORK( I, J ) = WORK( I, J ) + U( K, I )*AUKJ
   20       CONTINUE
   30    CONTINUE
         WORK( I, I ) = WORK( I, I ) - SD( I )
         IF( KBAND.EQ.1 ) THEN
            IF( I.NE.1 ) WORK( I, I-1 ) = WORK( I, I-1 ) - SE( I-1 )             IF( I.NE.N ) WORK( I, I+1 ) = WORK( I, I+1 ) - SE( I )
         END IF
   40 CONTINUE

      WNORM = SLANSY( '1', 'L', M, WORK, M, WORK( 1, M+1 ) )

      IF( ANORM.GT.WNORM ) THEN
         RESULT( 1 ) = ( WNORM / ANORM ) / ( M*ULP )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*ULP )
         ELSE
            RESULT( 1 ) = MIN( WNORM / ANORM, REAL( M ) ) / ( M*ULP )
         END IF
      END IF

      // Do Test 2

      // Compute  U'U - I

      CALL SGEMM( 'T', 'N', M, M, N, ONE, U, LDU, U, LDU, ZERO, WORK, M )

      DO 50 J = 1, M
         WORK( J, J ) = WORK( J, J ) - ONE
   50 CONTINUE

      RESULT( 2 ) = MIN( REAL( M ), SLANGE( '1', M, M, WORK, M, WORK( 1, M+1 ) ) ) / ( M*ULP )

      RETURN

      // End of SSTT22

      END
