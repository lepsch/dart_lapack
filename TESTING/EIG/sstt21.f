      SUBROUTINE SSTT21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KBAND, LDU, N;
      // ..
      // .. Array Arguments ..
      REAL               AD( * ), AE( * ), RESULT( 2 ), SD( * ), SE( * ), U( LDU, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLASET, SSYR, SSYR2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // 1)      Constants

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 ) RETURN

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )

      // Do Test 1

      // Copy A & Compute its 1-Norm:

      CALL SLASET( 'Full', N, N, ZERO, ZERO, WORK, N )

      ANORM = ZERO
      TEMP1 = ZERO

      DO 10 J = 1, N - 1
         WORK( ( N+1 )*( J-1 )+1 ) = AD( J )
         WORK( ( N+1 )*( J-1 )+2 ) = AE( J )
         TEMP2 = ABS( AE( J ) )
         ANORM = MAX( ANORM, ABS( AD( J ) )+TEMP1+TEMP2 )
         TEMP1 = TEMP2
   10 CONTINUE

      WORK( N**2 ) = AD( N )
      ANORM = MAX( ANORM, ABS( AD( N ) )+TEMP1, UNFL )

      // Norm of A - USU'

      DO 20 J = 1, N
         CALL SSYR( 'L', N, -SD( J ), U( 1, J ), 1, WORK, N )
   20 CONTINUE

      IF( N.GT.1 .AND. KBAND.EQ.1 ) THEN
         DO 30 J = 1, N - 1
            CALL SSYR2( 'L', N, -SE( J ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N )
   30    CONTINUE
      END IF

      WNORM = SLANSY( '1', 'L', N, WORK, N, WORK( N**2+1 ) )

      IF( ANORM.GT.WNORM ) THEN
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         ELSE
            RESULT( 1 ) = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
         END IF
      END IF

      // Do Test 2

      // Compute  UU' - I

      CALL SGEMM( 'N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N )

      DO 40 J = 1, N
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
   40 CONTINUE

      RESULT( 2 ) = MIN( REAL( N ), SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ) ) / ( N*ULP )

      RETURN

      // End of SSTT21

      END
