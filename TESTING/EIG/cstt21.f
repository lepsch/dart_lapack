      SUBROUTINE CSTT21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KBAND, LDU, N;
      // ..
      // .. Array Arguments ..
      REAL               AD( * ), AE( * ), RESULT( 2 ), RWORK( * ), SD( * ), SE( * )
      COMPLEX            U( LDU, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM
      // ..
      // .. External Functions ..
      REAL               CLANGE, CLANHE, SLAMCH
      // EXTERNAL CLANGE, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHER, CHER2, CLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX, MIN, REAL
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

      CALL CLASET( 'Full', N, N, CZERO, CZERO, WORK, N )

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

      // Norm of A - U S U**H

      DO 20 J = 1, N
         CALL CHER( 'L', N, -SD( J ), U( 1, J ), 1, WORK, N )
   20 CONTINUE

      IF( N.GT.1 .AND. KBAND.EQ.1 ) THEN
         DO 30 J = 1, N - 1
            CALL CHER2( 'L', N, -CMPLX( SE( J ) ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N )
   30    CONTINUE
      END IF

      WNORM = CLANHE( '1', 'L', N, WORK, N, RWORK )

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

      // Compute  U U**H - I

      CALL CGEMM( 'N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N )

      DO 40 J = 1, N
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
   40 CONTINUE

      RESULT( 2 ) = MIN( REAL( N ), CLANGE( '1', N, N, WORK, N, RWORK ) ) / ( N*ULP )

      RETURN

      // End of CSTT21

      }
