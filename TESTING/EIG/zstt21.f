      SUBROUTINE ZSTT21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KBAND, LDU, N;
      // ..
      // .. Array Arguments ..
      double             AD( * ), AE( * ), RESULT( 2 ), RWORK( * ), SD( * ), SE( * );
      COMPLEX*16         U( LDU, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANHE;
      // EXTERNAL DLAMCH, ZLANGE, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZHER, ZHER2, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // 1)      Constants

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 ) RETURN

      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Precision' )

      // Do Test 1

      // Copy A & Compute its 1-Norm:

      zlaset('Full', N, N, CZERO, CZERO, WORK, N );

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

      // Norm of A - USU*

      for (J = 1; J <= N; J++) { // 20
         zher('L', N, -SD( J ), U( 1, J ), 1, WORK, N );
   20 CONTINUE

      if ( N.GT.1 .AND. KBAND.EQ.1 ) {
         DO 30 J = 1, N - 1
            zher2('L', N, -DCMPLX( SE( J ) ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N );
   30    CONTINUE
      }

      WNORM = ZLANHE( '1', 'L', N, WORK, N, RWORK )

      if ( ANORM.GT.WNORM ) {
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      } else {
         if ( ANORM.LT.ONE ) {
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         } else {
            RESULT( 1 ) = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
         }
      }

      // Do Test 2

      // Compute  U U**H - I

      zgemm('N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N );

      for (J = 1; J <= N; J++) { // 40
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
   40 CONTINUE

      RESULT( 2 ) = MIN( DBLE( N ), ZLANGE( '1', N, N, WORK, N, RWORK ) ) / ( N*ULP )

      RETURN

      // End of ZSTT21

      }
