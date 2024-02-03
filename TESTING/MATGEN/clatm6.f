      SUBROUTINE CLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA, BETA, WX, WY, S, DIF )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDX, LDY, N, TYPE;
      COMPLEX            ALPHA, BETA, WX, WY
      // ..
      // .. Array Arguments ..
      REAL               DIF( * ), S( * )
      COMPLEX            A( LDA, * ), B( LDA, * ), X( LDX, * ), Y( LDY, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               RONE, TWO, THREE
      const              RONE = 1.0E+0, TWO = 2.0E+0, THREE = 3.0E+0 ;
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 50 )
      COMPLEX            WORK( 26 ), Z( 8, 8 )
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CABS, CMPLX, CONJG, REAL, SQRT
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGESVD, CLACPY, CLAKF2
      // ..
      // .. Executable Statements ..

      // Generate test problem ...
      // (Da, Db) ...

      for (I = 1; I <= N; I++) { // 20
         for (J = 1; J <= N; J++) { // 10

            if ( I.EQ.J ) {
               A( I, I ) = CMPLX( I ) + ALPHA
               B( I, I ) = ONE
            } else {
               A( I, J ) = ZERO
               B( I, J ) = ZERO
            }

         } // 10
      } // 20
      if ( TYPE.EQ.2 ) {
         A( 1, 1 ) = CMPLX( RONE, RONE )
         A( 2, 2 ) = CONJG( A( 1, 1 ) )
         A( 3, 3 ) = ONE
         A( 4, 4 ) = CMPLX( REAL( ONE+ALPHA ), REAL( ONE+BETA ) )
         A( 5, 5 ) = CONJG( A( 4, 4 ) )
      }

      // Form X and Y

      clacpy('F', N, N, B, LDA, Y, LDY );
      Y( 3, 1 ) = -CONJG( WY )
      Y( 4, 1 ) = CONJG( WY )
      Y( 5, 1 ) = -CONJG( WY )
      Y( 3, 2 ) = -CONJG( WY )
      Y( 4, 2 ) = CONJG( WY )
      Y( 5, 2 ) = -CONJG( WY )

      clacpy('F', N, N, B, LDA, X, LDX );
      X( 1, 3 ) = -WX
      X( 1, 4 ) = -WX
      X( 1, 5 ) = WX
      X( 2, 3 ) = WX
      X( 2, 4 ) = -WX
      X( 2, 5 ) = -WX

      // Form (A, B)

      B( 1, 3 ) = WX + WY
      B( 2, 3 ) = -WX + WY
      B( 1, 4 ) = WX - WY
      B( 2, 4 ) = WX - WY
      B( 1, 5 ) = -WX + WY
      B( 2, 5 ) = WX + WY
      A( 1, 3 ) = WX*A( 1, 1 ) + WY*A( 3, 3 )
      A( 2, 3 ) = -WX*A( 2, 2 ) + WY*A( 3, 3 )
      A( 1, 4 ) = WX*A( 1, 1 ) - WY*A( 4, 4 )
      A( 2, 4 ) = WX*A( 2, 2 ) - WY*A( 4, 4 )
      A( 1, 5 ) = -WX*A( 1, 1 ) + WY*A( 5, 5 )
      A( 2, 5 ) = WX*A( 2, 2 ) + WY*A( 5, 5 )

      // Compute condition numbers

      S( 1 ) = RONE / SQRT( ( RONE+THREE*CABS( WY )*CABS( WY ) ) / ( RONE+CABS( A( 1, 1 ) )*CABS( A( 1, 1 ) ) ) )       S( 2 ) = RONE / SQRT( ( RONE+THREE*CABS( WY )*CABS( WY ) ) / ( RONE+CABS( A( 2, 2 ) )*CABS( A( 2, 2 ) ) ) )       S( 3 ) = RONE / SQRT( ( RONE+TWO*CABS( WX )*CABS( WX ) ) / ( RONE+CABS( A( 3, 3 ) )*CABS( A( 3, 3 ) ) ) )       S( 4 ) = RONE / SQRT( ( RONE+TWO*CABS( WX )*CABS( WX ) ) / ( RONE+CABS( A( 4, 4 ) )*CABS( A( 4, 4 ) ) ) )       S( 5 ) = RONE / SQRT( ( RONE+TWO*CABS( WX )*CABS( WX ) ) / ( RONE+CABS( A( 5, 5 ) )*CABS( A( 5, 5 ) ) ) )

      clakf2(1, 4, A, LDA, A( 2, 2 ), B, B( 2, 2 ), Z, 8 );
      cgesvd('N', 'N', 8, 8, Z, 8, RWORK, WORK, 1, WORK( 2 ), 1, WORK( 3 ), 24, RWORK( 9 ), INFO );
      DIF( 1 ) = RWORK( 8 )

      clakf2(4, 1, A, LDA, A( 5, 5 ), B, B( 5, 5 ), Z, 8 );
      cgesvd('N', 'N', 8, 8, Z, 8, RWORK, WORK, 1, WORK( 2 ), 1, WORK( 3 ), 24, RWORK( 9 ), INFO );
      DIF( 5 ) = RWORK( 8 )

      RETURN

      // End of CLATM6

      }
