      SUBROUTINE ZLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA, BETA, WX, WY, S, DIF );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDX, LDY, N, TYPE;
      COMPLEX*16         ALPHA, BETA, WX, WY;
      // ..
      // .. Array Arguments ..
      double             DIF( * ), S( * );
      COMPLEX*16         A( LDA, * ), B( LDA, * ), X( LDX, * ), Y( LDY, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             RONE, TWO, THREE;
      const              RONE = 1.0, TWO = 2.0, THREE = 3.0 ;
      COMPLEX*16         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      double             RWORK( 50 );
      COMPLEX*16         WORK( 26 ), Z( 8, 8 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CDABS, DBLE, DCMPLX, DCONJG, SQRT
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGESVD, ZLACPY, ZLAKF2
      // ..
      // .. Executable Statements ..

      // Generate test problem ...
      // (Da, Db) ...

      for (I = 1; I <= N; I++) { // 20
         for (J = 1; J <= N; J++) { // 10

            if ( I == J ) {
               A( I, I ) = DCMPLX( I ) + ALPHA;
               B( I, I ) = ONE;
            } else {
               A( I, J ) = ZERO;
               B( I, J ) = ZERO;
            }

         } // 10
      } // 20
      if ( TYPE == 2 ) {
         A( 1, 1 ) = DCMPLX( RONE, RONE );
         A( 2, 2 ) = DCONJG( A( 1, 1 ) );
         A( 3, 3 ) = ONE;
         A( 4, 4 ) = DCMPLX( DBLE( ONE+ALPHA ), DBLE( ONE+BETA ) );
         A( 5, 5 ) = DCONJG( A( 4, 4 ) );
      }

      // Form X and Y

      zlacpy('F', N, N, B, LDA, Y, LDY );
      Y( 3, 1 ) = -DCONJG( WY );
      Y( 4, 1 ) = DCONJG( WY );
      Y( 5, 1 ) = -DCONJG( WY );
      Y( 3, 2 ) = -DCONJG( WY );
      Y( 4, 2 ) = DCONJG( WY );
      Y( 5, 2 ) = -DCONJG( WY );

      zlacpy('F', N, N, B, LDA, X, LDX );
      X( 1, 3 ) = -WX;
      X( 1, 4 ) = -WX;
      X( 1, 5 ) = WX;
      X( 2, 3 ) = WX;
      X( 2, 4 ) = -WX;
      X( 2, 5 ) = -WX;

      // Form (A, B)

      B( 1, 3 ) = WX + WY;
      B( 2, 3 ) = -WX + WY;
      B( 1, 4 ) = WX - WY;
      B( 2, 4 ) = WX - WY;
      B( 1, 5 ) = -WX + WY;
      B( 2, 5 ) = WX + WY;
      A( 1, 3 ) = WX*A( 1, 1 ) + WY*A( 3, 3 );
      A( 2, 3 ) = -WX*A( 2, 2 ) + WY*A( 3, 3 );
      A( 1, 4 ) = WX*A( 1, 1 ) - WY*A( 4, 4 );
      A( 2, 4 ) = WX*A( 2, 2 ) - WY*A( 4, 4 );
      A( 1, 5 ) = -WX*A( 1, 1 ) + WY*A( 5, 5 );
      A( 2, 5 ) = WX*A( 2, 2 ) + WY*A( 5, 5 );

      // Compute condition numbers

      S( 1 ) = RONE / SQRT( ( RONE+THREE*CDABS( WY )*CDABS( WY ) ) / ( RONE+CDABS( A( 1, 1 ) )*CDABS( A( 1, 1 ) ) ) )       S( 2 ) = RONE / SQRT( ( RONE+THREE*CDABS( WY )*CDABS( WY ) ) / ( RONE+CDABS( A( 2, 2 ) )*CDABS( A( 2, 2 ) ) ) )       S( 3 ) = RONE / SQRT( ( RONE+TWO*CDABS( WX )*CDABS( WX ) ) / ( RONE+CDABS( A( 3, 3 ) )*CDABS( A( 3, 3 ) ) ) )       S( 4 ) = RONE / SQRT( ( RONE+TWO*CDABS( WX )*CDABS( WX ) ) / ( RONE+CDABS( A( 4, 4 ) )*CDABS( A( 4, 4 ) ) ) )       S( 5 ) = RONE / SQRT( ( RONE+TWO*CDABS( WX )*CDABS( WX ) ) / ( RONE+CDABS( A( 5, 5 ) )*CDABS( A( 5, 5 ) ) ) );

      zlakf2(1, 4, A, LDA, A( 2, 2 ), B, B( 2, 2 ), Z, 8 );
      zgesvd('N', 'N', 8, 8, Z, 8, RWORK, WORK, 1, WORK( 2 ), 1, WORK( 3 ), 24, RWORK( 9 ), INFO );
      DIF( 1 ) = RWORK( 8 );

      zlakf2(4, 1, A, LDA, A( 5, 5 ), B, B( 5, 5 ), Z, 8 );
      zgesvd('N', 'N', 8, 8, Z, 8, RWORK, WORK, 1, WORK( 2 ), 1, WORK( 3 ), 24, RWORK( 9 ), INFO );
      DIF( 5 ) = RWORK( 8 );

      return;

      // End of ZLATM6

      }
