      void clargv(N, X, INCX, Y, INCY, C, INCC ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCC, INCX, INCY, N;
      // ..
      // .. Array Arguments ..
      REAL               C( * );
      COMPLEX            X( * ), Y( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               TWO, ONE, ZERO;
      const              TWO = 2.0, ONE = 1.0, ZERO = 0.0 ;
      COMPLEX            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      // LOGICAL            FIRST
      int                COUNT, I, IC, IX, IY, J;
      REAL               CS, D, DI, DR, EPS, F2, F2S, G2, G2S, SAFMIN, SAFMN2, SAFMX2, SCALE;
      COMPLEX            F, FF, FS, G, GS, R, SN;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLAPY2;
      // EXTERNAL SLAMCH, SLAPY2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, INT, LOG, MAX, REAL, SQRT
      // ..
      // .. Statement Functions ..
      REAL               ABS1, ABSSQ;
      // ..
      // .. Save statement ..
      // SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
      // ..
      // .. Data statements ..
      // DATA               FIRST / true /
      // ..
      // .. Statement Function definitions ..
      ABS1( FF ) = max( ( REAL( FF ) ).abs(), ( AIMAG( FF ) ) ).abs();
      ABSSQ( FF ) = REAL( FF )**2 + AIMAG( FF )**2;
      // ..
      // .. Executable Statements ..

      // IF( FIRST ) THEN
         // FIRST = false;
         SAFMIN = SLAMCH( 'S' );
         EPS = SLAMCH( 'E' );
         SAFMN2 = SLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / LOG( SLAMCH( 'B' ) ) / TWO );
         SAFMX2 = ONE / SAFMN2;
      // END IF
      IX = 1;
      IY = 1;
      IC = 1;
      for (I = 1; I <= N; I++) { // 60
         F = X( IX );
         G = Y( IY );

         // Use identical algorithm as in CLARTG

         SCALE = max( ABS1( F ), ABS1( G ) );
         FS = F;
         GS = G;
         COUNT = 0;
         if ( SCALE >= SAFMX2 ) {
            } // 10
            COUNT = COUNT + 1;
            FS = FS*SAFMN2;
            GS = GS*SAFMN2;
            SCALE = SCALE*SAFMN2;
            if (SCALE >= SAFMX2 && COUNT < 20) GO TO 10;
         } else if ( SCALE <= SAFMN2 ) {
            if ( G == CZERO ) {
               CS = ONE;
               SN = CZERO;
               R = F;
               GO TO 50;
            }
            } // 20
            COUNT = COUNT - 1;
            FS = FS*SAFMX2;
            GS = GS*SAFMX2;
            SCALE = SCALE*SAFMX2;
            if (SCALE <= SAFMN2) GO TO 20;
         }
         F2 = ABSSQ( FS );
         G2 = ABSSQ( GS );
         if ( F2 <= max( G2, ONE )*SAFMIN ) {

            // This is a rare case: F is very small.

            if ( F == CZERO ) {
               CS = ZERO;
               R = SLAPY2( REAL( G ), AIMAG( G ) );
               // Do complex/real division explicitly with two real
               // divisions
               D = SLAPY2( REAL( GS ), AIMAG( GS ) );
               SN = CMPLX( REAL( GS ) / D, -AIMAG( GS ) / D );
               GO TO 50;
            }
            F2S = SLAPY2( REAL( FS ), AIMAG( FS ) );
            // G2 and G2S are accurate
            // G2 is at least SAFMIN, and G2S is at least SAFMN2
            G2S = sqrt( G2 );
            // Error in CS from underflow in F2S is at most
            // UNFL / SAFMN2 < sqrt(UNFL*EPS) < EPS
            // If max(G2,ONE)=G2, then F2 < G2*SAFMIN,
            // and so CS < sqrt(SAFMIN)
            // If max(G2,ONE)=ONE, then F2 < SAFMIN
            // and so CS < sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
            // Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
            CS = F2S / G2S;
            // Make sure abs(FF) = 1
            // Do complex/real division explicitly with 2 real divisions
            if ( ABS1( F ) > ONE ) {
               D = SLAPY2( REAL( F ), AIMAG( F ) );
               FF = CMPLX( REAL( F ) / D, AIMAG( F ) / D );
            } else {
               DR = SAFMX2*REAL( F );
               DI = SAFMX2*AIMAG( F );
               D = SLAPY2( DR, DI );
               FF = CMPLX( DR / D, DI / D );
            }
            SN = FF*CMPLX( REAL( GS ) / G2S, -AIMAG( GS ) / G2S );
            R = CS*F + SN*G;
         } else {

            // This is the most common case.
            // Neither F2 nor F2/G2 are less than SAFMIN
            // F2S cannot overflow, and it is accurate

            F2S = sqrt( ONE+G2 / F2 );
            // Do the F2S(real)*FS(complex) multiply with two real
            // multiplies
            R = CMPLX( F2S*REAL( FS ), F2S*AIMAG( FS ) );
            CS = ONE / F2S;
            D = F2 + G2;
            // Do complex/real division explicitly with two real divisions
            SN = CMPLX( REAL( R ) / D, AIMAG( R ) / D );
            SN = SN*CONJG( GS );
            if ( COUNT != 0 ) {
               if ( COUNT > 0 ) {
                  for (J = 1; J <= COUNT; J++) { // 30
                     R = R*SAFMX2;
                  } // 30
               } else {
                  for (J = 1; J <= -COUNT; J++) { // 40
                     R = R*SAFMN2;
                  } // 40
               }
            }
         }
         } // 50
         C( IC ) = CS;
         Y( IY ) = SN;
         X( IX ) = R;
         IC = IC + INCC;
         IY = IY + INCY;
         IX = IX + INCX;
      } // 60
      return;
      }
