      void zlargv(N, X, INCX, Y, INCY, C, INCC ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCC, INCX, INCY, N;
      double             C( * );
      Complex         X( * ), Y( * );
      // ..

      double             TWO, ONE, ZERO;
      const              TWO = 2.0, ONE = 1.0, ZERO = 0.0 ;
      Complex         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // LOGICAL            FIRST

      int                COUNT, I, IC, IX, IY, J;
      double             CS, D, DI, DR, EPS, F2, F2S, G2, G2S, SAFMIN, SAFMN2, SAFMX2, SCALE;
      Complex         F, FF, FS, G, GS, R, SN;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLAPY2;
      // EXTERNAL DLAMCH, DLAPY2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, INT, LOG, MAX, SQRT
      // ..
      // .. Statement Functions ..
      double             ABS1, ABSSQ;
      // ..
      // .. Save statement ..
      // SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
      // ..
      // .. Data statements ..
      // DATA               FIRST / true /
      // ..
      // .. Statement Function definitions ..
      ABS1[FF] = max( ( FF.toDouble() ).abs(), ( DIMAG( FF ) ).abs() );
      ABSSQ[FF] = FF.toDouble()**2 + DIMAG( FF )**2;

      // IF( FIRST ) THEN
         // FIRST = false;
         SAFMIN = dlamch( 'S' );
         EPS = dlamch( 'E' );
         SAFMN2 = dlamch( 'B' )**INT( LOG( SAFMIN / EPS ) / LOG( dlamch( 'B' ) ) / TWO );
         SAFMX2 = ONE / SAFMN2;
      // END IF
      IX = 1;
      IY = 1;
      IC = 1;
      for (I = 1; I <= N; I++) { // 60
         F = X( IX );
         G = Y( IY );

         // Use identical algorithm as in ZLARTG

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
               R = dlapy2( G.toDouble(), DIMAG( G ) );
               // Do complex/real division explicitly with two real
               // divisions
               D = dlapy2( GS.toDouble(), DIMAG( GS ) );
               SN = DCMPLX( GS.toDouble() / D, -DIMAG( GS ) / D );
               GO TO 50;
            }
            F2S = dlapy2( FS.toDouble(), DIMAG( FS ) );
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
               D = dlapy2( F.toDouble(), DIMAG( F ) );
               FF = DCMPLX( F.toDouble() / D, DIMAG( F ) / D );
            } else {
               DR = SAFMX2*F.toDouble();
               DI = SAFMX2*DIMAG( F );
               D = dlapy2( DR, DI );
               FF = DCMPLX( DR / D, DI / D );
            }
            SN = FF*DCMPLX( GS.toDouble() / G2S, -DIMAG( GS ) / G2S );
            R = CS*F + SN*G;
         } else {

            // This is the most common case.
            // Neither F2 nor F2/G2 are less than SAFMIN
            // F2S cannot overflow, and it is accurate

            F2S = sqrt( ONE+G2 / F2 );
            // Do the F2S(real)*FS(complex) multiply with two real
            // multiplies
            R = DCMPLX( F2S*FS.toDouble(), F2S*DIMAG( FS ) );
            CS = ONE / F2S;
            D = F2 + G2;
            // Do complex/real division explicitly with two real divisions
            SN = DCMPLX( R.toDouble() / D, DIMAG( R ) / D );
            SN = SN*DCONJG( GS );
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
         C[IC] = CS;
         Y[IY] = SN;
         X[IX] = R;
         IC = IC + INCC;
         IY = IY + INCY;
         IX = IX + INCX;
      } // 60
      return;
      }
