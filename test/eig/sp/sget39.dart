      void sget39(RMAX, LMAX, NINFO, KNT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX, NINFO;
      double               RMAX;
      // ..

// =====================================================================

      // .. Parameters ..
      int                LDT, LDT2;
      const              LDT = 10, LDT2 = 2*LDT ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, IVM1, IVM2, IVM3, IVM4, IVM5, J, K, N, NDIM;
      double               BIGNUM, DOMIN, DUMM, EPS, NORM, NORMTB, RESID, SCALE, SMLNUM, W, XNORM;
      // ..
      // .. External Functions ..
      //- int                ISAMAX;
      //- REAL               SASUM, SDOT, SLAMCH, SLANGE;
      // EXTERNAL ISAMAX, SASUM, SDOT, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMV, SLAQTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, COS, MAX, REAL, SIN, SQRT
      // ..
      // .. Local Arrays ..
      int                IDIM( 6 ), IVAL( 5, 5, 6 );
      double               B( LDT ), D( LDT2 ), DUM( 1 ), T( LDT, LDT ), VM1( 5 ), VM2( 5 ), VM3( 5 ), VM4( 5 ), VM5( 3 ), WORK( LDT ), X( LDT2 ), Y( LDT2 );
      // ..
      // .. Data statements ..
      const IDIM = [ 4, 5, 5, 5, 5, 5,];
      const IVAL = [ 3, 0, 0, 0, 0, 1, 1, -1, 0, 0, 3, 2, 1, 0, 0, 4, 3, 2, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 3, 3, 4, 0, 0, 4, 2, 2, 3, 0, 1, 1, 1, 1, 5, 1, 0, 0, 0, 0, 2, 4, -2, 0, 0, 3, 3, 4, 0, 0, 4, 2, 2, 3, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 2, 1, -1, 0, 0, 9, 8, 1, 0, 0, 4, 9, 1, 2, -1, 2, 2, 2, 2, 2, 9, 0, 0, 0, 0, 6, 4, 0, 0, 0, 3, 2, 1, 1, 0, 5, 1, -1, 1, 0, 2, 2, 2, 2, 2, 4, 0, 0, 0, 0, 2, 2, 0, 0, 0, 1, 4, 4, 0, 0, 2, 4, 2, 2, -1, 2, 2, 2, 2, 2,];
      // ..
      // .. Executable Statements ..

      // Get machine parameters

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;

      // Set up test case parameters

      VM1[1] = ONE;
      VM1[2] = sqrt( SMLNUM );
      VM1[3] = sqrt( VM1( 2 ) );
      VM1[4] = sqrt( BIGNUM );
      VM1[5] = sqrt( VM1( 4 ) );

      VM2[1] = ONE;
      VM2[2] = sqrt( SMLNUM );
      VM2[3] = sqrt( VM2( 2 ) );
      VM2[4] = sqrt( BIGNUM );
      VM2[5] = sqrt( VM2( 4 ) );

      VM3[1] = ONE;
      VM3[2] = sqrt( SMLNUM );
      VM3[3] = sqrt( VM3( 2 ) );
      VM3[4] = sqrt( BIGNUM );
      VM3[5] = sqrt( VM3( 4 ) );

      VM4[1] = ONE;
      VM4[2] = sqrt( SMLNUM );
      VM4[3] = sqrt( VM4( 2 ) );
      VM4[4] = sqrt( BIGNUM );
      VM4[5] = sqrt( VM4( 4 ) );

      VM5[1] = ONE;
      VM5[2] = EPS;
      VM5[3] = sqrt( SMLNUM );

      // Initialization

      KNT = 0;
      RMAX = ZERO;
      NINFO = 0;
      SMLNUM = SMLNUM / EPS;

      // Begin test loop

      for (IVM5 = 1; IVM5 <= 3; IVM5++) { // 140
         for (IVM4 = 1; IVM4 <= 5; IVM4++) { // 130
            for (IVM3 = 1; IVM3 <= 5; IVM3++) { // 120
               for (IVM2 = 1; IVM2 <= 5; IVM2++) { // 110
                  for (IVM1 = 1; IVM1 <= 5; IVM1++) { // 100
                     for (NDIM = 1; NDIM <= 6; NDIM++) { // 90

                        N = IDIM( NDIM );
                        for (I = 1; I <= N; I++) { // 20
                           for (J = 1; J <= N; J++) { // 10
                              T[I][J] = double( IVAL( I, J, NDIM ) )* VM1( IVM1 )                               IF( I >= J ) T( I, J ) = T( I, J )*VM5( IVM5 );
                           } // 10
                        } // 20

                        W = ONE*VM2( IVM2 );

                        for (I = 1; I <= N; I++) { // 30
                           B[I] = COS( double( I ) )*VM3( IVM3 );
                        } // 30

                        for (I = 1; I <= 2*N; I++) { // 40
                           D[I] = SIN( double( I ) )*VM4( IVM4 );
                        } // 40

                        NORM = SLANGE( '1', N, N, T, LDT, WORK );
                        K = ISAMAX( N, B, 1 );
                        NORMTB = NORM + ( B( K ) ).abs() + ( W ).abs();

                        scopy(N, D, 1, X, 1 );
                        KNT = KNT + 1;
                        CALL SLAQTR( false , true , N, T, LDT, DUM, DUMM, SCALE, X, WORK, INFO )                         IF( INFO != 0 ) NINFO = NINFO + 1;

                        // || T*x - scale*d || /
                          // max(ulp*||T||*||x||,smlnum/ulp*||T||,smlnum)

                        scopy(N, D, 1, Y, 1 );
                        sgemv('No transpose', N, N, ONE, T, LDT, X, 1, -SCALE, Y, 1 );
                        XNORM = SASUM( N, X, 1 );
                        RESID = SASUM( N, Y, 1 );
                        DOMIN = max( SMLNUM, ( SMLNUM / EPS )*NORM, ( NORM*EPS )*XNORM );
                        RESID = RESID / DOMIN;
                        if ( RESID > RMAX ) {
                           RMAX = RESID;
                           LMAX = KNT;
                        }

                        scopy(N, D, 1, X, 1 );
                        KNT = KNT + 1;
                        CALL SLAQTR( true , true , N, T, LDT, DUM, DUMM, SCALE, X, WORK, INFO )                         IF( INFO != 0 ) NINFO = NINFO + 1;

                        // || T*x - scale*d || /
                          // max(ulp*||T||*||x||,smlnum/ulp*||T||,smlnum)

                        scopy(N, D, 1, Y, 1 );
                        sgemv('Transpose', N, N, ONE, T, LDT, X, 1, -SCALE, Y, 1 );
                        XNORM = SASUM( N, X, 1 );
                        RESID = SASUM( N, Y, 1 );
                        DOMIN = max( SMLNUM, ( SMLNUM / EPS )*NORM, ( NORM*EPS )*XNORM );
                        RESID = RESID / DOMIN;
                        if ( RESID > RMAX ) {
                           RMAX = RESID;
                           LMAX = KNT;
                        }

                        scopy(2*N, D, 1, X, 1 );
                        KNT = KNT + 1;
                        CALL SLAQTR( false , false , N, T, LDT, B, W, SCALE, X, WORK, INFO )                         IF( INFO != 0 ) NINFO = NINFO + 1;

                        // ||(T+i*B)*(x1+i*x2) - scale*(d1+i*d2)|| /
                           // max(ulp*(||T||+||B||)*(||x1||+||x2||),
                                   // smlnum/ulp * (||T||+||B||), smlnum )


                        scopy(2*N, D, 1, Y, 1 );
                        Y[1] = SDOT( N, B, 1, X( 1+N ), 1 ) + SCALE*Y( 1 );
                        for (I = 2; I <= N; I++) { // 50
                           Y[I] = W*X( I+N ) + SCALE*Y( I );
                        } // 50
                        sgemv('No transpose', N, N, ONE, T, LDT, X, 1, -ONE, Y, 1 );

                        Y[1+N] = SDOT( N, B, 1, X, 1 ) - SCALE*Y( 1+N );
                        for (I = 2; I <= N; I++) { // 60
                           Y[I+N] = W*X( I ) - SCALE*Y( I+N );
                        } // 60
                        sgemv('No transpose', N, N, ONE, T, LDT, X( 1+N ), 1, ONE, Y( 1+N ), 1 );

                        RESID = SASUM( 2*N, Y, 1 );
                        DOMIN = max( SMLNUM, ( SMLNUM / EPS )*NORMTB, EPS*( NORMTB*SASUM( 2*N, X, 1 ) ) );
                        RESID = RESID / DOMIN;
                        if ( RESID > RMAX ) {
                           RMAX = RESID;
                           LMAX = KNT;
                        }

                        scopy(2*N, D, 1, X, 1 );
                        KNT = KNT + 1;
                        CALL SLAQTR( true , false , N, T, LDT, B, W, SCALE, X, WORK, INFO )                         IF( INFO != 0 ) NINFO = NINFO + 1;

                        // ||(T+i*B)*(x1+i*x2) - scale*(d1+i*d2)|| /
                           // max(ulp*(||T||+||B||)*(||x1||+||x2||),
                                   // smlnum/ulp * (||T||+||B||), smlnum )

                        scopy(2*N, D, 1, Y, 1 );
                        Y[1] = B( 1 )*X( 1+N ) - SCALE*Y( 1 );
                        for (I = 2; I <= N; I++) { // 70
                           Y[I] = B( I )*X( 1+N ) + W*X( I+N ) - SCALE*Y( I );
                        } // 70
                        sgemv('Transpose', N, N, ONE, T, LDT, X, 1, ONE, Y, 1 );

                        Y[1+N] = B( 1 )*X( 1 ) + SCALE*Y( 1+N );
                        for (I = 2; I <= N; I++) { // 80
                           Y[I+N] = B( I )*X( 1 ) + W*X( I ) + SCALE*Y( I+N );
                        } // 80
                        sgemv('Transpose', N, N, ONE, T, LDT, X( 1+N ), 1, -ONE, Y( 1+N ), 1 );

                        RESID = SASUM( 2*N, Y, 1 );
                        DOMIN = max( SMLNUM, ( SMLNUM / EPS )*NORMTB, EPS*( NORMTB*SASUM( 2*N, X, 1 ) ) );
                        RESID = RESID / DOMIN;
                        if ( RESID > RMAX ) {
                           RMAX = RESID;
                           LMAX = KNT;
                        }

                     } // 90
                  } // 100
               } // 110
            } // 120
         } // 130
      } // 140

      return;
      }
