      void dget31(RMAX, LMAX, NINFO, KNT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX;
      double             RMAX;
      // ..
      // .. Array Arguments ..
      int                NINFO( 2 );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      double             TWO, THREE, FOUR;
      const              TWO = 2.0, THREE = 3.0, FOUR = 4.0 ;
      double             SEVEN, TEN;
      const              SEVEN = 7.0, TEN = 10.0 ;
      double             TWNONE;
      const              TWNONE = 21.0 ;
      // ..
      // .. Local Scalars ..
      int                IA, IB, ICA, ID1, ID2, INFO, ISMIN, ITRANS, IWI, IWR, NA, NW;
      double             BIGNUM, CA, D1, D2, DEN, EPS, RES, SCALE, SMIN, SMLNUM, TMP, UNFL, WI, WR, XNORM;
      // ..
      // .. Local Arrays ..
      bool               LTRANS( 0: 1 );
      double             A( 2, 2 ), B( 2, 2 ), VAB( 3 ), VCA( 5 ), VDD( 4 ), VSMIN( 4 ), VWI( 4 ), VWR( 4 ), X( 2, 2 );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLALN2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Data statements ..
      const LTRANS = [ false , true ];
      // ..
      // .. Executable Statements ..

      // Get machine parameters

      EPS = DLAMCH( 'P' );
      UNFL = DLAMCH( 'U' );
      SMLNUM = DLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Set up test case parameters

      VSMIN( 1 ) = SMLNUM;
      VSMIN( 2 ) = EPS;
      VSMIN( 3 ) = ONE / ( TEN*TEN );
      VSMIN( 4 ) = ONE / EPS;
      VAB( 1 ) = sqrt( SMLNUM );
      VAB( 2 ) = ONE;
      VAB( 3 ) = sqrt( BIGNUM );
      VWR( 1 ) = ZERO;
      VWR( 2 ) = HALF;
      VWR( 3 ) = TWO;
      VWR( 4 ) = ONE;
      VWI( 1 ) = SMLNUM;
      VWI( 2 ) = EPS;
      VWI( 3 ) = ONE;
      VWI( 4 ) = TWO;
      VDD( 1 ) = sqrt( SMLNUM );
      VDD( 2 ) = ONE;
      VDD( 3 ) = TWO;
      VDD( 4 ) = sqrt( BIGNUM );
      VCA( 1 ) = ZERO;
      VCA( 2 ) = sqrt( SMLNUM );
      VCA( 3 ) = EPS;
      VCA( 4 ) = HALF;
      VCA( 5 ) = ONE;

      KNT = 0;
      NINFO( 1 ) = 0;
      NINFO( 2 ) = 0;
      LMAX = 0;
      RMAX = ZERO;

      // Begin test loop

      for (ID1 = 1; ID1 <= 4; ID1++) { // 190
         D1 = VDD( ID1 );
         for (ID2 = 1; ID2 <= 4; ID2++) { // 180
            D2 = VDD( ID2 );
            for (ICA = 1; ICA <= 5; ICA++) { // 170
               CA = VCA( ICA );
               for (ITRANS = 0; ITRANS <= 1; ITRANS++) { // 160
                  for (ISMIN = 1; ISMIN <= 4; ISMIN++) { // 150
                     SMIN = VSMIN( ISMIN );

                     NA = 1;
                     NW = 1;
                     for (IA = 1; IA <= 3; IA++) { // 30
                        A( 1, 1 ) = VAB( IA );
                        for (IB = 1; IB <= 3; IB++) { // 20
                           B( 1, 1 ) = VAB( IB );
                           for (IWR = 1; IWR <= 4; IWR++) { // 10
                              if ( D1 == ONE && D2 == ONE && CA == ONE ) {
                                 WR = VWR( IWR )*A( 1, 1 );
                              } else {
                                 WR = VWR( IWR );
                              }
                              WI = ZERO;
                              dlaln2(LTRANS( ITRANS ), NA, NW, SMIN, CA, A, 2, D1, D2, B, 2, WR, WI, X, 2, SCALE, XNORM, INFO );
                              if (INFO < 0) NINFO( 1 ) = NINFO( 1 ) + 1;
                              if( INFO > 0 ) NINFO( 2 ) = NINFO( 2 ) + 1;
                              RES = ABS( ( CA*A( 1, 1 )-WR*D1 )* X( 1, 1 )-SCALE*B( 1, 1 ) );
                              if ( INFO == 0 ) {
                                 DEN = max( EPS*( ABS( ( CA*A( 1, 1 )-WR*D1 )*X( 1, 1 ) ) ), SMLNUM );
                              } else {
                                 DEN = max( SMIN*ABS( X( 1, 1 ) ), SMLNUM );
                              }
                              RES = RES / DEN;
                              if( ABS( X( 1, 1 ) ) < UNFL && ABS( B( 1, 1 ) ) <= SMLNUM* ABS( CA*A( 1, 1 )-WR*D1 ) )RES = ZERO;
                              if (SCALE > ONE) RES = RES + ONE / EPS;
                              RES = RES + ABS( XNORM-ABS( X( 1, 1 ) ) ) / max( SMLNUM, XNORM ) / EPS                               IF( INFO != 0 && INFO != 1 ) RES = RES + ONE / EPS;
                              KNT = KNT + 1;
                              if ( RES > RMAX ) {
                                 LMAX = KNT;
                                 RMAX = RES;
                              }
                           } // 10
                        } // 20
                     } // 30

                     NA = 1;
                     NW = 2;
                     for (IA = 1; IA <= 3; IA++) { // 70
                        A( 1, 1 ) = VAB( IA );
                        for (IB = 1; IB <= 3; IB++) { // 60
                           B( 1, 1 ) = VAB( IB );
                           B( 1, 2 ) = -HALF*VAB( IB );
                           for (IWR = 1; IWR <= 4; IWR++) { // 50
                              if ( D1 == ONE && D2 == ONE && CA == ONE ) {
                                 WR = VWR( IWR )*A( 1, 1 );
                              } else {
                                 WR = VWR( IWR );
                              }
                              for (IWI = 1; IWI <= 4; IWI++) { // 40
                                 if ( D1 == ONE && D2 == ONE && CA == ONE ) {
                                    WI = VWI( IWI )*A( 1, 1 );
                                 } else {
                                    WI = VWI( IWI );
                                 }
                                 dlaln2(LTRANS( ITRANS ), NA, NW, SMIN, CA, A, 2, D1, D2, B, 2, WR, WI, X, 2, SCALE, XNORM, INFO );
                                 if (INFO < 0) NINFO( 1 ) = NINFO( 1 ) + 1;
                                 if( INFO > 0 ) NINFO( 2 ) = NINFO( 2 ) + 1;
                                 RES = ABS( ( CA*A( 1, 1 )-WR*D1 )* X( 1, 1 )+( WI*D1 )*X( 1, 2 )- SCALE*B( 1, 1 ) )                                  RES = RES + ABS( ( -WI*D1 )*X( 1, 1 )+ ( CA*A( 1, 1 )-WR*D1 )*X( 1, 2 )- SCALE*B( 1, 2 ) );
                                 if ( INFO == 0 ) {
                                    DEN = max( EPS*( max( ABS( CA*A( 1, 1 )-WR*D1 ), ABS( D1*WI ) )* ( ABS( X( 1, 1 ) )+ABS( X( 1, 2 ) ) ) ), SMLNUM );
                                 } else {
                                    DEN = max( SMIN*( ABS( X( 1, 1 ) )+ABS( X( 1, 2 ) ) ), SMLNUM );
                                 }
                                 RES = RES / DEN;
                                 if( ABS( X( 1, 1 ) ) < UNFL && ABS( X( 1, 2 ) ) < UNFL && ABS( B( 1, 1 ) ) <= SMLNUM* ABS( CA*A( 1, 1 )-WR*D1 ) ) RES = ZERO;
                                 if (SCALE > ONE) RES = RES + ONE / EPS;
                                 RES = RES + ABS( XNORM- ABS( X( 1, 1 ) )- ABS( X( 1, 2 ) ) ) / max( SMLNUM, XNORM ) / EPS;
                                 if (INFO != 0 && INFO != 1) RES = RES + ONE / EPS;
                                 KNT = KNT + 1;
                                 if ( RES > RMAX ) {
                                    LMAX = KNT;
                                    RMAX = RES;
                                 }
                              } // 40
                           } // 50
                        } // 60
                     } // 70

                     NA = 2;
                     NW = 1;
                     for (IA = 1; IA <= 3; IA++) { // 100
                        A( 1, 1 ) = VAB( IA );
                        A( 1, 2 ) = -THREE*VAB( IA );
                        A( 2, 1 ) = -SEVEN*VAB( IA );
                        A( 2, 2 ) = TWNONE*VAB( IA );
                        for (IB = 1; IB <= 3; IB++) { // 90
                           B( 1, 1 ) = VAB( IB );
                           B( 2, 1 ) = -TWO*VAB( IB );
                           for (IWR = 1; IWR <= 4; IWR++) { // 80
                              if ( D1 == ONE && D2 == ONE && CA == ONE ) {
                                 WR = VWR( IWR )*A( 1, 1 );
                              } else {
                                 WR = VWR( IWR );
                              }
                              WI = ZERO;
                              dlaln2(LTRANS( ITRANS ), NA, NW, SMIN, CA, A, 2, D1, D2, B, 2, WR, WI, X, 2, SCALE, XNORM, INFO );
                              if (INFO < 0) NINFO( 1 ) = NINFO( 1 ) + 1;
                              IF( INFO > 0 ) NINFO( 2 ) = NINFO( 2 ) + 1;
                              if ( ITRANS == 1 ) {
                                 TMP = A( 1, 2 );
                                 A( 1, 2 ) = A( 2, 1 );
                                 A( 2, 1 ) = TMP;
                              }
                              RES = ABS( ( CA*A( 1, 1 )-WR*D1 )* X( 1, 1 )+( CA*A( 1, 2 ) )* X( 2, 1 )-SCALE*B( 1, 1 ) )                               RES = RES + ABS( ( CA*A( 2, 1 ) )* X( 1, 1 )+( CA*A( 2, 2 )-WR*D2 )* X( 2, 1 )-SCALE*B( 2, 1 ) );
                              if ( INFO == 0 ) {
                                 DEN = max( EPS*( max( ABS( CA*A( 1, 1 )-WR*D1 )+ABS( CA*A( 1, 2 ) ), ABS( CA*A( 2, 1 ) )+ABS( CA*A( 2, 2 )-WR*D2 ) )*max( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) ) ), SMLNUM );
                              } else {
                                 DEN = max( EPS*( max( SMIN / EPS, max( ABS( CA*A( 1, 1 )-WR*D1 )+ABS( CA*A( 1, 2 ) ), ABS( CA*A( 2, 1 ) )+ABS( CA*A( 2, 2 )-WR*D2 ) ) )*max( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) ) ), SMLNUM );
                              }
                              RES = RES / DEN;
                              if( ABS( X( 1, 1 ) ) < UNFL && ABS( X( 2, 1 ) ) < UNFL && ABS( B( 1, 1 ) )+ABS( B( 2, 1 ) ) <= SMLNUM*( ABS( CA*A( 1, 1 )-WR*D1 )+ABS( CA*A( 1, 2 ) )+ABS( CA*A( 2, 1 ) )+ABS( CA*A( 2, 2 )-WR*D2 ) ) ) RES = ZERO;
                              if (SCALE > ONE) RES = RES + ONE / EPS;
                              RES = RES + ABS( XNORM- max( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) ) ) / max( SMLNUM, XNORM ) / EPS;
                              if (INFO != 0 && INFO != 1) RES = RES + ONE / EPS;
                              KNT = KNT + 1;
                              if ( RES > RMAX ) {
                                 LMAX = KNT;
                                 RMAX = RES;
                              }
                           } // 80
                        } // 90
                     } // 100

                     NA = 2;
                     NW = 2;
                     for (IA = 1; IA <= 3; IA++) { // 140
                        A( 1, 1 ) = VAB( IA )*TWO;
                        A( 1, 2 ) = -THREE*VAB( IA );
                        A( 2, 1 ) = -SEVEN*VAB( IA );
                        A( 2, 2 ) = TWNONE*VAB( IA );
                        for (IB = 1; IB <= 3; IB++) { // 130
                           B( 1, 1 ) = VAB( IB );
                           B( 2, 1 ) = -TWO*VAB( IB );
                           B( 1, 2 ) = FOUR*VAB( IB );
                           B( 2, 2 ) = -SEVEN*VAB( IB );
                           for (IWR = 1; IWR <= 4; IWR++) { // 120
                              if ( D1 == ONE && D2 == ONE && CA == ONE ) {
                                 WR = VWR( IWR )*A( 1, 1 );
                              } else {
                                 WR = VWR( IWR );
                              }
                              for (IWI = 1; IWI <= 4; IWI++) { // 110
                                 if ( D1 == ONE && D2 == ONE && CA == ONE ) {
                                    WI = VWI( IWI )*A( 1, 1 );
                                 } else {
                                    WI = VWI( IWI );
                                 }
                                 dlaln2(LTRANS( ITRANS ), NA, NW, SMIN, CA, A, 2, D1, D2, B, 2, WR, WI, X, 2, SCALE, XNORM, INFO );
                                 if (INFO < 0) NINFO( 1 ) = NINFO( 1 ) + 1;
                                 IF( INFO > 0 ) NINFO( 2 ) = NINFO( 2 ) + 1;
                                 if ( ITRANS == 1 ) {
                                    TMP = A( 1, 2 );
                                    A( 1, 2 ) = A( 2, 1 );
                                    A( 2, 1 ) = TMP;
                                 }
                                 RES = ABS( ( CA*A( 1, 1 )-WR*D1 )* X( 1, 1 )+( CA*A( 1, 2 ) )* X( 2, 1 )+( WI*D1 )*X( 1, 2 )- SCALE*B( 1, 1 ) )                                  RES = RES + ABS( ( CA*A( 1, 1 )-WR*D1 )*X( 1, 2 )+ ( CA*A( 1, 2 ) )*X( 2, 2 )- ( WI*D1 )*X( 1, 1 )-SCALE* B( 1, 2 ) );
                                 RES = RES + ABS( ( CA*A( 2, 1 ) )* X( 1, 1 )+( CA*A( 2, 2 )-WR*D2 )* X( 2, 1 )+( WI*D2 )*X( 2, 2 )- SCALE*B( 2, 1 ) )                                  RES = RES + ABS( ( CA*A( 2, 1 ) )* X( 1, 2 )+( CA*A( 2, 2 )-WR*D2 )* X( 2, 2 )-( WI*D2 )*X( 2, 1 )- SCALE*B( 2, 2 ) );
                                 if ( INFO == 0 ) {
                                    DEN = max( EPS*( max( ABS( CA*A( 1, 1 )-WR*D1 )+ABS( CA*A( 1, 2 ) )+ABS( WI*D1 ), ABS( CA*A( 2, 1 ) )+ABS( CA*A( 2, 2 )-WR*D2 )+ABS( WI*D2 ) )* max( ABS( X( 1, 1 ) )+ABS( X( 2, 1 ) ), ABS( X( 1, 2 ) )+ABS( X( 2, 2 ) ) ) ), SMLNUM );
                                 } else {
                                    DEN = max( EPS*( max( SMIN / EPS, max( ABS( CA*A( 1, 1 )-WR*D1 )+ABS( CA*A( 1, 2 ) )+ABS( WI*D1 ), ABS( CA*A( 2, 1 ) )+ABS( CA*A( 2, 2 )-WR*D2 )+ABS( WI*D2 ) ) )* max( ABS( X( 1, 1 ) )+ABS( X( 2, 1 ) ), ABS( X( 1, 2 ) )+ABS( X( 2, 2 ) ) ) ), SMLNUM );
                                 }
                                 RES = RES / DEN;
                                 if( ABS( X( 1, 1 ) ) < UNFL && ABS( X( 2, 1 ) ) < UNFL && ABS( X( 1, 2 ) ) < UNFL && ABS( X( 2, 2 ) ) < UNFL && ABS( B( 1, 1 ) )+ ABS( B( 2, 1 ) ) <= SMLNUM* ( ABS( CA*A( 1, 1 )-WR*D1 )+ ABS( CA*A( 1, 2 ) )+ABS( CA*A( 2, 1 ) )+ABS( CA*A( 2, 2 )-WR*D2 )+ABS( WI*D2 )+ABS( WI* D1 ) ) )RES = ZERO;
                                 if (SCALE > ONE) RES = RES + ONE / EPS;
                                 RES = RES + ABS( XNORM- max( ABS( X( 1, 1 ) )+ABS( X( 1, 2 ) ), ABS( X( 2, 1 ) )+ABS( X( 2, 2 ) ) ) ) / max( SMLNUM, XNORM ) / EPS;
                                 if (INFO != 0 && INFO != 1) RES = RES + ONE / EPS;
                                 KNT = KNT + 1;
                                 if ( RES > RMAX ) {
                                    LMAX = KNT;
                                    RMAX = RES;
                                 }
                              } // 110
                           } // 120
                        } // 130
                     } // 140
                  } // 150
               } // 160
            } // 170
         } // 180
      } // 190

      return;

      // End of DGET31

      }
