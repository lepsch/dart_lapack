      SUBROUTINE DGET33( RMAX, LMAX, NINFO, KNT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX, NINFO;
      double             RMAX;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             TWO, FOUR;
      const              TWO = 2.0, FOUR = 4.0 ;
      // ..
      // .. Local Scalars ..
      int                I1, I2, I3, I4, IM1, IM2, IM3, IM4, J1, J2, J3;
      double             BIGNUM, CS, EPS, RES, SMLNUM, SN, SUM, TNRM, WI1, WI2, WR1, WR2;
      // ..
      // .. Local Arrays ..
      double             Q( 2, 2 ), T( 2, 2 ), T1( 2, 2 ), T2( 2, 2 ), VAL( 4 ), VM( 3 );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLANV2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SIGN
      // ..
      // .. Executable Statements ..

      // Get machine parameters

      EPS = DLAMCH( 'P' );
      SMLNUM = DLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Set up test case parameters

      VAL( 1 ) = ONE;
      VAL( 2 ) = ONE + TWO*EPS;
      VAL( 3 ) = TWO;
      VAL( 4 ) = TWO - FOUR*EPS;
      VM( 1 ) = SMLNUM;
      VM( 2 ) = ONE;
      VM( 3 ) = BIGNUM;

      KNT = 0;
      NINFO = 0;
      LMAX = 0;
      RMAX = ZERO;

      // Begin test loop

      for (I1 = 1; I1 <= 4; I1++) { // 150
         for (I2 = 1; I2 <= 4; I2++) { // 140
            for (I3 = 1; I3 <= 4; I3++) { // 130
               for (I4 = 1; I4 <= 4; I4++) { // 120
                  for (IM1 = 1; IM1 <= 3; IM1++) { // 110
                     for (IM2 = 1; IM2 <= 3; IM2++) { // 100
                        for (IM3 = 1; IM3 <= 3; IM3++) { // 90
                           for (IM4 = 1; IM4 <= 3; IM4++) { // 80
                              T( 1, 1 ) = VAL( I1 )*VM( IM1 );
                              T( 1, 2 ) = VAL( I2 )*VM( IM2 );
                              T( 2, 1 ) = -VAL( I3 )*VM( IM3 );
                              T( 2, 2 ) = VAL( I4 )*VM( IM4 );
                              TNRM = max( ABS( T( 1, 1 ) ), ABS( T( 1, 2 ) ), ABS( T( 2, 1 ) ), ABS( T( 2, 2 ) ) );
                              T1( 1, 1 ) = T( 1, 1 );
                              T1( 1, 2 ) = T( 1, 2 );
                              T1( 2, 1 ) = T( 2, 1 );
                              T1( 2, 2 ) = T( 2, 2 );
                              Q( 1, 1 ) = ONE;
                              Q( 1, 2 ) = ZERO;
                              Q( 2, 1 ) = ZERO;
                              Q( 2, 2 ) = ONE;

                              dlanv2(T( 1, 1 ), T( 1, 2 ), T( 2, 1 ), T( 2, 2 ), WR1, WI1, WR2, WI2, CS, SN );
                              for (J1 = 1; J1 <= 2; J1++) { // 10
                                 RES = Q( J1, 1 )*CS + Q( J1, 2 )*SN;
                                 Q( J1, 2 ) = -Q( J1, 1 )*SN + Q( J1, 2 )*CS;
                                 Q( J1, 1 ) = RES;
                              } // 10

                              RES = ZERO;
                              RES = RES + ABS( Q( 1, 1 )**2+ Q( 1, 2 )**2-ONE ) / EPS                               RES = RES + ABS( Q( 2, 2 )**2+ Q( 2, 1 )**2-ONE ) / EPS                               RES = RES + ABS( Q( 1, 1 )*Q( 2, 1 )+ Q( 1, 2 )*Q( 2, 2 ) ) / EPS;
                              for (J1 = 1; J1 <= 2; J1++) { // 40
                                 for (J2 = 1; J2 <= 2; J2++) { // 30
                                    T2( J1, J2 ) = ZERO;
                                    for (J3 = 1; J3 <= 2; J3++) { // 20
                                       T2( J1, J2 ) = T2( J1, J2 ) + T1( J1, J3 )* Q( J3, J2 );
                                    } // 20
                                 } // 30
                              } // 40
                              for (J1 = 1; J1 <= 2; J1++) { // 70
                                 for (J2 = 1; J2 <= 2; J2++) { // 60
                                    SUM = T( J1, J2 );
                                    for (J3 = 1; J3 <= 2; J3++) { // 50
                                       SUM = SUM - Q( J3, J1 )* T2( J3, J2 );
                                    } // 50
                                    RES = RES + ABS( SUM ) / EPS / TNRM;
                                 } // 60
                              } // 70
                              if( T( 2, 1 ) != ZERO && ( T( 1, 1 ) != T( 2, 2 ) || SIGN( ONE, T( 1, 2 ) )*SIGN( ONE, T( 2, 1 ) ) > ZERO ) )RES = RES + ONE / EPS;
                              KNT = KNT + 1;
                              if ( RES > RMAX ) {
                                 LMAX = KNT;
                                 RMAX = RES;
                              }
                           } // 80
                        } // 90
                     } // 100
                  } // 110
               } // 120
            } // 130
         } // 140
      } // 150

      return;

      // End of DGET33

      }
