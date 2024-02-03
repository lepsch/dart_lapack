      SUBROUTINE DGET32( RMAX, LMAX, NINFO, KNT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX, NINFO;
      double             RMAX;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      double             TWO, FOUR, EIGHT;
      const              TWO = 2.0D0, FOUR = 4.0D0, EIGHT = 8.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               LTRANL, LTRANR;
      int                IB, IB1, IB2, IB3, INFO, ISGN, ITL, ITLSCL, ITR, ITRANL, ITRANR, ITRSCL, N1, N2;
      double             BIGNUM, DEN, EPS, RES, SCALE, SGN, SMLNUM, TMP, TNRM, XNORM, XNRM;
      // ..
      // .. Local Arrays ..
      int                ITVAL( 2, 2, 8 );
      double             B( 2, 2 ), TL( 2, 2 ), TR( 2, 2 ), VAL( 3 ), X( 2, 2 );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASY2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               ITVAL / 8, 4, 2, 1, 4, 8, 1, 2, 2, 1, 8, 4, 1, 2, 4, 8, 9, 4, 2, 1, 4, 9, 1, 2, 2, 1, 9, 4, 1, 2, 4, 9 /
      // ..
      // .. Executable Statements ..

      // Get machine parameters

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      // Set up test case parameters

      VAL( 1 ) = SQRT( SMLNUM )
      VAL( 2 ) = ONE
      VAL( 3 ) = SQRT( BIGNUM )

      KNT = 0
      NINFO = 0
      LMAX = 0
      RMAX = ZERO

      // Begin test loop

      for (ITRANL = 0; ITRANL <= 1; ITRANL++) { // 230
         for (ITRANR = 0; ITRANR <= 1; ITRANR++) { // 220
            DO 210 ISGN = -1, 1, 2
               SGN = ISGN
               LTRANL = ITRANL == 1
               LTRANR = ITRANR == 1

               N1 = 1
               N2 = 1
               for (ITL = 1; ITL <= 3; ITL++) { // 30
                  for (ITR = 1; ITR <= 3; ITR++) { // 20
                     for (IB = 1; IB <= 3; IB++) { // 10
                        TL( 1, 1 ) = VAL( ITL )
                        TR( 1, 1 ) = VAL( ITR )
                        B( 1, 1 ) = VAL( IB )
                        KNT = KNT + 1
                        dlasy2(LTRANL, LTRANR, ISGN, N1, N2, TL, 2, TR, 2, B, 2, SCALE, X, 2, XNORM, INFO );
                        if (INFO != 0) NINFO = NINFO + 1                         RES = ABS( ( TL( 1, 1 )+SGN*TR( 1, 1 ) )* X( 1, 1 )-SCALE*B( 1, 1 ) );
                        if ( INFO == 0 ) {
                           DEN = MAX( EPS*( ( ABS( TR( 1, 1 ) )+ABS( TL( 1, 1 ) ) )*ABS( X( 1, 1 ) ) ), SMLNUM )
                        } else {
                           DEN = SMLNUM*MAX( ABS( X( 1, 1 ) ), ONE )
                        }
                        RES = RES / DEN
                        if (SCALE > ONE) RES = RES + ONE / EPS                         RES = RES + ABS( XNORM-ABS( X( 1, 1 ) ) ) / MAX( SMLNUM, XNORM ) / EPS                         IF( INFO != 0 && INFO != 1 ) RES = RES + ONE / EPS;
                        if ( RES > RMAX ) {
                           LMAX = KNT
                           RMAX = RES
                        }
                     } // 10
                  } // 20
               } // 30

               N1 = 2
               N2 = 1
               for (ITL = 1; ITL <= 8; ITL++) { // 80
                  for (ITLSCL = 1; ITLSCL <= 3; ITLSCL++) { // 70
                     for (ITR = 1; ITR <= 3; ITR++) { // 60
                        for (IB1 = 1; IB1 <= 3; IB1++) { // 50
                           for (IB2 = 1; IB2 <= 3; IB2++) { // 40
                              B( 1, 1 ) = VAL( IB1 )
                              B( 2, 1 ) = -FOUR*VAL( IB2 )
                              TL( 1, 1 ) = ITVAL( 1, 1, ITL )* VAL( ITLSCL )                               TL( 2, 1 ) = ITVAL( 2, 1, ITL )* VAL( ITLSCL )                               TL( 1, 2 ) = ITVAL( 1, 2, ITL )* VAL( ITLSCL )                               TL( 2, 2 ) = ITVAL( 2, 2, ITL )* VAL( ITLSCL )
                              TR( 1, 1 ) = VAL( ITR )
                              KNT = KNT + 1
                              dlasy2(LTRANL, LTRANR, ISGN, N1, N2, TL, 2, TR, 2, B, 2, SCALE, X, 2, XNORM, INFO );
                              if (INFO != 0) NINFO = NINFO + 1;
                              if ( LTRANL ) {
                                 TMP = TL( 1, 2 )
                                 TL( 1, 2 ) = TL( 2, 1 )
                                 TL( 2, 1 ) = TMP
                              }
                              RES = ABS( ( TL( 1, 1 )+SGN*TR( 1, 1 ) )* X( 1, 1 )+TL( 1, 2 )*X( 2, 1 )- SCALE*B( 1, 1 ) )                               RES = RES + ABS( ( TL( 2, 2 )+SGN*TR( 1, 1 ) )*X( 2, 1 )+TL( 2, 1 )* X( 1, 1 )-SCALE*B( 2, 1 ) )                               TNRM = ABS( TR( 1, 1 ) ) + ABS( TL( 1, 1 ) ) + ABS( TL( 1, 2 ) ) + ABS( TL( 2, 1 ) ) + ABS( TL( 2, 2 ) )
                              XNRM = MAX( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) )                               DEN = MAX( SMLNUM, SMLNUM*XNRM, ( TNRM*EPS )*XNRM )
                              RES = RES / DEN
                              if (SCALE > ONE) RES = RES + ONE / EPS                               RES = RES + ABS( XNORM-XNRM ) / MAX( SMLNUM, XNORM ) / EPS;
                              if ( RES > RMAX ) {
                                 LMAX = KNT
                                 RMAX = RES
                              }
                           } // 40
                        } // 50
                     } // 60
                  } // 70
               } // 80

               N1 = 1
               N2 = 2
               for (ITR = 1; ITR <= 8; ITR++) { // 130
                  for (ITRSCL = 1; ITRSCL <= 3; ITRSCL++) { // 120
                     for (ITL = 1; ITL <= 3; ITL++) { // 110
                        for (IB1 = 1; IB1 <= 3; IB1++) { // 100
                           for (IB2 = 1; IB2 <= 3; IB2++) { // 90
                              B( 1, 1 ) = VAL( IB1 )
                              B( 1, 2 ) = -TWO*VAL( IB2 )
                              TR( 1, 1 ) = ITVAL( 1, 1, ITR )* VAL( ITRSCL )                               TR( 2, 1 ) = ITVAL( 2, 1, ITR )* VAL( ITRSCL )                               TR( 1, 2 ) = ITVAL( 1, 2, ITR )* VAL( ITRSCL )                               TR( 2, 2 ) = ITVAL( 2, 2, ITR )* VAL( ITRSCL )
                              TL( 1, 1 ) = VAL( ITL )
                              KNT = KNT + 1
                              dlasy2(LTRANL, LTRANR, ISGN, N1, N2, TL, 2, TR, 2, B, 2, SCALE, X, 2, XNORM, INFO );
                              if (INFO != 0) NINFO = NINFO + 1;
                              if ( LTRANR ) {
                                 TMP = TR( 1, 2 )
                                 TR( 1, 2 ) = TR( 2, 1 )
                                 TR( 2, 1 ) = TMP
                              }
                              TNRM = ABS( TL( 1, 1 ) ) + ABS( TR( 1, 1 ) ) + ABS( TR( 1, 2 ) ) + ABS( TR( 2, 2 ) ) + ABS( TR( 2, 1 ) )
                              XNRM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
                              RES = ABS( ( ( TL( 1, 1 )+SGN*TR( 1, 1 ) ) )*( X( 1, 1 ) )+ ( SGN*TR( 2, 1 ) )*( X( 1, 2 ) )- ( SCALE*B( 1, 1 ) ) )                               RES = RES + ABS( ( ( TL( 1, 1 )+SGN*TR( 2, 2 ) ) )*( X( 1, 2 ) )+ ( SGN*TR( 1, 2 ) )*( X( 1, 1 ) )- ( SCALE*B( 1, 2 ) ) )
                              DEN = MAX( SMLNUM, SMLNUM*XNRM, ( TNRM*EPS )*XNRM )
                              RES = RES / DEN
                              if (SCALE > ONE) RES = RES + ONE / EPS                               RES = RES + ABS( XNORM-XNRM ) / MAX( SMLNUM, XNORM ) / EPS;
                              if ( RES > RMAX ) {
                                 LMAX = KNT
                                 RMAX = RES
                              }
                           } // 90
                        } // 100
                     } // 110
                  } // 120
               } // 130

               N1 = 2
               N2 = 2
               for (ITR = 1; ITR <= 8; ITR++) { // 200
                  for (ITRSCL = 1; ITRSCL <= 3; ITRSCL++) { // 190
                     for (ITL = 1; ITL <= 8; ITL++) { // 180
                        for (ITLSCL = 1; ITLSCL <= 3; ITLSCL++) { // 170
                           for (IB1 = 1; IB1 <= 3; IB1++) { // 160
                              for (IB2 = 1; IB2 <= 3; IB2++) { // 150
                                 for (IB3 = 1; IB3 <= 3; IB3++) { // 140
                                    B( 1, 1 ) = VAL( IB1 )
                                    B( 2, 1 ) = -FOUR*VAL( IB2 )
                                    B( 1, 2 ) = -TWO*VAL( IB3 )
                                    B( 2, 2 ) = EIGHT* MIN( VAL( IB1 ), VAL ( IB2 ), VAL( IB3 ) )
                                    TR( 1, 1 ) = ITVAL( 1, 1, ITR )* VAL( ITRSCL )                                     TR( 2, 1 ) = ITVAL( 2, 1, ITR )* VAL( ITRSCL )                                     TR( 1, 2 ) = ITVAL( 1, 2, ITR )* VAL( ITRSCL )                                     TR( 2, 2 ) = ITVAL( 2, 2, ITR )* VAL( ITRSCL )                                     TL( 1, 1 ) = ITVAL( 1, 1, ITL )* VAL( ITLSCL )                                     TL( 2, 1 ) = ITVAL( 2, 1, ITL )* VAL( ITLSCL )                                     TL( 1, 2 ) = ITVAL( 1, 2, ITL )* VAL( ITLSCL )                                     TL( 2, 2 ) = ITVAL( 2, 2, ITL )* VAL( ITLSCL )
                                    KNT = KNT + 1
                                    dlasy2(LTRANL, LTRANR, ISGN, N1, N2, TL, 2, TR, 2, B, 2, SCALE, X, 2, XNORM, INFO );
                                    if (INFO != 0) NINFO = NINFO + 1;
                                    if ( LTRANR ) {
                                       TMP = TR( 1, 2 )
                                       TR( 1, 2 ) = TR( 2, 1 )
                                       TR( 2, 1 ) = TMP
                                    }
                                    if ( LTRANL ) {
                                       TMP = TL( 1, 2 )
                                       TL( 1, 2 ) = TL( 2, 1 )
                                       TL( 2, 1 ) = TMP
                                    }
                                    TNRM = ABS( TR( 1, 1 ) ) + ABS( TR( 2, 1 ) ) + ABS( TR( 1, 2 ) ) + ABS( TR( 2, 2 ) ) + ABS( TL( 1, 1 ) ) + ABS( TL( 2, 1 ) ) + ABS( TL( 1, 2 ) ) + ABS( TL( 2, 2 ) )                                     XNRM = MAX( ABS( X( 1, 1 ) )+ ABS( X( 1, 2 ) ), ABS( X( 2, 1 ) )+ ABS( X( 2, 2 ) ) )                                     RES = ABS( ( ( TL( 1, 1 )+SGN*TR( 1, 1 ) ) )*( X( 1, 1 ) )+ ( SGN*TR( 2, 1 ) )* ( X( 1, 2 ) )+( TL( 1, 2 ) )* ( X( 2, 1 ) )- ( SCALE*B( 1, 1 ) ) )                                     RES = RES + ABS( ( TL( 1, 1 ) )* ( X( 1, 2 ) )+ ( SGN*TR( 1, 2 ) )* ( X( 1, 1 ) )+ ( SGN*TR( 2, 2 ) )* ( X( 1, 2 ) )+( TL( 1, 2 ) )* ( X( 2, 2 ) )- ( SCALE*B( 1, 2 ) ) )                                     RES = RES + ABS( ( TL( 2, 1 ) )* ( X( 1, 1 ) )+ ( SGN*TR( 1, 1 ) )* ( X( 2, 1 ) )+ ( SGN*TR( 2, 1 ) )* ( X( 2, 2 ) )+( TL( 2, 2 ) )* ( X( 2, 1 ) )- ( SCALE*B( 2, 1 ) ) )
                                    RES = RES + ABS( ( ( TL( 2, 2 )+SGN*TR( 2, 2 ) ) )* ( X( 2, 2 ) )+ ( SGN*TR( 1, 2 ) )* ( X( 2, 1 ) )+( TL( 2, 1 ) )* ( X( 1, 2 ) )- ( SCALE*B( 2, 2 ) ) )
                                    DEN = MAX( SMLNUM, SMLNUM*XNRM, ( TNRM*EPS )*XNRM )
                                    RES = RES / DEN
                                    if (SCALE > ONE) RES = RES + ONE / EPS                                     RES = RES + ABS( XNORM-XNRM ) / MAX( SMLNUM, XNORM ) / EPS;
                                    if ( RES > RMAX ) {
                                       LMAX = KNT
                                       RMAX = RES
                                    }
                                 } // 140
                              } // 150
                           } // 160
                        } // 170
                     } // 180
                  } // 190
               } // 200
            } // 210
         } // 220
      } // 230

      RETURN

      // End of DGET32

      }
