      SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR, LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LTRANL, LTRANR;
      int                INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2;
      double             SCALE, XNORM;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ), X( LDX, * );
      // ..

* =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      double             TWO, HALF, EIGHT;
      const              TWO = 2.0D+0, HALF = 0.5D+0, EIGHT = 8.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               BSWAP, XSWAP;
      int                I, IP, IPIV, IPSV, J, JP, JPSV, K;
      double             BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1, TEMP, U11, U12, U22, XMAX;
      // ..
      // .. Local Arrays ..
      bool               BSWPIV( 4 ), XSWPIV( 4 );
      int                JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ), LOCU22( 4 );
      double             BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 );
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL IDAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Data statements ..
      DATA               LOCU12 / 3, 4, 1, 2 / , LOCL21 / 2, 1, 4, 3 / , LOCU22 / 4, 3, 2, 1 /
      DATA               XSWPIV / false , false , true , true /
      DATA               BSWPIV / false , true , false , true /
      // ..
      // .. Executable Statements ..

      // Do not check the input parameters for errors

      INFO = 0

      // Quick return if possible

      if (N1 == 0 || N2 == 0) RETURN;

      // Set constants to control overflow

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SGN = ISGN

      K = N1 + N1 + N2 - 2
      GO TO ( 10, 20, 30, 50 )K

      // 1 by 1: TL11*X + SGN*X*TR11 = B11

      } // 10
      TAU1 = TL( 1, 1 ) + SGN*TR( 1, 1 )
      BET = ABS( TAU1 )
      if ( BET.LE.SMLNUM ) {
         TAU1 = SMLNUM
         BET = SMLNUM
         INFO = 1
      }

      SCALE = ONE
      GAM = ABS( B( 1, 1 ) )
      if (SMLNUM*GAM.GT.BET) SCALE = ONE / GAM;

      X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / TAU1
      XNORM = ABS( X( 1, 1 ) )
      RETURN

      // 1 by 2:
      // TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
                                        // [TR21 TR22]

      } // 20

      SMIN = MAX( EPS*MAX( ABS( TL( 1, 1 ) ), ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ), ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) ), SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      if ( LTRANR ) {
         TMP( 2 ) = SGN*TR( 2, 1 )
         TMP( 3 ) = SGN*TR( 1, 2 )
      } else {
         TMP( 2 ) = SGN*TR( 1, 2 )
         TMP( 3 ) = SGN*TR( 2, 1 )
      }
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 1, 2 )
      GO TO 40

      // 2 by 1:
           // op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
             // [TL21 TL22] [X21]         [X21]         [B21]

      } // 30
      SMIN = MAX( EPS*MAX( ABS( TR( 1, 1 ) ), ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ), ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) ), SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      if ( LTRANL ) {
         TMP( 2 ) = TL( 1, 2 )
         TMP( 3 ) = TL( 2, 1 )
      } else {
         TMP( 2 ) = TL( 2, 1 )
         TMP( 3 ) = TL( 1, 2 )
      }
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
      } // 40

      // Solve 2 by 2 system using complete pivoting.
      // Set pivots less than SMIN to SMIN.

      IPIV = IDAMAX( 4, TMP, 1 )
      U11 = TMP( IPIV )
      if ( ABS( U11 ).LE.SMIN ) {
         INFO = 1
         U11 = SMIN
      }
      U12 = TMP( LOCU12( IPIV ) )
      L21 = TMP( LOCL21( IPIV ) ) / U11
      U22 = TMP( LOCU22( IPIV ) ) - U12*L21
      XSWAP = XSWPIV( IPIV )
      BSWAP = BSWPIV( IPIV )
      if ( ABS( U22 ).LE.SMIN ) {
         INFO = 1
         U22 = SMIN
      }
      if ( BSWAP ) {
         TEMP = BTMP( 2 )
         BTMP( 2 ) = BTMP( 1 ) - L21*TEMP
         BTMP( 1 ) = TEMP
      } else {
         BTMP( 2 ) = BTMP( 2 ) - L21*BTMP( 1 )
      }
      SCALE = ONE
      if ( ( TWO*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( U22 ) || ( TWO*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( U11 ) ) {
         SCALE = HALF / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
      }
      X2( 2 ) = BTMP( 2 ) / U22
      X2( 1 ) = BTMP( 1 ) / U11 - ( U12 / U11 )*X2( 2 )
      if ( XSWAP ) {
         TEMP = X2( 2 )
         X2( 2 ) = X2( 1 )
         X2( 1 ) = TEMP
      }
      X( 1, 1 ) = X2( 1 )
      if ( N1 == 1 ) {
         X( 1, 2 ) = X2( 2 )
         XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
      } else {
         X( 2, 1 ) = X2( 2 )
         XNORM = MAX( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) )
      }
      RETURN

      // 2 by 2:
      // op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
        // [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]

      // Solve equivalent 4 by 4 system using complete pivoting.
      // Set pivots less than SMIN to SMIN.

      } // 50
      SMIN = MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ), ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )       SMIN = MAX( SMIN, ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ), ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )
      SMIN = MAX( EPS*SMIN, SMLNUM )
      BTMP( 1 ) = ZERO
      dcopy(16, BTMP, 0, T16, 1 );
      T16( 1, 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      T16( 2, 2 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      T16( 3, 3 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      T16( 4, 4 ) = TL( 2, 2 ) + SGN*TR( 2, 2 )
      if ( LTRANL ) {
         T16( 1, 2 ) = TL( 2, 1 )
         T16( 2, 1 ) = TL( 1, 2 )
         T16( 3, 4 ) = TL( 2, 1 )
         T16( 4, 3 ) = TL( 1, 2 )
      } else {
         T16( 1, 2 ) = TL( 1, 2 )
         T16( 2, 1 ) = TL( 2, 1 )
         T16( 3, 4 ) = TL( 1, 2 )
         T16( 4, 3 ) = TL( 2, 1 )
      }
      if ( LTRANR ) {
         T16( 1, 3 ) = SGN*TR( 1, 2 )
         T16( 2, 4 ) = SGN*TR( 1, 2 )
         T16( 3, 1 ) = SGN*TR( 2, 1 )
         T16( 4, 2 ) = SGN*TR( 2, 1 )
      } else {
         T16( 1, 3 ) = SGN*TR( 2, 1 )
         T16( 2, 4 ) = SGN*TR( 2, 1 )
         T16( 3, 1 ) = SGN*TR( 1, 2 )
         T16( 4, 2 ) = SGN*TR( 1, 2 )
      }
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
      BTMP( 3 ) = B( 1, 2 )
      BTMP( 4 ) = B( 2, 2 )

      // Perform elimination

      for (I = 1; I <= 3; I++) { // 100
         XMAX = ZERO
         for (IP = I; IP <= 4; IP++) { // 70
            for (JP = I; JP <= 4; JP++) { // 60
               if ( ABS( T16( IP, JP ) ).GE.XMAX ) {
                  XMAX = ABS( T16( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               }
            } // 60
         } // 70
         if ( IPSV != I ) {
            dswap(4, T16( IPSV, 1 ), 4, T16( I, 1 ), 4 );
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         }
         if (JPSV != I) CALL DSWAP( 4, T16( 1, JPSV ), 1, T16( 1, I ), 1 );
         JPIV( I ) = JPSV
         if ( ABS( T16( I, I ) ) < SMIN ) {
            INFO = 1
            T16( I, I ) = SMIN
         }
         for (J = I + 1; J <= 4; J++) { // 90
            T16( J, I ) = T16( J, I ) / T16( I, I )
            BTMP( J ) = BTMP( J ) - T16( J, I )*BTMP( I )
            for (K = I + 1; K <= 4; K++) { // 80
               T16( J, K ) = T16( J, K ) - T16( J, I )*T16( I, K )
            } // 80
         } // 90
      } // 100
      if ( ABS( T16( 4, 4 ) ) < SMIN ) {
         INFO = 1
         T16( 4, 4 ) = SMIN
      }
      SCALE = ONE
      if ( ( EIGHT*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T16( 1, 1 ) ) || ( EIGHT*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T16( 2, 2 ) ) || ( EIGHT*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T16( 3, 3 ) ) || ( EIGHT*SMLNUM )*ABS( BTMP( 4 ) ).GT.ABS( T16( 4, 4 ) ) ) {
         SCALE = ( ONE / EIGHT ) / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ), ABS( BTMP( 4 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
         BTMP( 4 ) = BTMP( 4 )*SCALE
      }
      for (I = 1; I <= 4; I++) { // 120
         K = 5 - I
         TEMP = ONE / T16( K, K )
         TMP( K ) = BTMP( K )*TEMP
         for (J = K + 1; J <= 4; J++) { // 110
            TMP( K ) = TMP( K ) - ( TEMP*T16( K, J ) )*TMP( J )
         } // 110
      } // 120
      for (I = 1; I <= 3; I++) { // 130
         if ( JPIV( 4-I ) != 4-I ) {
            TEMP = TMP( 4-I )
            TMP( 4-I ) = TMP( JPIV( 4-I ) )
            TMP( JPIV( 4-I ) ) = TEMP
         }
      } // 130
      X( 1, 1 ) = TMP( 1 )
      X( 2, 1 ) = TMP( 2 )
      X( 1, 2 ) = TMP( 3 )
      X( 2, 2 ) = TMP( 4 )
      XNORM = MAX( ABS( TMP( 1 ) )+ABS( TMP( 3 ) ), ABS( TMP( 2 ) )+ABS( TMP( 4 ) ) )
      RETURN

      // End of DLASY2

      }
