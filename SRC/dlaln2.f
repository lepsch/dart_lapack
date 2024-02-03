      SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, LDB, WR, WI, X, LDX, SCALE, XNORM, INFO );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LTRANS;
      int                INFO, LDA, LDB, LDX, NA, NW;
      double             CA, D1, D2, SCALE, SMIN, WI, WR, XNORM;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             TWO;
      const              TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      int                ICMAX, J;
      double             BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21, CI22, CMAX, CNORM, CR21, CR22, CSI, CSR, LI21, LR21, SMINI, SMLNUM, TEMP, U22ABS, UI11, UI11R, UI12, UI12S, UI22, UR11, UR11R, UR12, UR12S, UR22, XI1, XI2, XR1, XR2;
      // ..
      // .. Local Arrays ..
      bool               RSWAP( 4 ), ZSWAP( 4 );
      int                IPIVOT( 4, 4 );
      double             CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Equivalences ..
      EQUIVALENCE        ( CI( 1, 1 ), CIV( 1 ) ), ( CR( 1, 1 ), CRV( 1 ) );
      // ..
      // .. Data statements ..
      DATA               ZSWAP / false , false , true , true /;
      DATA               RSWAP / false , true , false , true /;
      DATA               IPIVOT / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1 /;
      // ..
      // .. Executable Statements ..

      // Compute BIGNUM

      SMLNUM = TWO*DLAMCH( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;
      SMINI = MAX( SMIN, SMLNUM );

      // Don't check for input errors

      INFO = 0;

      // Standard Initializations

      SCALE = ONE;

      if ( NA == 1 ) {

         // 1 x 1  (i.e., scalar) system   C X = B

         if ( NW == 1 ) {

            // Real 1x1 system.

            // C = ca A - w D

            CSR = CA*A( 1, 1 ) - WR*D1;
            CNORM = ABS( CSR );

            // If | C | < SMINI, use C = SMINI

            if ( CNORM < SMINI ) {
               CSR = SMINI;
               CNORM = SMINI;
               INFO = 1;
            }

            // Check scaling for  X = B / C

            BNORM = ABS( B( 1, 1 ) );
            if ( CNORM < ONE && BNORM > ONE ) {
               if (BNORM > BIGNUM*CNORM) SCALE = ONE / BNORM;
            }

            // Compute X

            X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / CSR;
            XNORM = ABS( X( 1, 1 ) );
         } else {

            // Complex 1x1 system (w is complex)

            // C = ca A - w D

            CSR = CA*A( 1, 1 ) - WR*D1;
            CSI = -WI*D1;
            CNORM = ABS( CSR ) + ABS( CSI );

            // If | C | < SMINI, use C = SMINI

            if ( CNORM < SMINI ) {
               CSR = SMINI;
               CSI = ZERO;
               CNORM = SMINI;
               INFO = 1;
            }

            // Check scaling for  X = B / C

            BNORM = ABS( B( 1, 1 ) ) + ABS( B( 1, 2 ) );
            if ( CNORM < ONE && BNORM > ONE ) {
               if (BNORM > BIGNUM*CNORM) SCALE = ONE / BNORM;
            }

            // Compute X

            dladiv(SCALE*B( 1, 1 ), SCALE*B( 1, 2 ), CSR, CSI, X( 1, 1 ), X( 1, 2 ) );
            XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) );
         }

      } else {

         // 2x2 System

         // Compute the real part of  C = ca A - w D  (or  ca A**T - w D )

         CR( 1, 1 ) = CA*A( 1, 1 ) - WR*D1;
         CR( 2, 2 ) = CA*A( 2, 2 ) - WR*D2;
         if ( LTRANS ) {
            CR( 1, 2 ) = CA*A( 2, 1 );
            CR( 2, 1 ) = CA*A( 1, 2 );
         } else {
            CR( 2, 1 ) = CA*A( 2, 1 );
            CR( 1, 2 ) = CA*A( 1, 2 );
         }

         if ( NW == 1 ) {

            // Real 2x2 system  (w is real)

            // Find the largest element in C

            CMAX = ZERO;
            ICMAX = 0;

            for (J = 1; J <= 4; J++) { // 10
               if ( ABS( CRV( J ) ) > CMAX ) {
                  CMAX = ABS( CRV( J ) );
                  ICMAX = J;
               }
            } // 10

            // If norm(C) < SMINI, use SMINI*identity.

            if ( CMAX < SMINI ) {
               BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 2, 1 ) ) );
               if ( SMINI < ONE && BNORM > ONE ) {
                  if (BNORM > BIGNUM*SMINI) SCALE = ONE / BNORM;
               }
               TEMP = SCALE / SMINI;
               X( 1, 1 ) = TEMP*B( 1, 1 );
               X( 2, 1 ) = TEMP*B( 2, 1 );
               XNORM = TEMP*BNORM;
               INFO = 1;
               return;
            }

            // Gaussian elimination with complete pivoting.

            UR11 = CRV( ICMAX );
            CR21 = CRV( IPIVOT( 2, ICMAX ) );
            UR12 = CRV( IPIVOT( 3, ICMAX ) );
            CR22 = CRV( IPIVOT( 4, ICMAX ) );
            UR11R = ONE / UR11;
            LR21 = UR11R*CR21;
            UR22 = CR22 - UR12*LR21;

            // If smaller pivot < SMINI, use SMINI

            if ( ABS( UR22 ) < SMINI ) {
               UR22 = SMINI;
               INFO = 1;
            }
            if ( RSWAP( ICMAX ) ) {
               BR1 = B( 2, 1 );
               BR2 = B( 1, 1 );
            } else {
               BR1 = B( 1, 1 );
               BR2 = B( 2, 1 );
            }
            BR2 = BR2 - LR21*BR1;
            BBND = MAX( ABS( BR1*( UR22*UR11R ) ), ABS( BR2 ) );
            if ( BBND > ONE && ABS( UR22 ) < ONE ) {
               IF( BBND >= BIGNUM*ABS( UR22 ) ) SCALE = ONE / BBND;
            }

            XR2 = ( BR2*SCALE ) / UR22;
            XR1 = ( SCALE*BR1 )*UR11R - XR2*( UR11R*UR12 );
            if ( ZSWAP( ICMAX ) ) {
               X( 1, 1 ) = XR2;
               X( 2, 1 ) = XR1;
            } else {
               X( 1, 1 ) = XR1;
               X( 2, 1 ) = XR2;
            }
            XNORM = MAX( ABS( XR1 ), ABS( XR2 ) );

            // Further scaling if  norm(A) norm(X) > overflow

            if ( XNORM > ONE && CMAX > ONE ) {
               if ( XNORM > BIGNUM / CMAX ) {
                  TEMP = CMAX / BIGNUM;
                  X( 1, 1 ) = TEMP*X( 1, 1 );
                  X( 2, 1 ) = TEMP*X( 2, 1 );
                  XNORM = TEMP*XNORM;
                  SCALE = TEMP*SCALE;
               }
            }
         } else {

            // Complex 2x2 system  (w is complex)

            // Find the largest element in C

            CI( 1, 1 ) = -WI*D1;
            CI( 2, 1 ) = ZERO;
            CI( 1, 2 ) = ZERO;
            CI( 2, 2 ) = -WI*D2;
            CMAX = ZERO;
            ICMAX = 0;

            for (J = 1; J <= 4; J++) { // 20
               if ( ABS( CRV( J ) )+ABS( CIV( J ) ) > CMAX ) {
                  CMAX = ABS( CRV( J ) ) + ABS( CIV( J ) );
                  ICMAX = J;
               }
            } // 20

            // If norm(C) < SMINI, use SMINI*identity.

            if ( CMAX < SMINI ) {
               BNORM = MAX( ABS( B( 1, 1 ) )+ABS( B( 1, 2 ) ), ABS( B( 2, 1 ) )+ABS( B( 2, 2 ) ) );
               if ( SMINI < ONE && BNORM > ONE ) {
                  if (BNORM > BIGNUM*SMINI) SCALE = ONE / BNORM;
               }
               TEMP = SCALE / SMINI;
               X( 1, 1 ) = TEMP*B( 1, 1 );
               X( 2, 1 ) = TEMP*B( 2, 1 );
               X( 1, 2 ) = TEMP*B( 1, 2 );
               X( 2, 2 ) = TEMP*B( 2, 2 );
               XNORM = TEMP*BNORM;
               INFO = 1;
               return;
            }

            // Gaussian elimination with complete pivoting.

            UR11 = CRV( ICMAX );
            UI11 = CIV( ICMAX );
            CR21 = CRV( IPIVOT( 2, ICMAX ) );
            CI21 = CIV( IPIVOT( 2, ICMAX ) );
            UR12 = CRV( IPIVOT( 3, ICMAX ) );
            UI12 = CIV( IPIVOT( 3, ICMAX ) );
            CR22 = CRV( IPIVOT( 4, ICMAX ) );
            CI22 = CIV( IPIVOT( 4, ICMAX ) );
            if ( ICMAX == 1 || ICMAX == 4 ) {

               // Code when off-diagonals of pivoted C are real

               if ( ABS( UR11 ) > ABS( UI11 ) ) {
                  TEMP = UI11 / UR11;
                  UR11R = ONE / ( UR11*( ONE+TEMP**2 ) );
                  UI11R = -TEMP*UR11R;
               } else {
                  TEMP = UR11 / UI11;
                  UI11R = -ONE / ( UI11*( ONE+TEMP**2 ) );
                  UR11R = -TEMP*UI11R;
               }
               LR21 = CR21*UR11R;
               LI21 = CR21*UI11R;
               UR12S = UR12*UR11R;
               UI12S = UR12*UI11R;
               UR22 = CR22 - UR12*LR21;
               UI22 = CI22 - UR12*LI21;
            } else {

               // Code when diagonals of pivoted C are real

               UR11R = ONE / UR11;
               UI11R = ZERO;
               LR21 = CR21*UR11R;
               LI21 = CI21*UR11R;
               UR12S = UR12*UR11R;
               UI12S = UI12*UR11R;
               UR22 = CR22 - UR12*LR21 + UI12*LI21;
               UI22 = -UR12*LI21 - UI12*LR21;
            }
            U22ABS = ABS( UR22 ) + ABS( UI22 );

            // If smaller pivot < SMINI, use SMINI

            if ( U22ABS < SMINI ) {
               UR22 = SMINI;
               UI22 = ZERO;
               INFO = 1;
            }
            if ( RSWAP( ICMAX ) ) {
               BR2 = B( 1, 1 );
               BR1 = B( 2, 1 );
               BI2 = B( 1, 2 );
               BI1 = B( 2, 2 );
            } else {
               BR1 = B( 1, 1 );
               BR2 = B( 2, 1 );
               BI1 = B( 1, 2 );
               BI2 = B( 2, 2 );
            }
            BR2 = BR2 - LR21*BR1 + LI21*BI1;
            BI2 = BI2 - LI21*BR1 - LR21*BI1;
            BBND = MAX( ( ABS( BR1 )+ABS( BI1 ) )* ( U22ABS*( ABS( UR11R )+ABS( UI11R ) ) ), ABS( BR2 )+ABS( BI2 ) );
            if ( BBND > ONE && U22ABS < ONE ) {
               if ( BBND >= BIGNUM*U22ABS ) {
                  SCALE = ONE / BBND;
                  BR1 = SCALE*BR1;
                  BI1 = SCALE*BI1;
                  BR2 = SCALE*BR2;
                  BI2 = SCALE*BI2;
               }
            }

            dladiv(BR2, BI2, UR22, UI22, XR2, XI2 );
            XR1 = UR11R*BR1 - UI11R*BI1 - UR12S*XR2 + UI12S*XI2;
            XI1 = UI11R*BR1 + UR11R*BI1 - UI12S*XR2 - UR12S*XI2;
            if ( ZSWAP( ICMAX ) ) {
               X( 1, 1 ) = XR2;
               X( 2, 1 ) = XR1;
               X( 1, 2 ) = XI2;
               X( 2, 2 ) = XI1;
            } else {
               X( 1, 1 ) = XR1;
               X( 2, 1 ) = XR2;
               X( 1, 2 ) = XI1;
               X( 2, 2 ) = XI2;
            }
            XNORM = MAX( ABS( XR1 )+ABS( XI1 ), ABS( XR2 )+ABS( XI2 ) );

            // Further scaling if  norm(A) norm(X) > overflow

            if ( XNORM > ONE && CMAX > ONE ) {
               if ( XNORM > BIGNUM / CMAX ) {
                  TEMP = CMAX / BIGNUM;
                  X( 1, 1 ) = TEMP*X( 1, 1 );
                  X( 2, 1 ) = TEMP*X( 2, 1 );
                  X( 1, 2 ) = TEMP*X( 1, 2 );
                  X( 2, 2 ) = TEMP*X( 2, 2 );
                  XNORM = TEMP*XNORM;
                  SCALE = TEMP*SCALE;
               }
            }
         }
      }

      return;

      // End of DLALN2

      }
