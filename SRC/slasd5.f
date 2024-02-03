      void slasd5(I, D, Z, DELTA, RHO, DSIGMA, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                I;
      REAL               DSIGMA, RHO;
      // ..
      // .. Array Arguments ..
      REAL               D( 2 ), DELTA( 2 ), WORK( 2 ), Z( 2 );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, THREE, FOUR;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0, FOUR = 4.0 ;
      // ..
      // .. Local Scalars ..
      REAL               B, C, DEL, DELSQ, TAU, W;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      DEL = D( 2 ) - D( 1 );
      DELSQ = DEL*( D( 2 )+D( 1 ) );
      if ( I == 1 ) {
         W = ONE + FOUR*RHO*( Z( 2 )*Z( 2 ) / ( D( 1 )+THREE*D( 2 ) )- Z( 1 )*Z( 1 ) / ( THREE*D( 1 )+D( 2 ) ) ) / DEL;
         if ( W > ZERO ) {
            B = DELSQ + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) );
            C = RHO*Z( 1 )*Z( 1 )*DELSQ;

            // B > ZERO, always

            // The following TAU is DSIGMA * DSIGMA - D( 1 ) * D( 1 )

            TAU = TWO*C / ( B+sqrt( ABS( B*B-FOUR*C ) ) );

            // The following TAU is DSIGMA - D( 1 )

            TAU = TAU / ( D( 1 )+sqrt( D( 1 )*D( 1 )+TAU ) );
            DSIGMA = D( 1 ) + TAU;
            DELTA( 1 ) = -TAU;
            DELTA( 2 ) = DEL - TAU;
            WORK( 1 ) = TWO*D( 1 ) + TAU;
            WORK( 2 ) = ( D( 1 )+TAU ) + D( 2 );
            // DELTA( 1 ) = -Z( 1 ) / TAU
            // DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
         } else {
            B = -DELSQ + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) );
            C = RHO*Z( 2 )*Z( 2 )*DELSQ;

            // The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )

            if ( B > ZERO ) {
               TAU = -TWO*C / ( B+sqrt( B*B+FOUR*C ) );
            } else {
               TAU = ( B-sqrt( B*B+FOUR*C ) ) / TWO;
            }

            // The following TAU is DSIGMA - D( 2 )

            TAU = TAU / ( D( 2 )+sqrt( ABS( D( 2 )*D( 2 )+TAU ) ) );
            DSIGMA = D( 2 ) + TAU;
            DELTA( 1 ) = -( DEL+TAU );
            DELTA( 2 ) = -TAU;
            WORK( 1 ) = D( 1 ) + TAU + D( 2 );
            WORK( 2 ) = TWO*D( 2 ) + TAU;
            // DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
            // DELTA( 2 ) = -Z( 2 ) / TAU
         }
         // TEMP = sqrt( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         // DELTA( 1 ) = DELTA( 1 ) / TEMP
         // DELTA( 2 ) = DELTA( 2 ) / TEMP
      } else {

         // Now I=2

         B = -DELSQ + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) );
         C = RHO*Z( 2 )*Z( 2 )*DELSQ;

         // The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )

         if ( B > ZERO ) {
            TAU = ( B+sqrt( B*B+FOUR*C ) ) / TWO;
         } else {
            TAU = TWO*C / ( -B+sqrt( B*B+FOUR*C ) );
         }

         // The following TAU is DSIGMA - D( 2 )

         TAU = TAU / ( D( 2 )+sqrt( D( 2 )*D( 2 )+TAU ) );
         DSIGMA = D( 2 ) + TAU;
         DELTA( 1 ) = -( DEL+TAU );
         DELTA( 2 ) = -TAU;
         WORK( 1 ) = D( 1 ) + TAU + D( 2 );
         WORK( 2 ) = TWO*D( 2 ) + TAU;
         // DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
         // DELTA( 2 ) = -Z( 2 ) / TAU
         // TEMP = sqrt( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         // DELTA( 1 ) = DELTA( 1 ) / TEMP
         // DELTA( 2 ) = DELTA( 2 ) / TEMP
      }
      return;
      }
