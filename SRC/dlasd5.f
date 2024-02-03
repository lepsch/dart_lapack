      SUBROUTINE DLASD5( I, D, Z, DELTA, RHO, DSIGMA, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                I;
      double             DSIGMA, RHO;
      // ..
      // .. Array Arguments ..
      double             D( 2 ), DELTA( 2 ), WORK( 2 ), Z( 2 );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, THREE, FOUR;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0, THREE = 3.0D+0, FOUR = 4.0D+0 ;
      // ..
      // .. Local Scalars ..
      double             B, C, DEL, DELSQ, TAU, W;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      DEL = D( 2 ) - D( 1 )
      DELSQ = DEL*( D( 2 )+D( 1 ) )
      if ( I == 1 ) {
         W = ONE + FOUR*RHO*( Z( 2 )*Z( 2 ) / ( D( 1 )+THREE*D( 2 ) )- Z( 1 )*Z( 1 ) / ( THREE*D( 1 )+D( 2 ) ) ) / DEL
         if ( W > ZERO ) {
            B = DELSQ + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 1 )*Z( 1 )*DELSQ

            // B > ZERO, always

            // The following TAU is DSIGMA * DSIGMA - D( 1 ) * D( 1 )

            TAU = TWO*C / ( B+SQRT( ABS( B*B-FOUR*C ) ) )

            // The following TAU is DSIGMA - D( 1 )

            TAU = TAU / ( D( 1 )+SQRT( D( 1 )*D( 1 )+TAU ) )
            DSIGMA = D( 1 ) + TAU
            DELTA( 1 ) = -TAU
            DELTA( 2 ) = DEL - TAU
            WORK( 1 ) = TWO*D( 1 ) + TAU
            WORK( 2 ) = ( D( 1 )+TAU ) + D( 2 )
            // DELTA( 1 ) = -Z( 1 ) / TAU
            // DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
         } else {
            B = -DELSQ + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 2 )*Z( 2 )*DELSQ

            // The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )

            if ( B > ZERO ) {
               TAU = -TWO*C / ( B+SQRT( B*B+FOUR*C ) )
            } else {
               TAU = ( B-SQRT( B*B+FOUR*C ) ) / TWO
            }

            // The following TAU is DSIGMA - D( 2 )

            TAU = TAU / ( D( 2 )+SQRT( ABS( D( 2 )*D( 2 )+TAU ) ) )
            DSIGMA = D( 2 ) + TAU
            DELTA( 1 ) = -( DEL+TAU )
            DELTA( 2 ) = -TAU
            WORK( 1 ) = D( 1 ) + TAU + D( 2 )
            WORK( 2 ) = TWO*D( 2 ) + TAU
            // DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
            // DELTA( 2 ) = -Z( 2 ) / TAU
         }
         // TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         // DELTA( 1 ) = DELTA( 1 ) / TEMP
         // DELTA( 2 ) = DELTA( 2 ) / TEMP
      } else {

         // Now I=2

         B = -DELSQ + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
         C = RHO*Z( 2 )*Z( 2 )*DELSQ

         // The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )

         if ( B > ZERO ) {
            TAU = ( B+SQRT( B*B+FOUR*C ) ) / TWO
         } else {
            TAU = TWO*C / ( -B+SQRT( B*B+FOUR*C ) )
         }

         // The following TAU is DSIGMA - D( 2 )

         TAU = TAU / ( D( 2 )+SQRT( D( 2 )*D( 2 )+TAU ) )
         DSIGMA = D( 2 ) + TAU
         DELTA( 1 ) = -( DEL+TAU )
         DELTA( 2 ) = -TAU
         WORK( 1 ) = D( 1 ) + TAU + D( 2 )
         WORK( 2 ) = TWO*D( 2 ) + TAU
         // DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
         // DELTA( 2 ) = -Z( 2 ) / TAU
         // TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         // DELTA( 1 ) = DELTA( 1 ) / TEMP
         // DELTA( 2 ) = DELTA( 2 ) / TEMP
      }
      RETURN

      // End of DLASD5

      }
