      SUBROUTINE SLAED5( I, D, Z, DELTA, RHO, DLAM )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                I;
      REAL               DLAM, RHO
      // ..
      // .. Array Arguments ..
      REAL               D( 2 ), DELTA( 2 ), Z( 2 )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, FOUR
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, FOUR = 4.0E0 ;
      // ..
      // .. Local Scalars ..
      REAL               B, C, DEL, TAU, TEMP, W
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      DEL = D( 2 ) - D( 1 )
      if ( I == 1 ) {
         W = ONE + TWO*RHO*( Z( 2 )*Z( 2 )-Z( 1 )*Z( 1 ) ) / DEL
         if ( W.GT.ZERO ) {
            B = DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 1 )*Z( 1 )*DEL

            // B > ZERO, always

            TAU = TWO*C / ( B+SQRT( ABS( B*B-FOUR*C ) ) )
            DLAM = D( 1 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / TAU
            DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
         } else {
            B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 2 )*Z( 2 )*DEL
            if ( B.GT.ZERO ) {
               TAU = -TWO*C / ( B+SQRT( B*B+FOUR*C ) )
            } else {
               TAU = ( B-SQRT( B*B+FOUR*C ) ) / TWO
            }
            DLAM = D( 2 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
            DELTA( 2 ) = -Z( 2 ) / TAU
         }
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      } else {

      // Now I=2

         B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
         C = RHO*Z( 2 )*Z( 2 )*DEL
         if ( B.GT.ZERO ) {
            TAU = ( B+SQRT( B*B+FOUR*C ) ) / TWO
         } else {
            TAU = TWO*C / ( -B+SQRT( B*B+FOUR*C ) )
         }
         DLAM = D( 2 ) + TAU
         DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
         DELTA( 2 ) = -Z( 2 ) / TAU
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      }
      RETURN

      // End of SLAED5

      }
