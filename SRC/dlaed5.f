      SUBROUTINE DLAED5( I, D, Z, DELTA, RHO, DLAM )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                I;
      double             DLAM, RHO;
      // ..
      // .. Array Arguments ..
      double             D( 2 ), DELTA( 2 ), Z( 2 );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, FOUR;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, FOUR = 4.0D0 )
      // ..
      // .. Local Scalars ..
      double             B, C, DEL, TAU, TEMP, W;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      DEL = D( 2 ) - D( 1 )
      IF( I.EQ.1 ) THEN
         W = ONE + TWO*RHO*( Z( 2 )*Z( 2 )-Z( 1 )*Z( 1 ) ) / DEL
         IF( W.GT.ZERO ) THEN
            B = DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 1 )*Z( 1 )*DEL

            // B > ZERO, always

            TAU = TWO*C / ( B+SQRT( ABS( B*B-FOUR*C ) ) )
            DLAM = D( 1 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / TAU
            DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
         ELSE
            B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 2 )*Z( 2 )*DEL
            IF( B.GT.ZERO ) THEN
               TAU = -TWO*C / ( B+SQRT( B*B+FOUR*C ) )
            ELSE
               TAU = ( B-SQRT( B*B+FOUR*C ) ) / TWO
            END IF
            DLAM = D( 2 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
            DELTA( 2 ) = -Z( 2 ) / TAU
         END IF
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      ELSE

      // Now I=2

         B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
         C = RHO*Z( 2 )*Z( 2 )*DEL
         IF( B.GT.ZERO ) THEN
            TAU = ( B+SQRT( B*B+FOUR*C ) ) / TWO
         ELSE
            TAU = TWO*C / ( -B+SQRT( B*B+FOUR*C ) )
         END IF
         DLAM = D( 2 ) + TAU
         DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
         DELTA( 2 ) = -Z( 2 ) / TAU
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      END IF
      RETURN

      // End of DLAED5

      END
