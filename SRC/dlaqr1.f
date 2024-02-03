      SUBROUTINE DLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             SI1, SI2, SR1, SR2;
      int                LDH, N;
      // ..
      // .. Array Arguments ..
      double             H( LDH, * ), V( * );
      // ..

*  ================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0d0 ;
      // ..
      // .. Local Scalars ..
      double             H21S, H31S, S;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N != 2 .AND. N != 3 ) {
         RETURN
      }

      if ( N == 2 ) {
         S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) )
         if ( S == ZERO ) {
            V( 1 ) = ZERO
            V( 2 ) = ZERO
         } else {
            H21S = H( 2, 1 ) / S
            V( 1 ) = H21S*H( 1, 2 ) + ( H( 1, 1 )-SR1 )* ( ( H( 1, 1 )-SR2 ) / S ) - SI1*( SI2 / S )
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 )
         }
      } else {
         S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) ) + ABS( H( 3, 1 ) )
         if ( S == ZERO ) {
            V( 1 ) = ZERO
            V( 2 ) = ZERO
            V( 3 ) = ZERO
         } else {
            H21S = H( 2, 1 ) / S
            H31S = H( 3, 1 ) / S
            V( 1 ) = ( H( 1, 1 )-SR1 )*( ( H( 1, 1 )-SR2 ) / S ) - SI1*( SI2 / S ) + H( 1, 2 )*H21S + H( 1, 3 )*H31S             V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 ) + H( 2, 3 )*H31S             V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-SR1-SR2 ) + H21S*H( 3, 2 )
         }
      }
      }
