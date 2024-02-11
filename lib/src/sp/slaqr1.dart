      void slaqr1(final int N, final Matrix<double> H, final int LDH, final int SR1, final int SI1, final int SR2, final int SI2, final int V,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               SI1, SI2, SR1, SR2;
      int                LDH, N;
      double               H( LDH, * ), V( * );
      // ..

// ================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      double               H21S, H31S, S;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      // Quick return if possible

      if ( N != 2 && N != 3 ) {
         return;
      }

      if ( N == 2 ) {
         S = ABS( H( 1, 1 )-SR2 ) + ( SI2 ).abs() + ( H( 2, 1 ) ).abs();
         if ( S == ZERO ) {
            V[1] = ZERO;
            V[2] = ZERO;
         } else {
            H21S = H( 2, 1 ) / S;
            V[1] = H21S*H( 1, 2 ) + ( H( 1, 1 )-SR1 )* ( ( H( 1, 1 )-SR2 ) / S ) - SI1*( SI2 / S );
            V[2] = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 );
         }
      } else {
         S = ABS( H( 1, 1 )-SR2 ) + ( SI2 ).abs() + ( H( 2, 1 ) ).abs() + ( H( 3, 1 ) ).abs();
         if ( S == ZERO ) {
            V[1] = ZERO;
            V[2] = ZERO;
            V[3] = ZERO;
         } else {
            H21S = H( 2, 1 ) / S;
            H31S = H( 3, 1 ) / S;
            V[1] = ( H( 1, 1 )-SR1 )*( ( H( 1, 1 )-SR2 ) / S ) - SI1*( SI2 / S ) + H( 1, 2 )*H21S + H( 1, 3 )*H31S             V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 ) + H( 2, 3 )*H31S             V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-SR1-SR2 ) + H21S*H( 3, 2 );
         }
      }
      }
