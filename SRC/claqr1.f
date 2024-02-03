      SUBROUTINE CLAQR1( N, H, LDH, S1, S2, V )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX            S1, S2
      int                LDH, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            H( LDH, * ), V( * )
      // ..

*  ================================================================

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0e0, 0.0e0 ) ;
      REAL               RZERO
      const              RZERO = 0.0e0 ;
      // ..
      // .. Local Scalars ..
      COMPLEX            CDUM, H21S, H31S
      REAL               S
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.NE.2 .AND. N.NE.3 ) {
         RETURN
      }

      if ( N == 2 ) {
         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) )
         if ( S == RZERO ) {
            V( 1 ) = ZERO
            V( 2 ) = ZERO
         } else {
            H21S = H( 2, 1 ) / S
            V( 1 ) = H21S*H( 1, 2 ) + ( H( 1, 1 )-S1 )* ( ( H( 1, 1 )-S2 ) / S )
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 )
         }
      } else {
         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) ) + CABS1( H( 3, 1 ) )
         if ( S == ZERO ) {
            V( 1 ) = ZERO
            V( 2 ) = ZERO
            V( 3 ) = ZERO
         } else {
            H21S = H( 2, 1 ) / S
            H31S = H( 3, 1 ) / S
            V( 1 ) = ( H( 1, 1 )-S1 )*( ( H( 1, 1 )-S2 ) / S ) + H( 1, 2 )*H21S + H( 1, 3 )*H31S
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 ) + H( 2, 3 )*H31S
            V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-S1-S2 ) + H21S*H( 3, 2 )
         }
      }
      }
