      SUBROUTINE SLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      REAL               SI1, SI2, SR1, SR2
      int                LDH, N;
      // ..
      // .. Array Arguments ..
      REAL               H( LDH, * ), V( * )
      // ..
*
*  ================================================================
*
      // .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0e0 )
      // ..
      // .. Local Scalars ..
      REAL               H21S, H31S, S
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..
*
      // Quick return if possible
*
      IF( N.NE.2 .AND. N.NE.3 ) THEN
         RETURN
      END IF
*
      IF( N.EQ.2 ) THEN
         S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) )
         IF( S.EQ.ZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            V( 1 ) = H21S*H( 1, 2 ) + ( H( 1, 1 )-SR1 )* ( ( H( 1, 1 )-SR2 ) / S ) - SI1*( SI2 / S )
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 )
         END IF
      ELSE
         S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) ) + ABS( H( 3, 1 ) )
         IF( S.EQ.ZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
            V( 3 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            H31S = H( 3, 1 ) / S
            V( 1 ) = ( H( 1, 1 )-SR1 )*( ( H( 1, 1 )-SR2 ) / S ) - SI1*( SI2 / S ) + H( 1, 2 )*H21S + H( 1, 3 )*H31S             V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 ) + H( 2, 3 )*H31S             V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-SR1-SR2 ) + H21S*H( 3, 2 )
         END IF
      END IF
      END
