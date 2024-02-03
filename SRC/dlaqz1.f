      SUBROUTINE DLAQZ1( A, LDA, B, LDB, SR1, SR2, SI, BETA1, BETA2, V )
      IMPLICIT NONE

      // Arguments
      int    , INTENT( IN ) :: LDA, LDB;
      double          , INTENT( IN ) :: A( LDA, * ), B( LDB, * ), SR1, SR2, SI, BETA1, BETA2;
      double          , INTENT( OUT ) :: V( * );

      // Parameters
      double           :: ZERO, ONE, HALF;
      const    ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 ;

      // Local scalars
      double           :: W( 2 ), SAFMIN, SAFMAX, SCALE1, SCALE2;

      // External Functions
      double          , EXTERNAL :: DLAMCH;
      bool   , EXTERNAL :: DISNAN;

      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE/SAFMIN

      // Calculate first shifted vector

      W( 1 ) = BETA1*A( 1, 1 )-SR1*B( 1, 1 )
      W( 2 ) = BETA1*A( 2, 1 )-SR1*B( 2, 1 )
      SCALE1 = SQRT( ABS( W( 1 ) ) ) * SQRT( ABS( W( 2 ) ) )
      if ( SCALE1 .GE. SAFMIN && SCALE1 .LE. SAFMAX ) {
         W( 1 ) = W( 1 )/SCALE1
         W( 2 ) = W( 2 )/SCALE1
      }

      // Solve linear system

      W( 2 ) = W( 2 )/B( 2, 2 )
      W( 1 ) = ( W( 1 )-B( 1, 2 )*W( 2 ) )/B( 1, 1 )
      SCALE2 = SQRT( ABS( W( 1 ) ) ) * SQRT( ABS( W( 2 ) ) )
      if ( SCALE2 .GE. SAFMIN && SCALE2 .LE. SAFMAX ) {
         W( 1 ) = W( 1 )/SCALE2
         W( 2 ) = W( 2 )/SCALE2
      }

      // Apply second shift

      V( 1 ) = BETA2*( A( 1, 1 )*W( 1 )+A( 1, 2 )*W( 2 ) )-SR2*( B( 1, 1 )*W( 1 )+B( 1, 2 )*W( 2 ) )       V( 2 ) = BETA2*( A( 2, 1 )*W( 1 )+A( 2, 2 )*W( 2 ) )-SR2*( B( 2, 1 )*W( 1 )+B( 2, 2 )*W( 2 ) )       V( 3 ) = BETA2*( A( 3, 1 )*W( 1 )+A( 3, 2 )*W( 2 ) )-SR2*( B( 3, 1 )*W( 1 )+B( 3, 2 )*W( 2 ) )

      // Account for imaginary part

      V( 1 ) = V( 1 )+SI*SI*B( 1, 1 )/SCALE1/SCALE2

      // Check for overflow

      if ( ABS( V( 1 ) ).GT.SAFMAX .OR. ABS( V( 2 ) ) .GT. SAFMAX .OR. ABS( V( 3 ) ).GT.SAFMAX .OR. DISNAN( V( 1 ) ) .OR. DISNAN( V( 2 ) ) .OR. DISNAN( V( 3 ) ) ) {
         V( 1 ) = ZERO
         V( 2 ) = ZERO
         V( 3 ) = ZERO
      }

      // End of DLAQZ1

      END SUBROUTINE
