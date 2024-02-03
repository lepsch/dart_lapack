      void slaqz1(A, LDA, B, LDB, SR1, SR2, SI, BETA1, BETA2, V ) {
      // IMPLICIT NONE

      // Arguments
      int    , INTENT( IN ) :: LDA, LDB;
      REAL, INTENT( IN ) :: A( LDA, * ), B( LDB, * ), SR1, SR2, SI, BETA1, BETA2;
      REAL, INTENT( OUT ) :: V( * );

      // Parameters
      REAL :: ZERO, ONE, HALF;
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local scalars
      REAL :: W( 2 ), SAFMIN, SAFMAX, SCALE1, SCALE2;

      // External Functions
      REAL, EXTERNAL :: SLAMCH;
      bool   , EXTERNAL :: SISNAN;

      SAFMIN = SLAMCH( 'SAFE MINIMUM' );
      SAFMAX = ONE/SAFMIN;

      // Calculate first shifted vector

      W( 1 ) = BETA1*A( 1, 1 )-SR1*B( 1, 1 );
      W( 2 ) = BETA1*A( 2, 1 )-SR1*B( 2, 1 );
      SCALE1 = sqrt( ABS( W( 1 ) ) ) * sqrt( ABS( W( 2 ) ) );
      if ( SCALE1 >= SAFMIN && SCALE1 <= SAFMAX ) {
         W( 1 ) = W( 1 )/SCALE1;
         W( 2 ) = W( 2 )/SCALE1;
      }

      // Solve linear system

      W( 2 ) = W( 2 )/B( 2, 2 );
      W( 1 ) = ( W( 1 )-B( 1, 2 )*W( 2 ) )/B( 1, 1 );
      SCALE2 = sqrt( ABS( W( 1 ) ) ) * sqrt( ABS( W( 2 ) ) );
      if ( SCALE2 >= SAFMIN && SCALE2 <= SAFMAX ) {
         W( 1 ) = W( 1 )/SCALE2;
         W( 2 ) = W( 2 )/SCALE2;
      }

      // Apply second shift

      V( 1 ) = BETA2*( A( 1, 1 )*W( 1 )+A( 1, 2 )*W( 2 ) )-SR2*( B( 1, 1 )*W( 1 )+B( 1, 2 )*W( 2 ) )       V( 2 ) = BETA2*( A( 2, 1 )*W( 1 )+A( 2, 2 )*W( 2 ) )-SR2*( B( 2, 1 )*W( 1 )+B( 2, 2 )*W( 2 ) )       V( 3 ) = BETA2*( A( 3, 1 )*W( 1 )+A( 3, 2 )*W( 2 ) )-SR2*( B( 3, 1 )*W( 1 )+B( 3, 2 )*W( 2 ) );

      // Account for imaginary part

      V( 1 ) = V( 1 )+SI*SI*B( 1, 1 )/SCALE1/SCALE2;

      // Check for overflow

      if ( ABS( V( 1 ) ) > SAFMAX || ABS( V( 2 ) ) > SAFMAX || ABS( V( 3 ) ) > SAFMAX || SISNAN( V( 1 ) ) || SISNAN( V( 2 ) ) || SISNAN( V( 3 ) ) ) {
         V( 1 ) = ZERO;
         V( 2 ) = ZERO;
         V( 3 ) = ZERO;
      }

      // End of SLAQZ1

      END SUBROUTINE;
