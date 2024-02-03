      REAL FUNCTION SLARMM( ANORM, BNORM, CNORM )
      IMPLICIT NONE
      // .. Scalar Arguments ..
      REAL               ANORM, BNORM, CNORM
      // .. Parameters ..
      REAL               ONE, HALF, FOUR
      const              ONE = 1.0E0, HALF = 0.5E+0, FOUR = 4.0E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               BIGNUM, SMLNUM
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Executable Statements ..


      // Determine machine dependent parameters to control overflow.

      SMLNUM = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
      BIGNUM = ( ONE / SMLNUM ) / FOUR

      // Compute a scale factor.

      SLARMM = ONE
      if ( BNORM <= ONE ) {
         if ( ANORM * BNORM > BIGNUM - CNORM ) {
            SLARMM = HALF
         }
      } else {
         if ( ANORM > (BIGNUM - CNORM) / BNORM ) {
            SLARMM = HALF / BNORM
         }
      }
      RETURN

      // ==== End of SLARMM ====

      }
