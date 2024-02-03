      double dlarmm(ANORM, BNORM, CNORM ) {
      // IMPLICIT NONE
      // .. Scalar Arguments ..
      double             ANORM, BNORM, CNORM;
      // .. Parameters ..
      double             ONE, HALF, FOUR;
      const              ONE = 1.0, HALF = 0.5, FOUR = 4.0 ;
      // ..
      // .. Local Scalars ..
       double             BIGNUM, SMLNUM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Executable Statements ..


      // Determine machine dependent parameters to control overflow.

      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' );
      BIGNUM = ( ONE / SMLNUM ) / FOUR;

      // Compute a scale factor.

      DLARMM = ONE;
      if ( BNORM <= ONE ) {
         if ( ANORM * BNORM > BIGNUM - CNORM ) {
            DLARMM = HALF;
         }
      } else {
         if ( ANORM > (BIGNUM - CNORM) / BNORM ) {
            DLARMM = HALF / BNORM;
         }
      }
      return;

      // ==== End of DLARMM ====

      }
