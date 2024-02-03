      double           FUNCTION DLARMM( ANORM, BNORM, CNORM );
      IMPLICIT NONE
*     .. Scalar Arguments ..
      double             ANORM, BNORM, CNORM;
*     .. Parameters ..
      double             ONE, HALF, FOUR;
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D+0, FOUR = 4.0D0 )
*     ..
*     .. Local Scalars ..
       double             BIGNUM, SMLNUM;
*     ..
*     .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
*     ..
*     .. Executable Statements ..
*
*
*     Determine machine dependent parameters to control overflow.
*
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ( ONE / SMLNUM ) / FOUR
*
*     Compute a scale factor.
*
      DLARMM = ONE
      IF( BNORM .LE. ONE ) THEN
         IF( ANORM * BNORM .GT. BIGNUM - CNORM ) THEN
            DLARMM = HALF
         END IF
      ELSE
         IF( ANORM .GT. (BIGNUM - CNORM) / BNORM ) THEN
            DLARMM = HALF / BNORM
         END IF
      END IF
      RETURN
*
*     ==== End of DLARMM ====
*
      END
