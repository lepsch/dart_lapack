      SUBROUTINE CLARFGP( N, ALPHA, X, INCX, TAU )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      COMPLEX            ALPHA, TAU
      // ..
      // .. Array Arguments ..
      COMPLEX            X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               TWO, ONE, ZERO
      const              TWO = 2.0E+0, ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                J, KNT;
      REAL               ALPHI, ALPHR, BETA, BIGNUM, EPS, SMLNUM, XNORM
      COMPLEX            SAVEALPHA
      // ..
      // .. External Functions ..
      REAL               SCNRM2, SLAMCH, SLAPY3, SLAPY2
      COMPLEX            CLADIV
      // EXTERNAL SCNRM2, SLAMCH, SLAPY3, SLAPY2, CLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, REAL, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSCAL, CSSCAL
      // ..
      // .. Executable Statements ..

      if ( N.LE.0 ) {
         TAU = ZERO
         RETURN
      }

      EPS = SLAMCH( 'Precision' )
      XNORM = SCNRM2( N-1, X, INCX )
      ALPHR = REAL( ALPHA )
      ALPHI = AIMAG( ALPHA )

      if ( XNORM.LE.EPS*ABS(ALPHA) .AND. ALPHI.EQ.ZERO ) {

         // H  =  [1-alpha/abs(alpha) 0; 0 I], sign chosen so ALPHA >= 0.

         if ( ALPHR.GE.ZERO ) {
            // When TAU.eq.ZERO, the vector is special-cased to be
            // all zeros in the application routines.  We do not need
            // to clear it.
            TAU = ZERO
         } else {
            // However, the application routines rely on explicit
            // zero checks when TAU.ne.ZERO, and we must clear X.
            TAU = TWO
            DO J = 1, N-1
               X( 1 + (J-1)*INCX ) = ZERO
            END DO
            ALPHA = -ALPHA
         }
      } else {

         // general case

         BETA = SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'E' )
         BIGNUM = ONE / SMLNUM

         KNT = 0
         if ( ABS( BETA ).LT.SMLNUM ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

   10       CONTINUE
            KNT = KNT + 1
            CALL CSSCAL( N-1, BIGNUM, X, INCX )
            BETA = BETA*BIGNUM
            ALPHI = ALPHI*BIGNUM
            ALPHR = ALPHR*BIGNUM
            IF( (ABS( BETA ).LT.SMLNUM) .AND. (KNT .LT. 20) ) GO TO 10

            // New BETA is at most 1, at least SMLNUM

            XNORM = SCNRM2( N-1, X, INCX )
            ALPHA = CMPLX( ALPHR, ALPHI )
            BETA = SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         }
         SAVEALPHA = ALPHA
         ALPHA = ALPHA + BETA
         if ( BETA.LT.ZERO ) {
            BETA = -BETA
            TAU = -ALPHA / BETA
         } else {
            ALPHR = ALPHI * (ALPHI/REAL( ALPHA ))
            ALPHR = ALPHR + XNORM * (XNORM/REAL( ALPHA ))
            TAU = CMPLX( ALPHR/BETA, -ALPHI/BETA )
            ALPHA = CMPLX( -ALPHR, ALPHI )
         }
         ALPHA = CLADIV( CMPLX( ONE ), ALPHA )

         if ( ABS(TAU).LE.SMLNUM ) {

            // In the case where the computed TAU ends up being a denormalized number,
            // it loses relative accuracy. This is a BIG problem. Solution: flush TAU
            // to ZERO (or TWO or whatever makes a nonnegative real number for BETA).

            // (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
            // (Thanks Pat. Thanks MathWorks.)

            ALPHR = REAL( SAVEALPHA )
            ALPHI = AIMAG( SAVEALPHA )
            if ( ALPHI.EQ.ZERO ) {
               if ( ALPHR.GE.ZERO ) {
                  TAU = ZERO
               } else {
                  TAU = TWO
                  DO J = 1, N-1
                     X( 1 + (J-1)*INCX ) = ZERO
                  END DO
                  BETA = REAL( -SAVEALPHA )
               }
            } else {
               XNORM = SLAPY2( ALPHR, ALPHI )
               TAU = CMPLX( ONE - ALPHR / XNORM, -ALPHI / XNORM )
               DO J = 1, N-1
                  X( 1 + (J-1)*INCX ) = ZERO
               END DO
               BETA = XNORM
            }

         } else {

            // This is the general case.

            CALL CSCAL( N-1, ALPHA, X, INCX )

         }

         // If BETA is subnormal, it may lose relative accuracy

         DO 20 J = 1, KNT
            BETA = BETA*SMLNUM
 20      CONTINUE
         ALPHA = BETA
      }

      RETURN

      // End of CLARFGP

      }
