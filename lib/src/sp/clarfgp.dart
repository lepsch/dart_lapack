      void clarfgp(final int N, final int ALPHA, final int X, final int INCX, final int TAU) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, N;
      Complex            ALPHA, TAU;
      Complex            X( * );
      // ..

      double               TWO, ONE, ZERO;
      const              TWO = 2.0, ONE = 1.0, ZERO = 0.0 ;
      int                J, KNT;
      double               ALPHI, ALPHR, BETA, BIGNUM, EPS, SMLNUM, XNORM;
      Complex            SAVEALPHA;
      // ..
      // .. External Functions ..
      //- REAL               SCNRM2, SLAMCH, SLAPY3, SLAPY2;
      //- COMPLEX            CLADIV;
      // EXTERNAL SCNRM2, SLAMCH, SLAPY3, SLAPY2, CLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, REAL, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSCAL, CSSCAL

      if ( N <= 0 ) {
         TAU = ZERO;
         return;
      }

      EPS = SLAMCH( 'Precision' );
      XNORM = SCNRM2( N-1, X, INCX );
      ALPHR = double( ALPHA );
      ALPHI = AIMAG( ALPHA );

      if ( XNORM <= EPS*(ALPHA).abs() && ALPHI == ZERO ) {

         // H  =  [1-alpha/abs(alpha) 0; 0 I], sign chosen so ALPHA >= 0.

         if ( ALPHR >= ZERO ) {
            // When TAU == ZERO, the vector is special-cased to be
            // all zeros in the application routines.  We do not need
            // to clear it.
            TAU = ZERO;
         } else {
            // However, the application routines rely on explicit
            // zero checks when TAU != ZERO, and we must clear X.
            TAU = TWO;
            for (J = 1; J <= N-1; J++) {
               X[1 + (J-1)*INCX] = ZERO;
            }
            ALPHA = -ALPHA;
         }
      } else {

         // general case

         BETA = sign( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR );
         SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'E' );
         BIGNUM = ONE / SMLNUM;

         KNT = 0;
         if ( ( BETA ).abs() < SMLNUM ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

            } // 10
            KNT = KNT + 1;
            csscal(N-1, BIGNUM, X, INCX );
            BETA = BETA*BIGNUM;
            ALPHI = ALPHI*BIGNUM;
            ALPHR = ALPHR*BIGNUM;
            if( (( BETA ).abs() < SMLNUM) && (KNT < 20) ) GO TO 10;

            // New BETA is at most 1, at least SMLNUM

            XNORM = SCNRM2( N-1, X, INCX );
            ALPHA = CMPLX( ALPHR, ALPHI );
            BETA = sign( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR );
         }
         SAVEALPHA = ALPHA;
         ALPHA = ALPHA + BETA;
         if ( BETA < ZERO ) {
            BETA = -BETA;
            TAU = -ALPHA / BETA;
         } else {
            ALPHR = ALPHI * (ALPHI/REAL( ALPHA ));
            ALPHR = ALPHR + XNORM * (XNORM/REAL( ALPHA ));
            TAU = CMPLX( ALPHR/BETA, -ALPHI/BETA );
            ALPHA = CMPLX( -ALPHR, ALPHI );
         }
         ALPHA = CLADIV( CMPLX( ONE ), ALPHA );

         if ( (TAU).abs() <= SMLNUM ) {

            // In the case where the computed TAU ends up being a denormalized number,
            // it loses relative accuracy. This is a BIG problem. Solution: flush TAU
            // to ZERO (or TWO or whatever makes a nonnegative real number for BETA).

            // (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
            // (Thanks Pat. Thanks MathWorks.)

            ALPHR = double( SAVEALPHA );
            ALPHI = AIMAG( SAVEALPHA );
            if ( ALPHI == ZERO ) {
               if ( ALPHR >= ZERO ) {
                  TAU = ZERO;
               } else {
                  TAU = TWO;
                  for (J = 1; J <= N-1; J++) {
                     X[1 + (J-1)*INCX] = ZERO;
                  }
                  BETA = double( -SAVEALPHA );
               }
            } else {
               XNORM = SLAPY2( ALPHR, ALPHI );
               TAU = CMPLX( ONE - ALPHR / XNORM, -ALPHI / XNORM );
               for (J = 1; J <= N-1; J++) {
                  X[1 + (J-1)*INCX] = ZERO;
               }
               BETA = XNORM;
            }

         } else {

            // This is the general case.

            cscal(N-1, ALPHA, X, INCX );

         }

         // If BETA is subnormal, it may lose relative accuracy

         for (J = 1; J <= KNT; J++) { // 20
            BETA = BETA*SMLNUM;
         } // 20
         ALPHA = BETA;
      }

      }
