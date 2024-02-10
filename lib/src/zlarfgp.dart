      void zlarfgp(N, ALPHA, X, INCX, TAU ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, N;
      Complex         ALPHA, TAU;
      Complex         X( * );
      // ..

      double             TWO, ONE, ZERO;
      const              TWO = 2.0, ONE = 1.0, ZERO = 0.0 ;
      int                J, KNT;
      double             ALPHI, ALPHR, BETA, BIGNUM, EPS, SMLNUM, XNORM;
      Complex         SAVEALPHA;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLAPY3, DLAPY2, DZNRM2;
      //- Complex         ZLADIV;
      // EXTERNAL DLAMCH, DLAPY3, DLAPY2, DZNRM2, ZLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DIMAG, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZSCAL

      if ( N <= 0 ) {
         TAU = ZERO;
         return;
      }

      EPS = dlamch( 'Precision' );
      XNORM = DZNRM2( N-1, X, INCX );
      ALPHR = ALPHA.toDouble();
      ALPHI = DIMAG( ALPHA );

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

         BETA = sign( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR );
         SMLNUM = dlamch( 'S' ) / dlamch( 'E' );
         BIGNUM = ONE / SMLNUM;

         KNT = 0;
         if ( ( BETA ).abs() < SMLNUM ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

            } // 10
            KNT = KNT + 1;
            zdscal(N-1, BIGNUM, X, INCX );
            BETA = BETA*BIGNUM;
            ALPHI = ALPHI*BIGNUM;
            ALPHR = ALPHR*BIGNUM;
            if( (( BETA ).abs() < SMLNUM) && (KNT < 20) ) GO TO 10;

            // New BETA is at most 1, at least SMLNUM

            XNORM = DZNRM2( N-1, X, INCX );
            ALPHA = DCMPLX( ALPHR, ALPHI );
            BETA = sign( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR );
         }
         SAVEALPHA = ALPHA;
         ALPHA = ALPHA + BETA;
         if ( BETA < ZERO ) {
            BETA = -BETA;
            TAU = -ALPHA / BETA;
         } else {
            ALPHR = ALPHI * (ALPHI/ALPHA.toDouble());
            ALPHR = ALPHR + XNORM * (XNORM/ALPHA.toDouble());
            TAU = DCMPLX( ALPHR/BETA, -ALPHI/BETA );
            ALPHA = DCMPLX( -ALPHR, ALPHI );
         }
         ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA );

         if ( (TAU).abs() <= SMLNUM ) {

            // In the case where the computed TAU ends up being a denormalized number,
            // it loses relative accuracy. This is a BIG problem. Solution: flush TAU
            // to ZERO (or TWO or whatever makes a nonnegative real number for BETA).

            // (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
            // (Thanks Pat. Thanks MathWorks.)

            ALPHR = SAVEALPHA.toDouble();
            ALPHI = DIMAG( SAVEALPHA );
            if ( ALPHI == ZERO ) {
               if ( ALPHR >= ZERO ) {
                  TAU = ZERO;
               } else {
                  TAU = TWO;
                  for (J = 1; J <= N-1; J++) {
                     X[1 + (J-1)*INCX] = ZERO;
                  }
                  BETA = (-SAVEALPHA).toDouble();
               }
            } else {
               XNORM = dlapy2( ALPHR, ALPHI );
               TAU = DCMPLX( ONE - ALPHR / XNORM, -ALPHI / XNORM );
               for (J = 1; J <= N-1; J++) {
                  X[1 + (J-1)*INCX] = ZERO;
               }
               BETA = XNORM;
            }

         } else {

            // This is the general case.

            zscal(N-1, ALPHA, X, INCX );

         }

         // If BETA is subnormal, it may lose relative accuracy

         for (J = 1; J <= KNT; J++) { // 20
            BETA = BETA*SMLNUM;
         } // 20
         ALPHA = BETA;
      }

      }
