      void zlarfg(N, ALPHA, X, INCX, TAU ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      Complex         ALPHA, TAU;
      // ..
      // .. Array Arguments ..
      Complex         X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                J, KNT;
      double             ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLAPY3, DZNRM2;
      Complex         ZLADIV;
      // EXTERNAL DLAMCH, DLAPY3, DZNRM2, ZLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DIMAG, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZSCAL
      // ..
      // .. Executable Statements ..

      if ( N <= 0 ) {
         TAU = ZERO;
         return;
      }

      XNORM = DZNRM2( N-1, X, INCX );
      ALPHR = DBLE( ALPHA );
      ALPHI = DIMAG( ALPHA );

      if ( XNORM == ZERO && ALPHI == ZERO ) {

         // H  =  I

         TAU = ZERO;
      } else {

         // general case

         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR );
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' );
         RSAFMN = ONE / SAFMIN;

         KNT = 0;
         if ( ABS( BETA ) < SAFMIN ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

            } // 10
            KNT = KNT + 1;
            zdscal(N-1, RSAFMN, X, INCX );
            BETA = BETA*RSAFMN;
            ALPHI = ALPHI*RSAFMN;
            ALPHR = ALPHR*RSAFMN;
            if( (ABS( BETA ) < SAFMIN) && (KNT < 20) ) GO TO 10;

            // New BETA is at most 1, at least SAFMIN

            XNORM = DZNRM2( N-1, X, INCX );
            ALPHA = DCMPLX( ALPHR, ALPHI );
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR );
         }
         TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA );
         ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA );
         zscal(N-1, ALPHA, X, INCX );

         // If ALPHA is subnormal, it may lose relative accuracy

         for (J = 1; J <= KNT; J++) { // 20
            BETA = BETA*SAFMIN;
         } // 20
         ALPHA = BETA;
      }

      return;

      // End of ZLARFG

      }
