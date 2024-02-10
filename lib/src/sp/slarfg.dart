      void slarfg(N, ALPHA, X, INCX, final int TAU) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, N;
      double               ALPHA, TAU;
      double               X( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                J, KNT;
      double               BETA, RSAFMN, SAFMIN, XNORM;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLAPY2, SNRM2;
      // EXTERNAL SLAMCH, SLAPY2, SNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL

      if ( N <= 1 ) {
         TAU = ZERO;
         return;
      }

      XNORM = SNRM2( N-1, X, INCX );

      if ( XNORM == ZERO ) {

         // H  =  I

         TAU = ZERO;
      } else {

         // general case

         BETA = -sign( SLAPY2( ALPHA, XNORM ), ALPHA );
         SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' );
         KNT = 0;
         if ( ( BETA ).abs() < SAFMIN ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

            RSAFMN = ONE / SAFMIN;
            } // 10
            KNT = KNT + 1;
            sscal(N-1, RSAFMN, X, INCX );
            BETA = BETA*RSAFMN;
            ALPHA = ALPHA*RSAFMN;
            if( (( BETA ).abs() < SAFMIN) && (KNT < 20) ) GO TO 10;

            // New BETA is at most 1, at least SAFMIN

            XNORM = SNRM2( N-1, X, INCX );
            BETA = -sign( SLAPY2( ALPHA, XNORM ), ALPHA );
         }
         TAU = ( BETA-ALPHA ) / BETA;
         sscal(N-1, ONE / ( ALPHA-BETA ), X, INCX );

         // If ALPHA is subnormal, it may lose relative accuracy

         for (J = 1; J <= KNT; J++) { // 20
            BETA = BETA*SAFMIN;
         } // 20
         ALPHA = BETA;
      }

      }
