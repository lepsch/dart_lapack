      void slarfg(N, ALPHA, X, INCX, TAU ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      REAL               ALPHA, TAU;
      // ..
      // .. Array Arguments ..
      REAL               X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                J, KNT;
      REAL               BETA, RSAFMN, SAFMIN, XNORM;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLAPY2, SNRM2;
      // EXTERNAL SLAMCH, SLAPY2, SNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL
      // ..
      // .. Executable Statements ..

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

         BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA );
         SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' );
         KNT = 0;
         if ( ABS( BETA ) < SAFMIN ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

            RSAFMN = ONE / SAFMIN;
            } // 10
            KNT = KNT + 1;
            sscal(N-1, RSAFMN, X, INCX );
            BETA = BETA*RSAFMN;
            ALPHA = ALPHA*RSAFMN;
            if( (ABS( BETA ) < SAFMIN) && (KNT < 20) ) GO TO 10;

            // New BETA is at most 1, at least SAFMIN

            XNORM = SNRM2( N-1, X, INCX );
            BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA );
         }
         TAU = ( BETA-ALPHA ) / BETA;
         sscal(N-1, ONE / ( ALPHA-BETA ), X, INCX );

         // If ALPHA is subnormal, it may lose relative accuracy

         for (J = 1; J <= KNT; J++) { // 20
            BETA = BETA*SAFMIN;
         } // 20
         ALPHA = BETA;
      }

      return;

      // End of SLARFG

      }
