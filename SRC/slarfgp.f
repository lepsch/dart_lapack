      SUBROUTINE SLARFGP( N, ALPHA, X, INCX, TAU );

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
      REAL               TWO, ONE, ZERO;
      const              TWO = 2.0, ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                J, KNT;
      REAL               BETA, BIGNUM, EPS, SAVEALPHA, SMLNUM, XNORM;
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

      if ( N <= 0 ) {
         TAU = ZERO;
         return;
      }

      EPS = SLAMCH( 'Precision' );
      XNORM = SNRM2( N-1, X, INCX );

      if ( XNORM <= EPS*ABS(ALPHA) ) {

         // H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0.

         if ( ALPHA >= ZERO ) {
            // When TAU == ZERO, the vector is special-cased to be
            // all zeros in the application routines.  We do not need
            // to clear it.
            TAU = ZERO;
         } else {
            // However, the application routines rely on explicit
            // zero checks when TAU != ZERO, and we must clear X.
            TAU = TWO;
            for (J = 1; J <= N-1; J++) {
               X( 1 + (J-1)*INCX ) = 0;
            }
            ALPHA = -ALPHA;
         }
      } else {

         // general case

         BETA = SIGN( SLAPY2( ALPHA, XNORM ), ALPHA );
         SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'E' );
         KNT = 0;
         if ( ABS( BETA ) < SMLNUM ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

            BIGNUM = ONE / SMLNUM;
            } // 10
            KNT = KNT + 1;
            sscal(N-1, BIGNUM, X, INCX );
            BETA = BETA*BIGNUM;
            ALPHA = ALPHA*BIGNUM;
            IF( (ABS( BETA ) < SMLNUM) && (KNT < 20) ) GO TO 10;

            // New BETA is at most 1, at least SMLNUM

            XNORM = SNRM2( N-1, X, INCX );
            BETA = SIGN( SLAPY2( ALPHA, XNORM ), ALPHA );
         }
         SAVEALPHA = ALPHA;
         ALPHA = ALPHA + BETA;
         if ( BETA < ZERO ) {
            BETA = -BETA;
            TAU = -ALPHA / BETA;
         } else {
            ALPHA = XNORM * (XNORM/ALPHA);
            TAU = ALPHA / BETA;
            ALPHA = -ALPHA;
         }

         if ( ABS(TAU) <= SMLNUM ) {

            // In the case where the computed TAU ends up being a denormalized number,
            // it loses relative accuracy. This is a BIG problem. Solution: flush TAU
            // to ZERO. This explains the next IF statement.

            // (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
            // (Thanks Pat. Thanks MathWorks.)

            if ( SAVEALPHA >= ZERO ) {
               TAU = ZERO;
            } else {
               TAU = TWO;
               for (J = 1; J <= N-1; J++) {
                  X( 1 + (J-1)*INCX ) = 0;
               }
               BETA = -SAVEALPHA;
            }

         } else {

            // This is the general case.

            sscal(N-1, ONE / ALPHA, X, INCX );

         }

         // If BETA is subnormal, it may lose relative accuracy

         for (J = 1; J <= KNT; J++) { // 20
            BETA = BETA*SMLNUM;
         } // 20
         ALPHA = BETA;
      }

      return;

      // End of SLARFGP

      }
