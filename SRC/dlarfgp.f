      SUBROUTINE DLARFGP( N, ALPHA, X, INCX, TAU )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      double             ALPHA, TAU;
      // ..
      // .. Array Arguments ..
      double             X( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             TWO, ONE, ZERO;
      const              TWO = 2.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                J, KNT;
      double             BETA, BIGNUM, EPS, SAVEALPHA, SMLNUM, XNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLAPY2, DNRM2;
      // EXTERNAL DLAMCH, DLAPY2, DNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL
      // ..
      // .. Executable Statements ..

      if ( N.LE.0 ) {
         TAU = ZERO
         RETURN
      }

      EPS = DLAMCH( 'Precision' )
      XNORM = DNRM2( N-1, X, INCX )

      if ( XNORM.LE.EPS*ABS(ALPHA) ) {

         // H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0.

         if ( ALPHA.GE.ZERO ) {
            // When TAU == ZERO, the vector is special-cased to be
            // all zeros in the application routines.  We do not need
            // to clear it.
            TAU = ZERO
         } else {
            // However, the application routines rely on explicit
            // zero checks when TAU.ne.ZERO, and we must clear X.
            TAU = TWO
            for (J = 1; J <= N-1; J++) {
               X( 1 + (J-1)*INCX ) = 0
            }
            ALPHA = -ALPHA
         }
      } else {

         // general case

         BETA = SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         if ( ABS( BETA ).LT.SMLNUM ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

            BIGNUM = ONE / SMLNUM
            } // 10
            KNT = KNT + 1
            dscal(N-1, BIGNUM, X, INCX );
            BETA = BETA*BIGNUM
            ALPHA = ALPHA*BIGNUM
            IF( (ABS( BETA ).LT.SMLNUM) .AND. (KNT .LT. 20) ) GO TO 10

            // New BETA is at most 1, at least SMLNUM

            XNORM = DNRM2( N-1, X, INCX )
            BETA = SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         }
         SAVEALPHA = ALPHA
         ALPHA = ALPHA + BETA
         if ( BETA.LT.ZERO ) {
            BETA = -BETA
            TAU = -ALPHA / BETA
         } else {
            ALPHA = XNORM * (XNORM/ALPHA)
            TAU = ALPHA / BETA
            ALPHA = -ALPHA
         }

         if ( ABS(TAU).LE.SMLNUM ) {

            // In the case where the computed TAU ends up being a denormalized number,
            // it loses relative accuracy. This is a BIG problem. Solution: flush TAU
            // to ZERO. This explains the next IF statement.

            // (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
            // (Thanks Pat. Thanks MathWorks.)

            if ( SAVEALPHA.GE.ZERO ) {
               TAU = ZERO
            } else {
               TAU = TWO
               for (J = 1; J <= N-1; J++) {
                  X( 1 + (J-1)*INCX ) = 0
               }
               BETA = -SAVEALPHA
            }

         } else {

            // This is the general case.

            dscal(N-1, ONE / ALPHA, X, INCX );

         }

         // If BETA is subnormal, it may lose relative accuracy

         for (J = 1; J <= KNT; J++) { // 20
            BETA = BETA*SMLNUM
         } // 20
         ALPHA = BETA
      }

      RETURN

      // End of DLARFGP

      }
