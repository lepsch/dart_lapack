      SUBROUTINE ZLARFGP( N, ALPHA, X, INCX, TAU )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      COMPLEX*16         ALPHA, TAU
      // ..
      // .. Array Arguments ..
      COMPLEX*16         X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             TWO, ONE, ZERO;
      const              TWO = 2.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                J, KNT;
      double             ALPHI, ALPHR, BETA, BIGNUM, EPS, SMLNUM, XNORM;
      COMPLEX*16         SAVEALPHA
      // ..
      // .. External Functions ..
      double             DLAMCH, DLAPY3, DLAPY2, DZNRM2;
      COMPLEX*16         ZLADIV
      // EXTERNAL DLAMCH, DLAPY3, DLAPY2, DZNRM2, ZLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DIMAG, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZSCAL
      // ..
      // .. Executable Statements ..

      if ( N.LE.0 ) {
         TAU = ZERO
         RETURN
      }

      EPS = DLAMCH( 'Precision' )
      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )

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

         BETA = SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'E' )
         BIGNUM = ONE / SMLNUM

         KNT = 0
         if ( ABS( BETA ).LT.SMLNUM ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

   10       CONTINUE
            KNT = KNT + 1
            zdscal(N-1, BIGNUM, X, INCX );
            BETA = BETA*BIGNUM
            ALPHI = ALPHI*BIGNUM
            ALPHR = ALPHR*BIGNUM
            IF( (ABS( BETA ).LT.SMLNUM) .AND. (KNT .LT. 20) ) GO TO 10

            // New BETA is at most 1, at least SMLNUM

            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         }
         SAVEALPHA = ALPHA
         ALPHA = ALPHA + BETA
         if ( BETA.LT.ZERO ) {
            BETA = -BETA
            TAU = -ALPHA / BETA
         } else {
            ALPHR = ALPHI * (ALPHI/DBLE( ALPHA ))
            ALPHR = ALPHR + XNORM * (XNORM/DBLE( ALPHA ))
            TAU = DCMPLX( ALPHR/BETA, -ALPHI/BETA )
            ALPHA = DCMPLX( -ALPHR, ALPHI )
         }
         ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA )

         if ( ABS(TAU).LE.SMLNUM ) {

            // In the case where the computed TAU ends up being a denormalized number,
            // it loses relative accuracy. This is a BIG problem. Solution: flush TAU
            // to ZERO (or TWO or whatever makes a nonnegative real number for BETA).

            // (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
            // (Thanks Pat. Thanks MathWorks.)

            ALPHR = DBLE( SAVEALPHA )
            ALPHI = DIMAG( SAVEALPHA )
            if ( ALPHI.EQ.ZERO ) {
               if ( ALPHR.GE.ZERO ) {
                  TAU = ZERO
               } else {
                  TAU = TWO
                  DO J = 1, N-1
                     X( 1 + (J-1)*INCX ) = ZERO
                  END DO
                  BETA = DBLE( -SAVEALPHA )
               }
            } else {
               XNORM = DLAPY2( ALPHR, ALPHI )
               TAU = DCMPLX( ONE - ALPHR / XNORM, -ALPHI / XNORM )
               DO J = 1, N-1
                  X( 1 + (J-1)*INCX ) = ZERO
               END DO
               BETA = XNORM
            }

         } else {

            // This is the general case.

            zscal(N-1, ALPHA, X, INCX );

         }

         // If BETA is subnormal, it may lose relative accuracy

         DO 20 J = 1, KNT
            BETA = BETA*SMLNUM
 20      CONTINUE
         ALPHA = BETA
      }

      RETURN

      // End of ZLARFGP

      }
