      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )

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
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                J, KNT;
      double             BETA, RSAFMN, SAFMIN, XNORM;
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

      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF

      XNORM = DNRM2( N-1, X, INCX )

      IF( XNORM.EQ.ZERO ) THEN

         // H  =  I

         TAU = ZERO
      ELSE

         // general case

         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN

            // XNORM, BETA may be inaccurate; scale X and recompute them

            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) ) GO TO 10

            // New BETA is at most 1, at least SAFMIN

            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )

         // If ALPHA is subnormal, it may lose relative accuracy

         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF

      RETURN

      // End of DLARFG

      }
