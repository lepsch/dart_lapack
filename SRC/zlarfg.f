      SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )

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
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      int                J, KNT;
      double             ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLAPY3, DZNRM2;
      COMPLEX*16         ZLADIV
      // EXTERNAL DLAMCH, DLAPY3, DZNRM2, ZLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DIMAG, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZSCAL
      // ..
      // .. Executable Statements ..

      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF

      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )

      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN

         // H  =  I

         TAU = ZERO
      ELSE

         // general case

         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN

         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN

            // XNORM, BETA may be inaccurate; scale X and recompute them

   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) ) GO TO 10

            // New BETA is at most 1, at least SAFMIN

            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         END IF
         TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
         ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
         CALL ZSCAL( N-1, ALPHA, X, INCX )

         // If ALPHA is subnormal, it may lose relative accuracy

         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF

      RETURN

      // End of ZLARFG

      }
