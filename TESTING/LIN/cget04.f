      SUBROUTINE CGET04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDX, LDXACT, N, NRHS;
      REAL               RCOND, RESID
      // ..
      // .. Array Arguments ..
      COMPLEX            X( LDX, * ), XACT( LDXACT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IX, J;
      REAL               DIFFNM, EPS, XNORM
      COMPLEX            ZDUM
      // ..
      // .. External Functions ..
      int                ICAMAX;
      REAL               SLAMCH
      // EXTERNAL ICAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N.LE.0 || NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if RCOND is invalid.

      EPS = SLAMCH( 'Epsilon' )
      if ( RCOND.LT.ZERO ) {
         RESID = 1.0 / EPS
         RETURN
      }

      // Compute the maximum of
         // norm(X - XACT) / ( norm(XACT) * EPS )
      // over all the vectors X and XACT .

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 20
         IX = ICAMAX( N, XACT( 1, J ), 1 )
         XNORM = CABS1( XACT( IX, J ) )
         DIFFNM = ZERO
         for (I = 1; I <= N; I++) { // 10
            DIFFNM = MAX( DIFFNM, CABS1( X( I, J )-XACT( I, J ) ) )
         } // 10
         if ( XNORM.LE.ZERO ) {
            if (DIFFNM.GT.ZERO) RESID = 1.0 / EPS;
         } else {
            RESID = MAX( RESID, ( DIFFNM / XNORM )*RCOND )
         }
      } // 20
      if (RESID*EPS.LT.1.0) RESID = RESID / EPS;

      RETURN

      // End of CGET04

      }
