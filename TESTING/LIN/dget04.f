      SUBROUTINE DGET04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDX, LDXACT, N, NRHS;
      double             RCOND, RESID;
      // ..
      // .. Array Arguments ..
      double             X( LDX, * ), XACT( LDXACT, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IX, J;
      double             DIFFNM, EPS, XNORM;
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL IDAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N.LE.0 .OR. NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if RCOND is invalid.

      EPS = DLAMCH( 'Epsilon' )
      if ( RCOND.LT.ZERO ) {
         RESID = 1.0D0 / EPS
         RETURN
      }

      // Compute the maximum of
         // norm(X - XACT) / ( norm(XACT) * EPS )
      // over all the vectors X and XACT .

      RESID = ZERO
      DO 20 J = 1, NRHS
         IX = IDAMAX( N, XACT( 1, J ), 1 )
         XNORM = ABS( XACT( IX, J ) )
         DIFFNM = ZERO
         DO 10 I = 1, N
            DIFFNM = MAX( DIFFNM, ABS( X( I, J )-XACT( I, J ) ) )
   10    CONTINUE
         if ( XNORM.LE.ZERO ) {
            IF( DIFFNM.GT.ZERO ) RESID = 1.0D0 / EPS
         } else {
            RESID = MAX( RESID, ( DIFFNM / XNORM )*RCOND )
         }
   20 CONTINUE
      IF( RESID*EPS.LT.1.0D0 ) RESID = RESID / EPS

      RETURN

      // End of DGET04

      }
