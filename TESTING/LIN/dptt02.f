      SUBROUTINE DPTT02( N, NRHS, D, E, X, LDX, B, LDB, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), D( * ), E( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANST;
      // EXTERNAL DASUM, DLAMCH, DLANST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAPTM
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // Compute the 1-norm of the tridiagonal matrix A.

      ANORM = DLANST( '1', N, D, E )

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF

      // Compute B - A*X.

      CALL DLAPTM( N, NRHS, -ONE, D, E, X, LDX, ONE, B, LDB )

      // Compute the maximum over the number of right hand sides of
         // norm(B - A*X) / ( norm(A) * norm(X) * EPS ).

      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = DASUM( N, B( 1, J ), 1 )
         XNORM = DASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE

      RETURN

      // End of DPTT02

      }
