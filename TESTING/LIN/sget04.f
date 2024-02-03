      SUBROUTINE SGET04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDX, LDXACT, N, NRHS;
      REAL               RCOND, RESID
      // ..
      // .. Array Arguments ..
      REAL               X( LDX, * ), XACT( LDXACT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IX, J;
      REAL               DIFFNM, EPS, XNORM
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SLAMCH
      // EXTERNAL ISAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // Exit with RESID = 1/EPS if RCOND is invalid.

      EPS = SLAMCH( 'Epsilon' )
      IF( RCOND.LT.ZERO ) THEN
         RESID = 1.0 / EPS
         RETURN
      END IF

      // Compute the maximum of
         // norm(X - XACT) / ( norm(XACT) * EPS )
      // over all the vectors X and XACT .

      RESID = ZERO
      DO 20 J = 1, NRHS
         IX = ISAMAX( N, XACT( 1, J ), 1 )
         XNORM = ABS( XACT( IX, J ) )
         DIFFNM = ZERO
         DO 10 I = 1, N
            DIFFNM = MAX( DIFFNM, ABS( X( I, J )-XACT( I, J ) ) )
   10    CONTINUE
         IF( XNORM.LE.ZERO ) THEN
            IF( DIFFNM.GT.ZERO ) RESID = 1.0 / EPS
         } else {
            RESID = MAX( RESID, ( DIFFNM / XNORM )*RCOND )
         END IF
   20 CONTINUE
      IF( RESID*EPS.LT.1.0 ) RESID = RESID / EPS

      RETURN

      // End of SGET04

      }
