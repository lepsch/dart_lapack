      SUBROUTINE STPT02( UPLO, TRANS, DIAG, N, NRHS, AP, X, LDX, B, LDB, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, BNORM, EPS, XNORM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SASUM, SLAMCH, SLANTP
      // EXTERNAL LSAME, SASUM, SLAMCH, SLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, STPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // Compute the 1-norm of op(A).

      IF( LSAME( TRANS, 'N' ) ) THEN
         ANORM = SLANTP( '1', UPLO, DIAG, N, AP, WORK )
      ELSE
         ANORM = SLANTP( 'I', UPLO, DIAG, N, AP, WORK )
      END IF

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS ).

      RESID = ZERO
      DO 10 J = 1, NRHS
         CALL SCOPY( N, X( 1, J ), 1, WORK, 1 )
         CALL STPMV( UPLO, TRANS, DIAG, N, AP, WORK, 1 )
         CALL SAXPY( N, -ONE, B( 1, J ), 1, WORK, 1 )
         BNORM = SASUM( N, WORK, 1 )
         XNORM = SASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE

      RETURN

      // End of STPT02

      }
