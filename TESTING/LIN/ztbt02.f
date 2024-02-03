      SUBROUTINE ZTBT02( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, X, LDX, B, LDB, WORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         AB( LDAB, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DZASUM, ZLANTB;
      // EXTERNAL LSAME, DLAMCH, DZASUM, ZLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZTBMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // Compute the 1-norm of op(A).

      IF( LSAME( TRANS, 'N' ) ) THEN
         ANORM = ZLANTB( '1', UPLO, DIAG, N, KD, AB, LDAB, RWORK )
      ELSE
         ANORM = ZLANTB( 'I', UPLO, DIAG, N, KD, AB, LDAB, RWORK )
      END IF

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*x - b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO
      DO 10 J = 1, NRHS
         CALL ZCOPY( N, X( 1, J ), 1, WORK, 1 )
         CALL ZTBMV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 )
         CALL ZAXPY( N, DCMPLX( -ONE ), B( 1, J ), 1, WORK, 1 )
         BNORM = DZASUM( N, WORK, 1 )
         XNORM = DZASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE

      RETURN

      // End of ZTBT02

      END
