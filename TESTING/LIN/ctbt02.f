      SUBROUTINE CTBT02( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, X, LDX, B, LDB, WORK, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      int                KD, LDAB, LDB, LDX, N, NRHS
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            AB( LDAB, * ), B( LDB, * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      int                J
      REAL               ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANTB, SCASUM, SLAMCH
      EXTERNAL           LSAME, CLANTB, SCASUM, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CCOPY, CTBMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Compute the 1-norm of op(A).
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         ANORM = CLANTB( '1', UPLO, DIAG, N, KD, AB, LDAB, RWORK )
      ELSE
         ANORM = CLANTB( 'I', UPLO, DIAG, N, KD, AB, LDAB, RWORK )
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = SLAMCH( 'Epsilon' )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute the maximum over the number of right hand sides of
*        norm(op(A)*x - b) / ( norm(op(A)) * norm(x) * EPS ).
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         CALL CCOPY( N, X( 1, J ), 1, WORK, 1 )
         CALL CTBMV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 )
         CALL CAXPY( N, CMPLX( -ONE ), B( 1, J ), 1, WORK, 1 )
         BNORM = SCASUM( N, WORK, 1 )
         XNORM = SCASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of CTBT02
*
      END
