      SUBROUTINE ZPOT06( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double             RESID;
*     ..
*     .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CONE, NEGCONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
      PARAMETER          ( NEGCONE = ( -1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      int                IFAIL, J;
      double             ANORM, BNORM, EPS, XNORM;
      COMPLEX*16         ZDUM
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH, ZLANSY;
      EXTERNAL           LSAME, IZAMAX, DLAMCH, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZHEMM
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
*     ..
*     .. Statement Functions ..
      double             CABS1;
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0
*
      IF( N.LE.0 .OR. NRHS.EQ.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANSY( 'I', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute  B - A*X  and store in B.
      IFAIL=0
*
      CALL ZHEMM( 'Left', UPLO, N, NRHS, NEGCONE, A, LDA, X, LDX, CONE, B, LDB )
*
*     Compute the maximum over the number of right hand sides of
*        norm(B - A*X) / ( norm(A) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = CABS1(B(IZAMAX( N, B( 1, J ), 1 ),J))
         XNORM = CABS1(X(IZAMAX( N, X( 1, J ), 1 ),J))
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of ZPOT06
*
      END
