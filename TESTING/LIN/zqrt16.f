      SUBROUTINE ZQRT16( TRANS, M, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDA, LDB, LDX, M, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                J, N1, N2;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DZASUM, ZLANGE;
      // EXTERNAL LSAME, DLAMCH, DZASUM, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if M = 0 or N = 0 or NRHS = 0

      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.EQ.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      IF( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) THEN
         ANORM = ZLANGE( 'I', M, N, A, LDA, RWORK )
         N1 = N
         N2 = M
      ELSE
         ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
         N1 = M
         N2 = N
      END IF

      EPS = DLAMCH( 'Epsilon' )

      // Compute  B - A*X  (or  B - A'*X ) and store in B.

      CALL ZGEMM( TRANS, 'No transpose', N1, NRHS, N2, -CONE, A, LDA, X, LDX, CONE, B, LDB )

      // Compute the maximum over the number of right hand sides of
         // norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ) .

      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = DZASUM( N1, B( 1, J ), 1 )
         XNORM = DZASUM( N2, X( 1, J ), 1 )
         IF( ANORM.EQ.ZERO .AND. BNORM.EQ.ZERO ) THEN
            RESID = ZERO
         ELSE IF( ANORM.LE.ZERO .OR. XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / ( MAX( M, N )*EPS ) )
         END IF
   10 CONTINUE

      RETURN

      // End of ZQRT16

      }
