      SUBROUTINE ZGET08( TRANS, M, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

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
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                J, N1, N2;
      double             ANORM, BNORM, EPS, XNORM;
      COMPLEX*16         ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, IZAMAX, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Quick exit if M = 0 or N = 0 or NRHS = 0

      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.EQ.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      IF( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) THEN
         N1 = N
         N2 = M
      ELSE
         N1 = M
         N2 = N
      END IF

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANGE( 'I', N1, N2, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF

      // Compute  B - A*X  (or  B - A'*X ) and store in B.

      CALL ZGEMM( TRANS, 'No transpose', N1, NRHS, N2, -CONE, A, LDA, X, LDX, CONE, B, LDB )

      // Compute the maximum over the number of right hand sides of
         // norm(B - A*X) / ( norm(A) * norm(X) * EPS ) .

      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = CABS1( B( IZAMAX( N1, B( 1, J ), 1 ), J ) )
         XNORM = CABS1( X( IZAMAX( N2, X( 1, J ), 1 ), J ) )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE

      RETURN

      // End of ZGET08

      END
