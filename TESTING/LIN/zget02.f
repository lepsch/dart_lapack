      SUBROUTINE ZGET02( TRANS, M, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

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

      if ( M.LE.0 .OR. N.LE.0 .OR. NRHS.EQ.0 ) {
         RESID = ZERO
         RETURN
      }

      if ( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) {
         N1 = N
         N2 = M
      } else {
         N1 = M
         N2 = N
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      if ( LSAME( TRANS, 'N' ) ) {
         ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      } else {
         ANORM = ZLANGE( 'I', M, N, A, LDA, RWORK )
      }
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute B - op(A)*X and store in B.

      zgemm(TRANS, 'No transpose', N1, NRHS, N2, -CONE, A, LDA, X, LDX, CONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
         // norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ) .

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = DZASUM( N1, B( 1, J ), 1 )
         XNORM = DZASUM( N2, X( 1, J ), 1 )
         if ( XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 10

      RETURN

      // End of ZGET02

      }
