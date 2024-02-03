      SUBROUTINE CQRT16( TRANS, M, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDA, LDB, LDX, M, N, NRHS;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                J, N1, N2;
      REAL               ANORM, BNORM, EPS, XNORM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SCASUM, SLAMCH
      // EXTERNAL LSAME, CLANGE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if M = 0 or N = 0 or NRHS = 0

      if ( M.LE.0 .OR. N.LE.0 .OR. NRHS == 0 ) {
         RESID = ZERO
         RETURN
      }

      if ( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) {
         ANORM = CLANGE( 'I', M, N, A, LDA, RWORK )
         N1 = N
         N2 = M
      } else {
         ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
         N1 = M
         N2 = N
      }

      EPS = SLAMCH( 'Epsilon' )

      // Compute  B - A*X  (or  B - A'*X ) and store in B.

      cgemm(TRANS, 'No transpose', N1, NRHS, N2, -CONE, A, LDA, X, LDX, CONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
         // norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ) .

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = SCASUM( N1, B( 1, J ), 1 )
         XNORM = SCASUM( N2, X( 1, J ), 1 )
         if ( ANORM == ZERO && BNORM == ZERO ) {
            RESID = ZERO
         } else if ( ANORM.LE.ZERO .OR. XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / ( MAX( M, N )*EPS ) )
         }
      } // 10

      RETURN

      // End of CQRT16

      }
