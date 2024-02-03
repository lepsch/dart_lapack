      SUBROUTINE ZTRT01( UPLO, DIAG, N, A, LDA, AINV, LDAINV, RCOND, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                LDA, LDAINV, N;
      double             RCOND, RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AINV( LDAINV, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANTR;
      // EXTERNAL LSAME, DLAMCH, ZLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0

      if ( N.LE.0 ) {
         RCOND = ONE
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANTR( '1', UPLO, DIAG, N, N, A, LDA, RWORK )
      AINVNM = ZLANTR( '1', UPLO, DIAG, N, N, AINV, LDAINV, RWORK )
      if ( ANORM.LE.ZERO || AINVNM.LE.ZERO ) {
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      }
      RCOND = ( ONE / ANORM ) / AINVNM

      // Set the diagonal of AINV to 1 if AINV has unit diagonal.

      if ( LSAME( DIAG, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 10
            AINV( J, J ) = ONE
         } // 10
      }

      // Compute A * AINV, overwriting AINV.

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 20
            ztrmv('Upper', 'No transpose', DIAG, J, A, LDA, AINV( 1, J ), 1 );
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 30
            ztrmv('Lower', 'No transpose', DIAG, N-J+1, A( J, J ), LDA, AINV( J, J ), 1 );
         } // 30
      }

      // Subtract 1 from each diagonal element to form A*AINV - I.

      for (J = 1; J <= N; J++) { // 40
         AINV( J, J ) = AINV( J, J ) - ONE
      } // 40

      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

      RESID = ZLANTR( '1', UPLO, 'Non-unit', N, N, AINV, LDAINV, RWORK )

      RESID = ( ( RESID*RCOND ) / DBLE( N ) ) / EPS

      RETURN

      // End of ZTRT01

      }
