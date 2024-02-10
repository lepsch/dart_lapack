      void dget03(N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK, RCOND, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDAINV, LDWORK, N;
      double             RCOND, RESID;
      double             A( LDA, * ), AINV( LDAINV, * ), RWORK( * ), WORK( LDWORK, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I;
      double             AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = dlamch( 'Epsilon' );
      ANORM = dlange( '1', N, N, A, LDA, RWORK );
      AINVNM = dlange( '1', N, N, AINV, LDAINV, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE / ANORM ) / AINVNM;

      // Compute I - A * AINV

      dgemm('No transpose', 'No transpose', N, N, N, -ONE, AINV, LDAINV, A, LDA, ZERO, WORK, LDWORK );
      for (I = 1; I <= N; I++) { // 10
         WORK[I][I] = ONE + WORK( I, I );
      } // 10

      // Compute norm(I - AINV*A) / (N * norm(A) * norm(AINV) * EPS)

      RESID = dlange( '1', N, N, WORK, LDWORK, RWORK );

      RESID = ( ( RESID*RCOND ) / EPS ) / N.toDouble();

      }
