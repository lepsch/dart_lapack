      void sget03(N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK, RCOND, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDAINV, LDWORK, N;
      REAL               RCOND, RESID;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), AINV( LDAINV, * ), RWORK( * ), WORK( LDWORK, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = SLANGE( '1', N, N, A, LDA, RWORK );
      AINVNM = SLANGE( '1', N, N, AINV, LDAINV, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE / ANORM ) / AINVNM;

      // Compute I - A * AINV

      sgemm('No transpose', 'No transpose', N, N, N, -ONE, AINV, LDAINV, A, LDA, ZERO, WORK, LDWORK );
      for (I = 1; I <= N; I++) { // 10
         WORK( I, I ) = ONE + WORK( I, I );
      } // 10

      // Compute norm(I - AINV*A) / (N * norm(A) * norm(AINV) * EPS)

      RESID = SLANGE( '1', N, N, WORK, LDWORK, RWORK );

      RESID = ( ( RESID*RCOND ) / EPS ) / REAL( N );

      return;
      }
