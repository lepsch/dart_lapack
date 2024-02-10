      void cget03(N, final Matrix<double> A, final int LDA, final Matrix<double> AINV, final int LDAINV, final Matrix<double> WORK, final int LDWORK, final Array<double> RWORK, RCOND, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDAINV, LDWORK, N;
      double               RCOND, RESID;
      double               RWORK( * );
      Complex            A( LDA, * ), AINV( LDAINV, * ), WORK( LDWORK, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I;
      double               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANGE( '1', N, N, A, LDA, RWORK );
      AINVNM = CLANGE( '1', N, N, AINV, LDAINV, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE/ANORM ) / AINVNM;

      // Compute I - A * AINV

      cgemm('No transpose', 'No transpose', N, N, N, -CONE, AINV, LDAINV, A, LDA, CZERO, WORK, LDWORK );
      for (I = 1; I <= N; I++) { // 10
         WORK[I][I] = CONE + WORK( I, I );
      } // 10

      // Compute norm(I - AINV*A) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANGE( '1', N, N, WORK, LDWORK, RWORK );

      RESID = ( ( RESID*RCOND )/EPS ) / REAL( N );

      }
