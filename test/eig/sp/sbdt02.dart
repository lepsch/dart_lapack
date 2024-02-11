      void sbdt02(final int M, final int N, final Matrix<double> B, final int LDB, final Matrix<double> C, final int LDC, final Matrix<double> U, final int LDU, final Array<double> _WORK, final int RESID,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDB, LDC, LDU, M, N;
      double               RESID;
      double               B( LDB, * ), C( LDC, * ), U( LDU, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               BNORM, EPS, REALMN;
      // ..
      // .. External Functions ..
      //- REAL               SASUM, SLAMCH, SLANGE;
      // EXTERNAL SASUM, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      // Quick return if possible

      RESID = ZERO;
      if (M <= 0 || N <= 0) return;
      REALMN = double( max( M, N ) );
      EPS = SLAMCH( 'Precision' );

      // Compute norm(B - U * C)

      for (J = 1; J <= N; J++) { // 10
         scopy(M, B( 1, J ), 1, WORK, 1 );
         sgemv('No transpose', M, M, -ONE, U, LDU, C( 1, J ), 1, ONE, WORK, 1 );
         RESID = max( RESID, SASUM( M, WORK, 1 ) );
      } // 10

      // Compute norm of B.

      BNORM = SLANGE( '1', M, N, B, LDB, WORK );

      if ( BNORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( BNORM >= RESID ) {
            RESID = ( RESID / BNORM ) / ( REALMN*EPS );
         } else {
            if ( BNORM < ONE ) {
               RESID = ( min( RESID, REALMN*BNORM ) / BNORM ) / ( REALMN*EPS );
            } else {
               RESID = min( RESID / BNORM, REALMN ) / ( REALMN*EPS );
            }
         }
      }
      }
