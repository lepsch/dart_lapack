      void cbdt02(M, N, final Matrix<double> B, final int LDB, final Matrix<double> C, final int LDC, final Matrix<double> U, final int LDU, WORK, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDB, LDC, LDU, M, N;
      double               RESID;
      double               RWORK( * );
      Complex            B( LDB, * ), C( LDC, * ), U( LDU, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               BNORM, EPS, REALMN;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SCASUM, SLAMCH;
      // EXTERNAL CLANGE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL

      // Quick return if possible

      RESID = ZERO;
      if (M <= 0 || N <= 0) return;
      REALMN = double( max( M, N ) );
      EPS = SLAMCH( 'Precision' );

      // Compute norm(B - U * C)

      for (J = 1; J <= N; J++) { // 10
         ccopy(M, B( 1, J ), 1, WORK, 1 );
         cgemv('No transpose', M, M, -CMPLX( ONE ), U, LDU, C( 1, J ), 1, CMPLX( ONE ), WORK, 1 );
         RESID = max( RESID, SCASUM( M, WORK, 1 ) );
      } // 10

      // Compute norm of B.

      BNORM = CLANGE( '1', M, N, B, LDB, RWORK );

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
