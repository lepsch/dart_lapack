      void dbdt02(M, N, B, LDB, C, LDC, U, LDU, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDB, LDC, LDU, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), C( LDC, * ), U( LDU, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             BNORM, EPS, REALMN;
      // ..
      // .. External Functions ..
      //- double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      RESID = ZERO;
      if (M <= 0 || N <= 0) return;
      REALMN = (max( M, N )).toDouble();
      EPS = dlamch( 'Precision' );

      // Compute norm(B - U * C)

      for (J = 1; J <= N; J++) { // 10
         dcopy(M, B( 1, J ), 1, WORK, 1 );
         dgemv('No transpose', M, M, -ONE, U, LDU, C( 1, J ), 1, ONE, WORK, 1 );
         RESID = max( RESID, DASUM( M, WORK, 1 ) );
      } // 10

      // Compute norm of B.

      BNORM = DLANGE( '1', M, N, B, LDB, WORK );

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
      return;
      }
