      void zbdt02(M, N, B, LDB, C, LDC, U, LDU, WORK, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDB, LDC, LDU, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      Complex         B( LDB, * ), C( LDC, * ), U( LDU, * ), WORK( * );
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
      //- double             DLAMCH, DZASUM, ZLANGE;
      // EXTERNAL DLAMCH, DZASUM, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      RESID = ZERO;
      if (M <= 0 || N <= 0) return;
      REALMN = (max( M, N )).toDouble();
      EPS = dlamch( 'Precision' );

      // Compute norm(B - U * C)

      for (J = 1; J <= N; J++) { // 10
         zcopy(M, B( 1, J ), 1, WORK, 1 );
         zgemv('No transpose', M, M, -DCMPLX( ONE ), U, LDU, C( 1, J ), 1, DCMPLX( ONE ), WORK, 1 );
         RESID = max( RESID, DZASUM( M, WORK, 1 ) );
      } // 10

      // Compute norm of B.

      BNORM = ZLANGE( '1', M, N, B, LDB, RWORK );

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
