      void cget54(N, A, LDA, B, LDB, S, LDS, T, LDT, U, LDU, V, LDV, WORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDS, LDT, LDU, LDV, N;
      double               RESULT;
      // ..
      // .. Array Arguments ..
      Complex            A( LDA, * ), B( LDB, * ), S( LDS, * ), T( LDT, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      double               ABNORM, ULP, UNFL, WNORM;
      // ..
      // .. Local Arrays ..
      double               DUM( 1 );
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      RESULT = ZERO;
      if (N <= 0) return;

      // Constants

      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );

      // compute the norm of (A,B)

      clacpy('Full', N, N, A, LDA, WORK, N );
      clacpy('Full', N, N, B, LDB, WORK( N*N+1 ), N );
      ABNORM = max( CLANGE( '1', N, 2*N, WORK, N, DUM ), UNFL );

      // Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)

      clacpy(' ', N, N, A, LDA, WORK, N );
      cgemm('N', 'N', N, N, N, CONE, U, LDU, S, LDS, CZERO, WORK( N*N+1 ), N );

      cgemm('N', 'C', N, N, N, -CONE, WORK( N*N+1 ), N, V, LDV, CONE, WORK, N );

      // Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)

      clacpy(' ', N, N, B, LDB, WORK( N*N+1 ), N );
      cgemm('N', 'N', N, N, N, CONE, U, LDU, T, LDT, CZERO, WORK( 2*N*N+1 ), N );

      cgemm('N', 'C', N, N, N, -CONE, WORK( 2*N*N+1 ), N, V, LDV, CONE, WORK( N*N+1 ), N );

      // Compute norm(W)/ ( ulp*norm((A,B)) )

      WNORM = CLANGE( '1', N, 2*N, WORK, N, DUM );

      if ( ABNORM > WNORM ) {
         RESULT = ( WNORM / ABNORM ) / ( 2*N*ULP );
      } else {
         if ( ABNORM < ONE ) {
            RESULT = ( min( WNORM, 2*N*ABNORM ) / ABNORM ) / ( 2*N*ULP );
         } else {
            RESULT = min( WNORM / ABNORM, REAL( 2*N ) ) / ( 2*N*ULP );
         }
      }

      return;
      }
