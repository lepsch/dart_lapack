      void zget54(N, A, LDA, B, LDB, S, LDS, T, LDT, U, LDU, V, LDV, WORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LDS, LDT, LDU, LDV, N;
      double             RESULT;
      Complex         A( LDA, * ), B( LDB, * ), S( LDS, * ), T( LDT, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      double             ABNORM, ULP, UNFL, WNORM;
      double             DUM( 1 );
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN

      RESULT = ZERO;
      if (N <= 0) return;

      // Constants

      UNFL = dlamch( 'Safe minimum' );
      ULP = dlamch( 'Epsilon' )*dlamch( 'Base' );

      // compute the norm of (A,B)

      zlacpy('Full', N, N, A, LDA, WORK, N );
      zlacpy('Full', N, N, B, LDB, WORK( N*N+1 ), N );
      ABNORM = max( ZLANGE( '1', N, 2*N, WORK, N, DUM ), UNFL );

      // Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)

      zlacpy(' ', N, N, A, LDA, WORK, N );
      zgemm('N', 'N', N, N, N, CONE, U, LDU, S, LDS, CZERO, WORK( N*N+1 ), N );

      zgemm('N', 'C', N, N, N, -CONE, WORK( N*N+1 ), N, V, LDV, CONE, WORK, N );

      // Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)

      zlacpy(' ', N, N, B, LDB, WORK( N*N+1 ), N );
      zgemm('N', 'N', N, N, N, CONE, U, LDU, T, LDT, CZERO, WORK( 2*N*N+1 ), N );

      zgemm('N', 'C', N, N, N, -CONE, WORK( 2*N*N+1 ), N, V, LDV, CONE, WORK( N*N+1 ), N );

      // Compute norm(W)/ ( ulp*norm((A,B)) )

      WNORM = ZLANGE( '1', N, 2*N, WORK, N, DUM );

      if ( ABNORM > WNORM ) {
         RESULT = ( WNORM / ABNORM ) / ( 2*N*ULP );
      } else {
         if ( ABNORM < ONE ) {
            RESULT = ( min( WNORM, 2*N*ABNORM ) / ABNORM ) / ( 2*N*ULP );
         } else {
            RESULT = min( WNORM / ABNORM, (2*N).toDouble() ) / ( 2*N*ULP );
         }
      }

      return;
      }
