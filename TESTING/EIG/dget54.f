      SUBROUTINE DGET54( N, A, LDA, B, LDB, S, LDS, T, LDT, U, LDU, V, LDV, WORK, RESULT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDS, LDT, LDU, LDV, N;
      double             RESULT;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), S( LDS, * ), T( LDT, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double             ABNORM, ULP, UNFL, WNORM;
      // ..
      // .. Local Arrays ..
      double             DUM( 1 );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      RESULT = ZERO;
      if (N <= 0) return;

      // Constants

      UNFL = DLAMCH( 'Safe minimum' );
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' );

      // compute the norm of (A,B)

      dlacpy('Full', N, N, A, LDA, WORK, N );
      dlacpy('Full', N, N, B, LDB, WORK( N*N+1 ), N );
      ABNORM = max( DLANGE( '1', N, 2*N, WORK, N, DUM ), UNFL );

      // Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)

      dlacpy(' ', N, N, A, LDA, WORK, N );
      dgemm('N', 'N', N, N, N, ONE, U, LDU, S, LDS, ZERO, WORK( N*N+1 ), N );

      dgemm('N', 'C', N, N, N, -ONE, WORK( N*N+1 ), N, V, LDV, ONE, WORK, N );

      // Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)

      dlacpy(' ', N, N, B, LDB, WORK( N*N+1 ), N );
      dgemm('N', 'N', N, N, N, ONE, U, LDU, T, LDT, ZERO, WORK( 2*N*N+1 ), N );

      dgemm('N', 'C', N, N, N, -ONE, WORK( 2*N*N+1 ), N, V, LDV, ONE, WORK( N*N+1 ), N );

      // Compute norm(W)/ ( ulp*norm((A,B)) )

      WNORM = DLANGE( '1', N, 2*N, WORK, N, DUM );

      if ( ABNORM > WNORM ) {
         RESULT = ( WNORM / ABNORM ) / ( 2*N*ULP );
      } else {
         if ( ABNORM < ONE ) {
            RESULT = ( min( WNORM, 2*N*ABNORM ) / ABNORM ) / ( 2*N*ULP );
         } else {
            RESULT = min( WNORM / ABNORM, DBLE( 2*N ) ) / ( 2*N*ULP );
         }
      }

      return;

      // End of DGET54

      }
