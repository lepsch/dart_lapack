      void dptt01(N, D, E, DF, EF, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                N;
      double             RESID;
      double             D( * ), DF( * ), E( * ), EF( * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I;
      double             ANORM, DE, EPS;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX

      // Quick return if possible

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      EPS = dlamch( 'Epsilon' );

      // Construct the difference L*D*L' - A.

      WORK[1] = DF( 1 ) - D( 1 );
      for (I = 1; I <= N - 1; I++) { // 10
         DE = DF( I )*EF( I );
         WORK[N+I] = DE - E( I );
         WORK[1+I] = DE*EF( I ) + DF( I+1 ) - D( I+1 );
      } // 10

      // Compute the 1-norms of the tridiagonal matrices A and WORK.

      if ( N == 1 ) {
         ANORM = D( 1 );
         RESID = ( WORK( 1 ) ).abs();
      } else {
         ANORM = max( D( 1 )+( E( 1 ) ).abs(), D( N )+( E( N-1 ) ).abs() );
         RESID = max( ( WORK( 1 ) ).abs()+( WORK( N+1 ) ).abs(), ( WORK( N ) ).abs()+( WORK( 2*N-1 ) ).abs() );
         for (I = 2; I <= N - 1; I++) { // 20
            ANORM = max( ANORM, D( I )+( E( I ) ).abs()+( E( I-1 ) ).abs() );
            RESID = max( RESID, ( WORK( I ) ).abs()+( WORK( N+I-1 ) ).abs()+ ( WORK( N+I ) ).abs() );
         } // 20
      }

      // Compute norm(L*D*L' - A) / (n * norm(A) * EPS)

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / N.toDouble() ) / ANORM ) / EPS;
      }

      return;
      }
