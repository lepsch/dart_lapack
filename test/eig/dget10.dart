      void dget10(M, N, A, LDA, B, LDB, WORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, M, N;
      double             RESULT;
      double             A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                J;
      double             ANORM, EPS, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN

      // Quick return if possible

      if ( M <= 0 || N <= 0 ) {
         RESULT = ZERO;
         return;
      }

      UNFL = dlamch( 'Safe minimum' );
      EPS = dlamch( 'Precision' );

      WNORM = ZERO;
      for (J = 1; J <= N; J++) { // 10
         dcopy(M, A( 1, J ), 1, WORK, 1 );
         daxpy(M, -ONE, B( 1, J ), 1, WORK, 1 );
         WNORM = max( WNORM, dasum( N, WORK, 1 ) );
      } // 10

      ANORM = max( dlange( '1', M, N, A, LDA, WORK ), UNFL );

      if ( ANORM > WNORM ) {
         RESULT = ( WNORM / ANORM ) / ( M*EPS );
      } else {
         if ( ANORM < ONE ) {
            RESULT = ( min( WNORM, M*ANORM ) / ANORM ) / ( M*EPS );
         } else {
            RESULT = min( WNORM / ANORM, M.toDouble() ) / ( M*EPS );
         }
      }

      return;
      }
