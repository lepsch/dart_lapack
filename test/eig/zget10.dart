      void zget10(M, N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Array<double> _WORK, final Array<double> RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, M, N;
      double             RESULT;
      double             RWORK( * );
      Complex         A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                J;
      double             ANORM, EPS, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DZASUM, ZLANGE;
      // EXTERNAL DLAMCH, DZASUM, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN

      // Quick return if possible

      if ( M <= 0 || N <= 0 ) {
         RESULT = ZERO;
         return;
      }

      UNFL = dlamch( 'Safe minimum' );
      EPS = dlamch( 'Precision' );

      WNORM = ZERO;
      for (J = 1; J <= N; J++) { // 10
         zcopy(M, A( 1, J ), 1, WORK, 1 );
         zaxpy(M, DCMPLX( -ONE ), B( 1, J ), 1, WORK, 1 );
         WNORM = max( WNORM, DZASUM( N, WORK, 1 ) );
      } // 10

      ANORM = max( ZLANGE( '1', M, N, A, LDA, RWORK ), UNFL );

      if ( ANORM > WNORM ) {
         RESULT = ( WNORM / ANORM ) / ( M*EPS );
      } else {
         if ( ANORM < ONE ) {
            RESULT = ( min( WNORM, M*ANORM ) / ANORM ) / ( M*EPS );
         } else {
            RESULT = min( WNORM / ANORM, M.toDouble() ) / ( M*EPS );
         }
      }

      }
