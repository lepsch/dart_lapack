      void cget10(M, N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, WORK, final Array<double> RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, M, N;
      double               RESULT;
      double               RWORK( * );
      Complex            A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                J;
      double               ANORM, EPS, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- REAL               SCASUM, SLAMCH, CLANGE;
      // EXTERNAL SCASUM, SLAMCH, CLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      // Quick return if possible

      if ( M <= 0 || N <= 0 ) {
         RESULT = ZERO;
         return;
      }

      UNFL = SLAMCH( 'Safe minimum' );
      EPS = SLAMCH( 'Precision' );

      WNORM = ZERO;
      for (J = 1; J <= N; J++) { // 10
         ccopy(M, A( 1, J ), 1, WORK, 1 );
         caxpy(M, CMPLX( -ONE ), B( 1, J ), 1, WORK, 1 );
         WNORM = max( WNORM, SCASUM( N, WORK, 1 ) );
      } // 10

      ANORM = max( CLANGE( '1', M, N, A, LDA, RWORK ), UNFL );

      if ( ANORM > WNORM ) {
         RESULT = ( WNORM / ANORM ) / ( M*EPS );
      } else {
         if ( ANORM < ONE ) {
            RESULT = ( min( WNORM, M*ANORM ) / ANORM ) / ( M*EPS );
         } else {
            RESULT = min( WNORM / ANORM, REAL( M ) ) / ( M*EPS );
         }
      }

      }
