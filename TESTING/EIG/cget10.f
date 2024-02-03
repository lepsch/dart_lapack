      SUBROUTINE CGET10( M, N, A, LDA, B, LDB, WORK, RWORK, RESULT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, M, N;
      REAL               RESULT;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, EPS, UNFL, WNORM;
      // ..
      // .. External Functions ..
      REAL               SCASUM, SLAMCH, CLANGE;
      // EXTERNAL SCASUM, SLAMCH, CLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( M <= 0 || N <= 0 ) {
         RESULT = ZERO;
         RETURN;
      }

      UNFL = SLAMCH( 'Safe minimum' );
      EPS = SLAMCH( 'Precision' );

      WNORM = ZERO;
      for (J = 1; J <= N; J++) { // 10
         ccopy(M, A( 1, J ), 1, WORK, 1 );
         caxpy(M, CMPLX( -ONE ), B( 1, J ), 1, WORK, 1 );
         WNORM = MAX( WNORM, SCASUM( N, WORK, 1 ) );
      } // 10

      ANORM = MAX( CLANGE( '1', M, N, A, LDA, RWORK ), UNFL );

      if ( ANORM > WNORM ) {
         RESULT = ( WNORM / ANORM ) / ( M*EPS );
      } else {
         if ( ANORM < ONE ) {
            RESULT = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*EPS );
         } else {
            RESULT = MIN( WNORM / ANORM, REAL( M ) ) / ( M*EPS );
         }
      }

      RETURN;

      // End of CGET10

      }
