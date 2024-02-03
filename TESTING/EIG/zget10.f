      SUBROUTINE ZGET10( M, N, A, LDA, B, LDB, WORK, RWORK, RESULT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, M, N;
      double             RESULT;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, EPS, UNFL, WNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DZASUM, ZLANGE;
      // EXTERNAL DLAMCH, DZASUM, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( M <= 0 || N <= 0 ) {
         RESULT = ZERO;
         return;
      }

      UNFL = DLAMCH( 'Safe minimum' );
      EPS = DLAMCH( 'Precision' );

      WNORM = ZERO;
      for (J = 1; J <= N; J++) { // 10
         zcopy(M, A( 1, J ), 1, WORK, 1 );
         zaxpy(M, DCMPLX( -ONE ), B( 1, J ), 1, WORK, 1 );
         WNORM = MAX( WNORM, DZASUM( N, WORK, 1 ) );
      } // 10

      ANORM = MAX( ZLANGE( '1', M, N, A, LDA, RWORK ), UNFL );

      if ( ANORM > WNORM ) {
         RESULT = ( WNORM / ANORM ) / ( M*EPS );
      } else {
         if ( ANORM < ONE ) {
            RESULT = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*EPS );
         } else {
            RESULT = MIN( WNORM / ANORM, DBLE( M ) ) / ( M*EPS );
         }
      }

      return;

      // End of ZGET10

      }
