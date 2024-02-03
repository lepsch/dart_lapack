      SUBROUTINE ZHET01_AA( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE  = ( 1.0, 0.0 ) ;
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHE;
      // EXTERNAL LSAME, DLAMCH, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASET, ZLAVHE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' );
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK );

      // Initialize C to the tridiagonal matrix T.

      zlaset('Full', N, N, CZERO, CZERO, C, LDC );
      zlacpy('F', 1, N, AFAC( 1, 1 ), LDAFAC+1, C( 1, 1 ), LDC+1 );
      if ( N > 1 ) {
         if ( LSAME( UPLO, 'U' ) ) {
            zlacpy('F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 1, 2 ), LDC+1 );
            zlacpy('F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 2, 1 ), LDC+1 );
            zlacgv(N-1, C( 2, 1 ), LDC+1 );
         } else {
            zlacpy('F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 1, 2 ), LDC+1 );
            zlacpy('F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 2, 1 ), LDC+1 );
            zlacgv(N-1, C( 1, 2 ), LDC+1 );
         }

         // Call ZTRMM to form the product U' * D (or L * D ).

         if ( LSAME( UPLO, 'U' ) ) {
            ztrmm('Left', UPLO, 'Conjugate transpose', 'Unit', N-1, N, CONE, AFAC( 1, 2 ), LDAFAC, C( 2, 1 ), LDC );
         } else {
            ztrmm('Left', UPLO, 'No transpose', 'Unit', N-1, N, CONE, AFAC( 2, 1 ), LDAFAC, C( 2, 1 ), LDC );
         }

         // Call ZTRMM again to multiply by U (or L ).

         if ( LSAME( UPLO, 'U' ) ) {
            ztrmm('Right', UPLO, 'No transpose', 'Unit', N, N-1, CONE, AFAC( 1, 2 ), LDAFAC, C( 1, 2 ), LDC );
         } else {
            ztrmm('Right', UPLO, 'Conjugate transpose', 'Unit', N, N-1, CONE, AFAC( 2, 1 ), LDAFAC, C( 1, 2 ), LDC );
         }

         // Apply hermitian pivots

         DO J = N, 1, -1;
            I = IPIV( J );
            if (I != J) CALL ZSWAP( N, C( J, 1 ), LDC, C( I, 1 ), LDC );
         }
         DO J = N, 1, -1;
            I = IPIV( J );
            if (I != J) CALL ZSWAP( N, C( 1, J ), 1, C( 1, I ), 1 );
         }
      }


      // Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               C( I, J ) = C( I, J ) - A( I, J );
            }
         }
      } else {
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               C( I, J ) = C( I, J ) - A( I, J );
            }
         }
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = ZLANHE( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS;
      }

      return;

      // End of ZHET01_AA

      }
