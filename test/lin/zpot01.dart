      void zpot01(UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      Complex         A( LDA, * ), AFAC( LDAFAC, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      double             ANORM, EPS, TR;
      Complex         TC;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- double             DLAMCH, ZLANHE;
      //- Complex         ZDOTC;
      // EXTERNAL LSAME, DLAMCH, ZLANHE, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHER, ZSCAL, ZTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' );
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) { // 10
         if ( DIMAG( AFAC( J, J ) ) != ZERO ) {
            RESID = ONE / EPS;
            return;
         }
      } // 10

      // Compute the product U**H * U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {
         for (K = N; K >= 1; K--) { // 20

            // Compute the (K,K) element of the result.

            TR = DBLE( ZDOTC( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 ) );
            AFAC[K, K] = TR;

            // Compute the rest of column K.

            ztrmv('Upper', 'Conjugate', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );

         } // 20

      // Compute the product L * L**H, overwriting L.

      } else {
         for (K = N; K >= 1; K--) { // 30

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (K+1 <= N) zher( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC );

            // Scale column K by the diagonal element.

            TC = AFAC( K, K );
            zscal(N-K+1, TC, AFAC( K, K ), 1 );

         } // 30
      }

      // Compute the difference L * L**H - A (or U**H * U - A).

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 50
            for (I = 1; I <= J - 1; I++) { // 40
               AFAC[I, J] = AFAC( I, J ) - A( I, J );
            } // 40
            AFAC[J, J] = AFAC( J, J ) - DBLE( A( J, J ) );
         } // 50
      } else {
         for (J = 1; J <= N; J++) { // 70
            AFAC[J, J] = AFAC( J, J ) - DBLE( A( J, J ) );
            for (I = J + 1; I <= N; I++) { // 60
               AFAC[I, J] = AFAC( I, J ) - A( I, J );
            } // 60
         } // 70
      }

      // Compute norm(L*U - A) / ( N * norm(A) * EPS )

      RESID = ZLANHE( '1', UPLO, N, AFAC, LDAFAC, RWORK );

      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS;

      return;
      }
