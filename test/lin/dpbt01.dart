      void dpbt01(UPLO, N, KD, final Matrix<double> A, final int LDA, final Matrix<double> AFAC, final int LDAFAC, final Array<double> RWORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                KD, LDA, LDAFAC, N;
      double             RESID;
      double             A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * );
      // ..

// =====================================================================


      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, J, K, KC, KLEN, ML, MU;
      double             ANORM, EPS, T;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DDOT, DLAMCH, DLANSB;
      // EXTERNAL lsame, DDOT, DLAMCH, DLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSYR, DTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = dlamch( 'Epsilon' );
      ANORM = dlansb( '1', UPLO, N, KD, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the product U'*U, overwriting U.

      if ( lsame( UPLO, 'U' ) ) {
         for (K = N; K >= 1; K--) { // 10
            KC = max( 1, KD+2-K );
            KLEN = KD + 1 - KC;

            // Compute the (K,K) element of the result.

            T = ddot( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 );
            AFAC[KD+1][K] = T;

            // Compute the rest of column K.

            if (KLEN > 0) dtrmv( 'Upper', 'Transpose', 'Non-unit', KLEN, AFAC( KD+1, K-KLEN ), LDAFAC-1, AFAC( KC, K ), 1 );

         } // 10

      // UPLO = 'L':  Compute the product L*L', overwriting L.

      } else {
         for (K = N; K >= 1; K--) { // 20
            KLEN = min( KD, N-K );

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (KLEN > 0) dsyr( 'Lower', KLEN, ONE, AFAC( 2, K ), 1, AFAC( 1, K+1 ), LDAFAC-1 );

            // Scale column K by the diagonal element.

            T = AFAC( 1, K );
            dscal(KLEN+1, T, AFAC( 1, K ), 1 );

         } // 20
      }

      // Compute the difference  L*L' - A  or  U'*U - A.

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 40
            MU = max( 1, KD+2-J );
            for (I = MU; I <= KD + 1; I++) { // 30
               AFAC[I][J] = AFAC( I, J ) - A( I, J );
            } // 30
         } // 40
      } else {
         for (J = 1; J <= N; J++) { // 60
            ML = min( KD+1, N-J+1 );
            for (I = 1; I <= ML; I++) { // 50
               AFAC[I][J] = AFAC( I, J ) - A( I, J );
            } // 50
         } // 60
      }

      // Compute norm( L*L' - A ) / ( N * norm(A) * EPS )

      RESID = dlansb( 'I', UPLO, N, KD, AFAC, LDAFAC, RWORK );

      RESID = ( ( RESID / N.toDouble() ) / ANORM ) / EPS;

      }
