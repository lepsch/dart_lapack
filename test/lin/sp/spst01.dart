      void spst01(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> AFAC_, final int LDAFAC, final Matrix<double> PERM_, final int LDPERM, final int PIV, final Array<double> RWORK_, final int RESID, final int RANK,) {
  final A = A_.dim();
  final AFAC = AFAC_.dim();
  final PERM = PERM_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               RESID;
      int                LDA, LDAFAC, LDPERM, N, RANK;
      String             UPLO;
      double               A( LDA, * ), AFAC( LDAFAC, * ), PERM( LDPERM, * ), RWORK( * );
      int                PIV( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ANORM, EPS, T;
      int                I, J, K;
      // ..
      // .. External Functions ..
      //- REAL               SDOT, SLAMCH, SLANSY;
      //- bool               lsame;
      // EXTERNAL SDOT, SLAMCH, SLANSY, lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSYR, STRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the product U'*U, overwriting U.

      if ( lsame( UPLO, 'U' ) ) {

         if ( RANK < N ) {
            for (J = RANK + 1; J <= N; J++) { // 110
               for (I = RANK + 1; I <= J; I++) { // 100
                  AFAC[I][J] = ZERO;
               } // 100
            } // 110
         }

         for (K = N; K >= 1; K--) { // 120

            // Compute the (K,K) element of the result.

            T = SDOT( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 );
            AFAC[K][K] = T;

            // Compute the rest of column K.

            strmv('Upper', 'Transpose', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );

         } // 120

      // Compute the product L*L', overwriting L.

      } else {

         if ( RANK < N ) {
            for (J = RANK + 1; J <= N; J++) { // 140
               for (I = J; I <= N; I++) { // 130
                  AFAC[I][J] = ZERO;
               } // 130
            } // 140
         }

         for (K = N; K >= 1; K--) { // 150
            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (K+1 <= N) ssyr( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC );

            // Scale column K by the diagonal element.

            T = AFAC( K, K );
            sscal(N-K+1, T, AFAC( K, K ), 1 );
         } // 150

      }

         // Form P*L*L'*P' or P*U'*U*P'

      if ( lsame( UPLO, 'U' ) ) {

         for (J = 1; J <= N; J++) { // 170
            for (I = 1; I <= N; I++) { // 160
               if ( PIV( I ) <= PIV( J ) ) {
                  if ( I <= J ) {
                     PERM[PIV( I )][PIV( J )] = AFAC( I, J );
                  } else {
                     PERM[PIV( I )][PIV( J )] = AFAC( J, I );
                  }
               }
            } // 160
         } // 170


      } else {

         for (J = 1; J <= N; J++) { // 190
            for (I = 1; I <= N; I++) { // 180
               if ( PIV( I ) >= PIV( J ) ) {
                  if ( I >= J ) {
                     PERM[PIV( I )][PIV( J )] = AFAC( I, J );
                  } else {
                     PERM[PIV( I )][PIV( J )] = AFAC( J, I );
                  }
               }
            } // 180
         } // 190

      }

      // Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 210
            for (I = 1; I <= J; I++) { // 200
               PERM[I][J] = PERM( I, J ) - A( I, J );
            } // 200
         } // 210
      } else {
         for (J = 1; J <= N; J++) { // 230
            for (I = J; I <= N; I++) { // 220
               PERM[I][J] = PERM( I, J ) - A( I, J );
            } // 220
         } // 230
      }

      // Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
      // ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).

      RESID = SLANSY( '1', UPLO, N, PERM, LDAFAC, RWORK );

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;

      }