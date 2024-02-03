      SUBROUTINE DPST01( UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM, PIV, RWORK, RESID, RANK );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             RESID;
      int                LDA, LDAFAC, LDPERM, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AFAC( LDAFAC, * ), PERM( LDPERM, * ), RWORK( * );
      int                PIV( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double             ANORM, EPS, T;
      int                I, J, K;
      // ..
      // .. External Functions ..
      double             DDOT, DLAMCH, DLANSY;
      bool               LSAME;
      // EXTERNAL DDOT, DLAMCH, DLANSY, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSYR, DTRMV
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

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' );
      ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the product U'*U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {

         if ( RANK < N ) {
            for (J = RANK + 1; J <= N; J++) { // 110
               for (I = RANK + 1; I <= J; I++) { // 100
                  AFAC( I, J ) = ZERO;
               } // 100
            } // 110
         }

         DO 120 K = N, 1, -1;

            // Compute the (K,K) element of the result.

            T = DDOT( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 );
            AFAC( K, K ) = T;

            // Compute the rest of column K.

            dtrmv('Upper', 'Transpose', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );

         } // 120

      // Compute the product L*L', overwriting L.

      } else {

         if ( RANK < N ) {
            for (J = RANK + 1; J <= N; J++) { // 140
               for (I = J; I <= N; I++) { // 130
                  AFAC( I, J ) = ZERO;
               } // 130
            } // 140
         }

         DO 150 K = N, 1, -1;
            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (K+1 <= N) CALL DSYR( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC );

            // Scale column K by the diagonal element.

            T = AFAC( K, K );
            dscal(N-K+1, T, AFAC( K, K ), 1 );
         } // 150

      }

         // Form P*L*L'*P' or P*U'*U*P'

      if ( LSAME( UPLO, 'U' ) ) {

         for (J = 1; J <= N; J++) { // 170
            for (I = 1; I <= N; I++) { // 160
               if ( PIV( I ) <= PIV( J ) ) {
                  if ( I <= J ) {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J );
                  } else {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( J, I );
                  }
               }
            } // 160
         } // 170


      } else {

         for (J = 1; J <= N; J++) { // 190
            for (I = 1; I <= N; I++) { // 180
               if ( PIV( I ) >= PIV( J ) ) {
                  if ( I >= J ) {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J );
                  } else {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( J, I );
                  }
               }
            } // 180
         } // 190

      }

      // Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 210
            for (I = 1; I <= J; I++) { // 200
               PERM( I, J ) = PERM( I, J ) - A( I, J );
            } // 200
         } // 210
      } else {
         for (J = 1; J <= N; J++) { // 230
            for (I = J; I <= N; I++) { // 220
               PERM( I, J ) = PERM( I, J ) - A( I, J );
            } // 220
         } // 230
      }

      // Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
      // ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).

      RESID = DLANSY( '1', UPLO, N, PERM, LDAFAC, RWORK );

      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS;

      return;

      // End of DPST01

      }
