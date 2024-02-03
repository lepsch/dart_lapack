      void cpst01(UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM, PIV, RWORK, RESID, RANK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               RESID;
      int                LDA, LDAFAC, LDPERM, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), PERM( LDPERM, * );
      REAL               RWORK( * );
      int                PIV( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      COMPLEX            TC;
      REAL               ANORM, EPS, TR;
      int                I, J, K;
      // ..
      // .. External Functions ..
      //- COMPLEX            CDOTC;
      //- REAL               CLANHE, SLAMCH;
      //- bool               LSAME;
      // EXTERNAL CDOTC, CLANHE, SLAMCH, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHER, CSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CONJG, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) { // 100
         if ( AIMAG( AFAC( J, J ) ) != ZERO ) {
            RESID = ONE / EPS;
            return;
         }
      } // 100

      // Compute the product U'*U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {

         if ( RANK < N ) {
            for (J = RANK + 1; J <= N; J++) { // 120
               for (I = RANK + 1; I <= J; I++) { // 110
                  AFAC( I, J ) = CZERO;
               } // 110
            } // 120
         }

         for (K = N; K >= 1; K--) { // 130

            // Compute the (K,K) element of the result.

            TR = REAL( CDOTC( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 ) );
            AFAC( K, K ) = TR;

            // Compute the rest of column K.

            ctrmv('Upper', 'Conjugate', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );

         } // 130

      // Compute the product L*L', overwriting L.

      } else {

         if ( RANK < N ) {
            for (J = RANK + 1; J <= N; J++) { // 150
               for (I = J; I <= N; I++) { // 140
                  AFAC( I, J ) = CZERO;
               } // 140
            } // 150
         }

         for (K = N; K >= 1; K--) { // 160
            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (K+1 <= N) cher( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC );

            // Scale column K by the diagonal element.

            TC = AFAC( K, K );
            cscal(N-K+1, TC, AFAC( K, K ), 1 );
         } // 160

      }

         // Form P*L*L'*P' or P*U'*U*P'

      if ( LSAME( UPLO, 'U' ) ) {

         for (J = 1; J <= N; J++) { // 180
            for (I = 1; I <= N; I++) { // 170
               if ( PIV( I ) <= PIV( J ) ) {
                  if ( I <= J ) {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J );
                  } else {
                     PERM( PIV( I ), PIV( J ) ) = CONJG( AFAC( J, I ) );
                  }
               }
            } // 170
         } // 180


      } else {

         for (J = 1; J <= N; J++) { // 200
            for (I = 1; I <= N; I++) { // 190
               if ( PIV( I ) >= PIV( J ) ) {
                  if ( I >= J ) {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J );
                  } else {
                     PERM( PIV( I ), PIV( J ) ) = CONJG( AFAC( J, I ) );
                  }
               }
            } // 190
         } // 200

      }

      // Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 220
            for (I = 1; I <= J - 1; I++) { // 210
               PERM( I, J ) = PERM( I, J ) - A( I, J );
            } // 210
            PERM( J, J ) = PERM( J, J ) - REAL( A( J, J ) );
         } // 220
      } else {
         for (J = 1; J <= N; J++) { // 240
            PERM( J, J ) = PERM( J, J ) - REAL( A( J, J ) );
            for (I = J + 1; I <= N; I++) { // 230
               PERM( I, J ) = PERM( I, J ) - A( I, J );
            } // 230
         } // 240
      }

      // Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
      // ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).

      RESID = CLANHE( '1', UPLO, N, PERM, LDAFAC, RWORK );

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;

      return;
      }
