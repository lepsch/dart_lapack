      void cpbt01(final int UPLO, final int N, final int KD, final Matrix<double> A_, final int LDA, final Matrix<double> AFAC_, final int LDAFAC, final Array<double> RWORK_, final int RESID,) {
  final A = A_.dim();
  final AFAC = AFAC_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                KD, LDA, LDAFAC, N;
      double               RESID;
      double               RWORK( * );
      Complex            A( LDA, * ), AFAC( LDAFAC, * );
      // ..

// =====================================================================


      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, J, K, KC, KLEN, ML, MU;
      double               AKK, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANHB, SLAMCH;
      //- COMPLEX            CDOTC;
      // EXTERNAL lsame, CLANHB, SLAMCH, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHER, CSSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, MAX, MIN, REAL

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANHB( '1', UPLO, N, KD, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 10
            if ( AIMAG( AFAC( KD+1, J ) ) != ZERO ) {
               RESID = ONE / EPS;
               return;
            }
         } // 10
      } else {
         for (J = 1; J <= N; J++) { // 20
            if ( AIMAG( AFAC( 1, J ) ) != ZERO ) {
               RESID = ONE / EPS;
               return;
            }
         } // 20
      }

      // Compute the product U'*U, overwriting U.

      if ( lsame( UPLO, 'U' ) ) {
         for (K = N; K >= 1; K--) { // 30
            KC = max( 1, KD+2-K );
            KLEN = KD + 1 - KC;

            // Compute the (K,K) element of the result.

            AKK = double( CDOTC( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 ) );
            AFAC[KD+1][K] = AKK;

            // Compute the rest of column K.

            if (KLEN > 0) ctrmv( 'Upper', 'Conjugate', 'Non-unit', KLEN, AFAC( KD+1, K-KLEN ), LDAFAC-1, AFAC( KC, K ), 1 );

         } // 30

      // UPLO = 'L':  Compute the product L*L', overwriting L.

      } else {
         for (K = N; K >= 1; K--) { // 40
            KLEN = min( KD, N-K );

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (KLEN > 0) cher( 'Lower', KLEN, ONE, AFAC( 2, K ), 1, AFAC( 1, K+1 ), LDAFAC-1 );

            // Scale column K by the diagonal element.

            AKK = double( AFAC( 1, K ) );
            csscal(KLEN+1, AKK, AFAC( 1, K ), 1 );

         } // 40
      }

      // Compute the difference  L*L' - A  or  U'*U - A.

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 60
            MU = max( 1, KD+2-J );
            for (I = MU; I <= KD + 1; I++) { // 50
               AFAC[I][J] = AFAC( I, J ) - A( I, J );
            } // 50
         } // 60
      } else {
         for (J = 1; J <= N; J++) { // 80
            ML = min( KD+1, N-J+1 );
            for (I = 1; I <= ML; I++) { // 70
               AFAC[I][J] = AFAC( I, J ) - A( I, J );
            } // 70
         } // 80
      }

      // Compute norm( L*L' - A ) / ( N * norm(A) * EPS )

      RESID = CLANHB( '1', UPLO, N, KD, AFAC, LDAFAC, RWORK );

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;

      }