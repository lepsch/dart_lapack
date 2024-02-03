      SUBROUTINE DPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK, RESID );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDA, LDAFAC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * );
      // ..

// =====================================================================


      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K, KC, KLEN, ML, MU;
      double             ANORM, EPS, T;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DDOT, DLAMCH, DLANSB;
      // EXTERNAL LSAME, DDOT, DLAMCH, DLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSYR, DTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' );
      ANORM = DLANSB( '1', UPLO, N, KD, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the product U'*U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 10 K = N, 1, -1;
            KC = MAX( 1, KD+2-K );
            KLEN = KD + 1 - KC;

            // Compute the (K,K) element of the result.

            T = DDOT( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 );
            AFAC( KD+1, K ) = T;

            // Compute the rest of column K.

            if (KLEN > 0) CALL DTRMV( 'Upper', 'Transpose', 'Non-unit', KLEN, AFAC( KD+1, K-KLEN ), LDAFAC-1, AFAC( KC, K ), 1 );

         } // 10

      // UPLO = 'L':  Compute the product L*L', overwriting L.

      } else {
         DO 20 K = N, 1, -1;
            KLEN = MIN( KD, N-K );

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (KLEN > 0) CALL DSYR( 'Lower', KLEN, ONE, AFAC( 2, K ), 1, AFAC( 1, K+1 ), LDAFAC-1 );

            // Scale column K by the diagonal element.

            T = AFAC( 1, K );
            dscal(KLEN+1, T, AFAC( 1, K ), 1 );

         } // 20
      }

      // Compute the difference  L*L' - A  or  U'*U - A.

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 40
            MU = MAX( 1, KD+2-J );
            for (I = MU; I <= KD + 1; I++) { // 30
               AFAC( I, J ) = AFAC( I, J ) - A( I, J );
            } // 30
         } // 40
      } else {
         for (J = 1; J <= N; J++) { // 60
            ML = MIN( KD+1, N-J+1 );
            for (I = 1; I <= ML; I++) { // 50
               AFAC( I, J ) = AFAC( I, J ) - A( I, J );
            } // 50
         } // 60
      }

      // Compute norm( L*L' - A ) / ( N * norm(A) * EPS )

      RESID = DLANSB( 'I', UPLO, N, KD, AFAC, LDAFAC, RWORK );

      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS;

      return;

      // End of DPBT01

      }
