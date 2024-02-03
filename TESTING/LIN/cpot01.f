      SUBROUTINE CPOT01( UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      REAL               ANORM, EPS, TR
      COMPLEX            TC
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHE, SLAMCH
      COMPLEX            CDOTC
      // EXTERNAL LSAME, CLANHE, SLAMCH, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHER, CSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) { // 10
         if ( AIMAG( AFAC( J, J ) ).NE.ZERO ) {
            RESID = ONE / EPS
            RETURN
         }
      } // 10

      // Compute the product U**H * U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 20 K = N, 1, -1

            // Compute the (K,K) element of the result.

            TR = REAL( CDOTC( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 ) )
            AFAC( K, K ) = TR

            // Compute the rest of column K.

            ctrmv('Upper', 'Conjugate', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );

         } // 20

      // Compute the product L * L**H, overwriting L.

      } else {
         DO 30 K = N, 1, -1

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            IF( K+1.LE.N ) CALL CHER( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC )

            // Scale column K by the diagonal element.

            TC = AFAC( K, K )
            cscal(N-K+1, TC, AFAC( K, K ), 1 );

         } // 30
      }

      // Compute the difference L * L**H - A (or U**H * U - A).

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 50
            for (I = 1; I <= J - 1; I++) { // 40
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
            } // 40
            AFAC( J, J ) = AFAC( J, J ) - REAL( A( J, J ) )
         } // 50
      } else {
         for (J = 1; J <= N; J++) { // 70
            AFAC( J, J ) = AFAC( J, J ) - REAL( A( J, J ) )
            for (I = J + 1; I <= N; I++) { // 60
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
            } // 60
         } // 70
      }

      // Compute norm(L*U - A) / ( N * norm(A) * EPS )

      RESID = CLANHE( '1', UPLO, N, AFAC, LDAFAC, RWORK )

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS

      RETURN

      // End of CPOT01

      }
