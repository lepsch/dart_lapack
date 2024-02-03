      SUBROUTINE DPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDA, LDAFAC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * );
      // ..

*  =====================================================================


      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
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

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANSB( '1', UPLO, N, KD, A, LDA, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute the product U'*U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 10 K = N, 1, -1
            KC = MAX( 1, KD+2-K )
            KLEN = KD + 1 - KC

            // Compute the (K,K) element of the result.

            T = DDOT( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 )
            AFAC( KD+1, K ) = T

            // Compute the rest of column K.

            IF( KLEN.GT.0 ) CALL DTRMV( 'Upper', 'Transpose', 'Non-unit', KLEN, AFAC( KD+1, K-KLEN ), LDAFAC-1, AFAC( KC, K ), 1 )

   10    CONTINUE

      // UPLO = 'L':  Compute the product L*L', overwriting L.

      } else {
         DO 20 K = N, 1, -1
            KLEN = MIN( KD, N-K )

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            IF( KLEN.GT.0 ) CALL DSYR( 'Lower', KLEN, ONE, AFAC( 2, K ), 1, AFAC( 1, K+1 ), LDAFAC-1 )

            // Scale column K by the diagonal element.

            T = AFAC( 1, K )
            dscal(KLEN+1, T, AFAC( 1, K ), 1 );

   20    CONTINUE
      }

      // Compute the difference  L*L' - A  or  U'*U - A.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 40 J = 1, N
            MU = MAX( 1, KD+2-J )
            DO 30 I = MU, KD + 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   30       CONTINUE
   40    CONTINUE
      } else {
         DO 60 J = 1, N
            ML = MIN( KD+1, N-J+1 )
            DO 50 I = 1, ML
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   50       CONTINUE
   60    CONTINUE
      }

      // Compute norm( L*L' - A ) / ( N * norm(A) * EPS )

      RESID = DLANSB( 'I', UPLO, N, KD, AFAC, LDAFAC, RWORK )

      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS

      RETURN

      // End of DPBT01

      }
