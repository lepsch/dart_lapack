      SUBROUTINE ZPPT01( UPLO, N, A, AFAC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( * ), AFAC( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, K, KC;
      double             ANORM, EPS, TR;
      COMPLEX*16         TC
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHP;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, DLAMCH, ZLANHP, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHPR, ZSCAL, ZTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHP( '1', UPLO, N, A, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      KC = 1
      if ( LSAME( UPLO, 'U' ) ) {
         DO 10 K = 1, N
            if ( DIMAG( AFAC( KC ) ).NE.ZERO ) {
               RESID = ONE / EPS
               RETURN
            }
            KC = KC + K + 1
   10    CONTINUE
      } else {
         DO 20 K = 1, N
            if ( DIMAG( AFAC( KC ) ).NE.ZERO ) {
               RESID = ONE / EPS
               RETURN
            }
            KC = KC + N - K + 1
   20    CONTINUE
      }

      // Compute the product U'*U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {
         KC = ( N*( N-1 ) ) / 2 + 1
         DO 30 K = N, 1, -1

            // Compute the (K,K) element of the result.

            TR = DBLE( ZDOTC( K, AFAC( KC ), 1, AFAC( KC ), 1 ) )
            AFAC( KC+K-1 ) = TR

            // Compute the rest of column K.

            if ( K.GT.1 ) {
               CALL ZTPMV( 'Upper', 'Conjugate', 'Non-unit', K-1, AFAC, AFAC( KC ), 1 )
               KC = KC - ( K-1 )
            }
   30    CONTINUE

         // Compute the difference  L*L' - A

         KC = 1
         DO 50 K = 1, N
            DO 40 I = 1, K - 1
               AFAC( KC+I-1 ) = AFAC( KC+I-1 ) - A( KC+I-1 )
   40       CONTINUE
            AFAC( KC+K-1 ) = AFAC( KC+K-1 ) - DBLE( A( KC+K-1 ) )
            KC = KC + K
   50    CONTINUE

      // Compute the product L*L', overwriting L.

      } else {
         KC = ( N*( N+1 ) ) / 2
         DO 60 K = N, 1, -1

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            IF( K.LT.N ) CALL ZHPR( 'Lower', N-K, ONE, AFAC( KC+1 ), 1, AFAC( KC+N-K+1 ) )

            // Scale column K by the diagonal element.

            TC = AFAC( KC )
            CALL ZSCAL( N-K+1, TC, AFAC( KC ), 1 )

            KC = KC - ( N-K+2 )
   60    CONTINUE

         // Compute the difference  U'*U - A

         KC = 1
         DO 80 K = 1, N
            AFAC( KC ) = AFAC( KC ) - DBLE( A( KC ) )
            DO 70 I = K + 1, N
               AFAC( KC+I-K ) = AFAC( KC+I-K ) - A( KC+I-K )
   70       CONTINUE
            KC = KC + N - K + 1
   80    CONTINUE
      }

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = ZLANHP( '1', UPLO, N, AFAC, RWORK )

      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS

      RETURN

      // End of ZPPT01

      }
