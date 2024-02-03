      SUBROUTINE ZPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDA, LDAFAC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * )
      // ..

*  =====================================================================


      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K, KC, KLEN, ML, MU;
      double             AKK, ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHB;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, DLAMCH, ZLANHB, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZHER, ZTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHB( '1', UPLO, N, KD, A, LDA, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 10 J = 1, N
            if ( DIMAG( AFAC( KD+1, J ) ).NE.ZERO ) {
               RESID = ONE / EPS
               RETURN
            }
   10    CONTINUE
      } else {
         DO 20 J = 1, N
            if ( DIMAG( AFAC( 1, J ) ).NE.ZERO ) {
               RESID = ONE / EPS
               RETURN
            }
   20    CONTINUE
      }

      // Compute the product U'*U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 30 K = N, 1, -1
            KC = MAX( 1, KD+2-K )
            KLEN = KD + 1 - KC

            // Compute the (K,K) element of the result.

            AKK = DBLE( ZDOTC( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 ) )
            AFAC( KD+1, K ) = AKK

            // Compute the rest of column K.

            IF( KLEN.GT.0 ) CALL ZTRMV( 'Upper', 'Conjugate', 'Non-unit', KLEN, AFAC( KD+1, K-KLEN ), LDAFAC-1, AFAC( KC, K ), 1 )

   30    CONTINUE

      // UPLO = 'L':  Compute the product L*L', overwriting L.

      } else {
         DO 40 K = N, 1, -1
            KLEN = MIN( KD, N-K )

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            IF( KLEN.GT.0 ) CALL ZHER( 'Lower', KLEN, ONE, AFAC( 2, K ), 1, AFAC( 1, K+1 ), LDAFAC-1 )

            // Scale column K by the diagonal element.

            AKK = DBLE( AFAC( 1, K ) )
            zdscal(KLEN+1, AKK, AFAC( 1, K ), 1 );

   40    CONTINUE
      }

      // Compute the difference  L*L' - A  or  U'*U - A.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 60 J = 1, N
            MU = MAX( 1, KD+2-J )
            DO 50 I = MU, KD + 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   50       CONTINUE
   60    CONTINUE
      } else {
         DO 80 J = 1, N
            ML = MIN( KD+1, N-J+1 )
            DO 70 I = 1, ML
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   70       CONTINUE
   80    CONTINUE
      }

      // Compute norm( L*L' - A ) / ( N * norm(A) * EPS )

      RESID = ZLANHB( '1', UPLO, N, KD, AFAC, LDAFAC, RWORK )

      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS

      RETURN

      // End of ZPBT01

      }
