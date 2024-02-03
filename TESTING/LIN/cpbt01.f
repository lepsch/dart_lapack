      SUBROUTINE CPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDA, LDAFAC, N;
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
      int                I, J, K, KC, KLEN, ML, MU;
      REAL               AKK, ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHB, SLAMCH
      COMPLEX            CDOTC
      // EXTERNAL LSAME, CLANHB, SLAMCH, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHER, CSSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANHB( '1', UPLO, N, KD, A, LDA, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 10
            if ( AIMAG( AFAC( KD+1, J ) ).NE.ZERO ) {
               RESID = ONE / EPS
               RETURN
            }
   10    CONTINUE
      } else {
         for (J = 1; J <= N; J++) { // 20
            if ( AIMAG( AFAC( 1, J ) ).NE.ZERO ) {
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

            AKK = REAL( CDOTC( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 ) )
            AFAC( KD+1, K ) = AKK

            // Compute the rest of column K.

            IF( KLEN.GT.0 ) CALL CTRMV( 'Upper', 'Conjugate', 'Non-unit', KLEN, AFAC( KD+1, K-KLEN ), LDAFAC-1, AFAC( KC, K ), 1 )

   30    CONTINUE

      // UPLO = 'L':  Compute the product L*L', overwriting L.

      } else {
         DO 40 K = N, 1, -1
            KLEN = MIN( KD, N-K )

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            IF( KLEN.GT.0 ) CALL CHER( 'Lower', KLEN, ONE, AFAC( 2, K ), 1, AFAC( 1, K+1 ), LDAFAC-1 )

            // Scale column K by the diagonal element.

            AKK = REAL( AFAC( 1, K ) )
            csscal(KLEN+1, AKK, AFAC( 1, K ), 1 );

   40    CONTINUE
      }

      // Compute the difference  L*L' - A  or  U'*U - A.

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 60
            MU = MAX( 1, KD+2-J )
            DO 50 I = MU, KD + 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   50       CONTINUE
   60    CONTINUE
      } else {
         for (J = 1; J <= N; J++) { // 80
            ML = MIN( KD+1, N-J+1 )
            for (I = 1; I <= ML; I++) { // 70
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   70       CONTINUE
   80    CONTINUE
      }

      // Compute norm( L*L' - A ) / ( N * norm(A) * EPS )

      RESID = CLANHB( '1', UPLO, N, KD, AFAC, LDAFAC, RWORK )

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS

      RETURN

      // End of CPBT01

      }
