      SUBROUTINE SPST01( UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM, PIV, RWORK, RESID, RANK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               RESID
      int                LDA, LDAFAC, LDPERM, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), AFAC( LDAFAC, * ), PERM( LDPERM, * ), RWORK( * )
      int                PIV( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               ANORM, EPS, T
      int                I, J, K;
      // ..
      // .. External Functions ..
      REAL               SDOT, SLAMCH, SLANSY
      bool               LSAME;
      // EXTERNAL SDOT, SLAMCH, SLANSY, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSYR, STRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute the product U'*U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {

         if ( RANK.LT.N ) {
            DO 110 J = RANK + 1, N
               DO 100 I = RANK + 1, J
                  AFAC( I, J ) = ZERO
  100          CONTINUE
  110       CONTINUE
         }

         DO 120 K = N, 1, -1

            // Compute the (K,K) element of the result.

            T = SDOT( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 )
            AFAC( K, K ) = T

            // Compute the rest of column K.

            strmv('Upper', 'Transpose', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );

  120    CONTINUE

      // Compute the product L*L', overwriting L.

      } else {

         if ( RANK.LT.N ) {
            DO 140 J = RANK + 1, N
               for (I = J; I <= N; I++) { // 130
                  AFAC( I, J ) = ZERO
  130          CONTINUE
  140       CONTINUE
         }

         DO 150 K = N, 1, -1
            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            IF( K+1.LE.N ) CALL SSYR( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC )

            // Scale column K by the diagonal element.

            T = AFAC( K, K )
            sscal(N-K+1, T, AFAC( K, K ), 1 );
  150    CONTINUE

      }

         // Form P*L*L'*P' or P*U'*U*P'

      if ( LSAME( UPLO, 'U' ) ) {

         for (J = 1; J <= N; J++) { // 170
            for (I = 1; I <= N; I++) { // 160
               if ( PIV( I ).LE.PIV( J ) ) {
                  if ( I.LE.J ) {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J )
                  } else {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( J, I )
                  }
               }
  160       CONTINUE
  170    CONTINUE


      } else {

         for (J = 1; J <= N; J++) { // 190
            for (I = 1; I <= N; I++) { // 180
               if ( PIV( I ).GE.PIV( J ) ) {
                  if ( I.GE.J ) {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J )
                  } else {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( J, I )
                  }
               }
  180       CONTINUE
  190    CONTINUE

      }

      // Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 210
            for (I = 1; I <= J; I++) { // 200
               PERM( I, J ) = PERM( I, J ) - A( I, J )
  200       CONTINUE
  210    CONTINUE
      } else {
         for (J = 1; J <= N; J++) { // 230
            for (I = J; I <= N; I++) { // 220
               PERM( I, J ) = PERM( I, J ) - A( I, J )
  220       CONTINUE
  230    CONTINUE
      }

      // Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
      // ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).

      RESID = SLANSY( '1', UPLO, N, PERM, LDAFAC, RWORK )

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS

      RETURN

      // End of SPST01

      }
