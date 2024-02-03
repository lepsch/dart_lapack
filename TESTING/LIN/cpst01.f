      SUBROUTINE CPST01( UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM, PIV, RWORK, RESID, RANK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               RESID
      int                LDA, LDAFAC, LDPERM, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), PERM( LDPERM, * )
      REAL               RWORK( * )
      int                PIV( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO
      const              CZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      COMPLEX            TC
      REAL               ANORM, EPS, TR
      int                I, J, K;
      // ..
      // .. External Functions ..
      COMPLEX            CDOTC
      REAL               CLANHE, SLAMCH
      bool               LSAME;
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

      DO 100 J = 1, N
         if ( AIMAG( AFAC( J, J ) ).NE.ZERO ) {
            RESID = ONE / EPS
            RETURN
         }
  100 CONTINUE

      // Compute the product U'*U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {

         if ( RANK.LT.N ) {
            DO 120 J = RANK + 1, N
               DO 110 I = RANK + 1, J
                  AFAC( I, J ) = CZERO
  110          CONTINUE
  120       CONTINUE
         }

         DO 130 K = N, 1, -1

            // Compute the (K,K) element of the result.

            TR = REAL( CDOTC( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 ) )
            AFAC( K, K ) = TR

            // Compute the rest of column K.

            ctrmv('Upper', 'Conjugate', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );

  130    CONTINUE

      // Compute the product L*L', overwriting L.

      } else {

         if ( RANK.LT.N ) {
            DO 150 J = RANK + 1, N
               DO 140 I = J, N
                  AFAC( I, J ) = CZERO
  140          CONTINUE
  150       CONTINUE
         }

         DO 160 K = N, 1, -1
            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            IF( K+1.LE.N ) CALL CHER( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC )

            // Scale column K by the diagonal element.

            TC = AFAC( K, K )
            cscal(N-K+1, TC, AFAC( K, K ), 1 );
  160    CONTINUE

      }

         // Form P*L*L'*P' or P*U'*U*P'

      if ( LSAME( UPLO, 'U' ) ) {

         DO 180 J = 1, N
            DO 170 I = 1, N
               if ( PIV( I ).LE.PIV( J ) ) {
                  if ( I.LE.J ) {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J )
                  } else {
                     PERM( PIV( I ), PIV( J ) ) = CONJG( AFAC( J, I ) )
                  }
               }
  170       CONTINUE
  180    CONTINUE


      } else {

         DO 200 J = 1, N
            DO 190 I = 1, N
               if ( PIV( I ).GE.PIV( J ) ) {
                  if ( I.GE.J ) {
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J )
                  } else {
                     PERM( PIV( I ), PIV( J ) ) = CONJG( AFAC( J, I ) )
                  }
               }
  190       CONTINUE
  200    CONTINUE

      }

      // Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).

      if ( LSAME( UPLO, 'U' ) ) {
         DO 220 J = 1, N
            DO 210 I = 1, J - 1
               PERM( I, J ) = PERM( I, J ) - A( I, J )
  210       CONTINUE
            PERM( J, J ) = PERM( J, J ) - REAL( A( J, J ) )
  220    CONTINUE
      } else {
         DO 240 J = 1, N
            PERM( J, J ) = PERM( J, J ) - REAL( A( J, J ) )
            DO 230 I = J + 1, N
               PERM( I, J ) = PERM( I, J ) - A( I, J )
  230       CONTINUE
  240    CONTINUE
      }

      // Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
      // ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).

      RESID = CLANHE( '1', UPLO, N, PERM, LDAFAC, RWORK )

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS

      RETURN

      // End of CPST01

      }
